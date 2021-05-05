
#include "wgriddingmsgridder.h"

#include "wgriddinggridder_simple.h"

#include "../msproviders/msreaders/msreader.h"

#include "../io/logger.h"

#include "../math/fftresampler.h"

#include "../msproviders/msprovider.h"

#include "../system/buffered_lane.h"

#include "../structures/imageweights.h"
#include "../structures/image.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

WGriddingMSGridder::WGriddingMSGridder(const Settings& settings)
    : MSGridderBase(settings),
      _cpuCount(_settings.threadCount),
      _accuracy(_settings.wgridderAccuracy) {
  _memSize = getAvailableMemory(_settings.memFraction, _settings.absMemLimit);
  // It may happen that several FFTResamplers are created concurrently, so we
  // must make sure that the FFTW planner can deal with this.
  fftwf_make_planner_thread_safe();
}

size_t WGriddingMSGridder::calculateMaxNRowsInMemory(
    size_t channelCount) const {
  size_t constantMem, perVisMem;
  _gridder->memUsage(constantMem, perVisMem);
  if (int64_t(constantMem) >= _memSize) {
    constantMem = _memSize / 2;
    Logger::Warn << "Not enough memory available for doing the gridding:\n"
                    "swapping might occur!\n";
  }
  uint64_t memForBuffers = _memSize - constantMem;

  uint64_t memPerRow = (perVisMem + sizeof(std::complex<float>)) *
                           channelCount       // vis themselves
                       + sizeof(double) * 3;  // uvw
  size_t maxNRows = std::max(memForBuffers / memPerRow, uint64_t(100));
  if (maxNRows < 1000) {
    Logger::Warn << "Less than 1000 data rows fit in memory: this probably "
                    "means performance is going to be very poor!\n";
  }

  return maxNRows;
}

void WGriddingMSGridder::gridMeasurementSet(MSData& msData) {
  const MultiBandData selectedBands(msData.SelectedBand());

  aocommon::UVector<std::complex<float>> modelBuffer(
      selectedBands.MaxChannels());
  aocommon::UVector<float> weightBuffer(selectedBands.MaxChannels());
  aocommon::UVector<bool> isSelected(selectedBands.MaxChannels(), true);

  size_t totalNRows = 0;
  for (size_t dataDescId = 0; dataDescId != selectedBands.DataDescCount();
       ++dataDescId) {
    const BandData& band = selectedBands[dataDescId];
    aocommon::UVector<double> frequencies(band.ChannelCount());
    for (size_t i = 0; i != frequencies.size(); ++i)
      frequencies[i] = band.Channel(i).Frequency();

    size_t maxNRows = calculateMaxNRowsInMemory(band.ChannelCount());

    aocommon::UVector<std::complex<float>> visBuffer(maxNRows *
                                                     band.ChannelCount());
    aocommon::UVector<double> uvwBuffer(maxNRows * 3);

    std::unique_ptr<MSReader> msReader = msData.msProvider->MakeReader();
    aocommon::UVector<std::complex<float>> newItemData(band.ChannelCount());
    InversionRow newRowData;
    newRowData.data = newItemData.data();

    // Iterate over chunks until all data has been gridded
    while (msReader->CurrentRowAvailable()) {
      Logger::Debug << "Max " << maxNRows << " rows fit in memory.\n";
      Logger::Info << "Loading data in memory...\n";

      size_t nRows = 0;

      // Read / fill the chunk
      while (msReader->CurrentRowAvailable() && nRows < maxNRows) {
        size_t rowDataDescId;
        double uInMeters, vInMeters, wInMeters;
        msReader->ReadMeta(uInMeters, vInMeters, wInMeters, rowDataDescId);
        if (rowDataDescId == dataDescId) {
          newRowData.uvw[0] = uInMeters;
          newRowData.uvw[1] = vInMeters;
          newRowData.uvw[2] = wInMeters;
          newRowData.dataDescId = dataDescId;
          readAndWeightVisibilities<1>(*msReader, newRowData, band,
                                       weightBuffer.data(), modelBuffer.data(),
                                       isSelected.data());

          std::copy_n(newRowData.data, band.ChannelCount(),
                      &visBuffer[nRows * band.ChannelCount()]);
          std::copy_n(newRowData.uvw, 3, &uvwBuffer[nRows * 3]);

          ++nRows;
        }
        msReader->NextInputRow();
      }

      Logger::Info << "Gridding " << nRows << " rows...\n";
      _gridder->AddInversionData(nRows, band.ChannelCount(), uvwBuffer.data(),
                                 frequencies.data(), visBuffer.data());

      totalNRows += nRows;
    }  // end of chunk
  }    // finished all chunks

  msData.totalRowsProcessed += totalNRows;
}

void WGriddingMSGridder::predictMeasurementSet(MSData& msData) {
  msData.msProvider->ReopenRW();
  const MultiBandData selectedBands(msData.SelectedBand());

  size_t totalNRows = 0;
  for (size_t dataDescId = 0; dataDescId != selectedBands.DataDescCount();
       ++dataDescId) {
    const BandData& band = selectedBands[dataDescId];
    aocommon::UVector<double> frequencies(band.ChannelCount());
    for (size_t i = 0; i != frequencies.size(); ++i)
      frequencies[i] = band.Channel(i).Frequency();

    size_t maxNRows = calculateMaxNRowsInMemory(band.ChannelCount());

    aocommon::UVector<double> uvwBuffer(maxNRows * 3);
    // Iterate over chunks until all data has been gridded
    msData.msProvider->ResetWritePosition();
    std::unique_ptr<MSReader> msReader = msData.msProvider->MakeReader();
    while (msReader->CurrentRowAvailable()) {
      size_t nRows = 0;
      // Read / fill the chunk
      while (msReader->CurrentRowAvailable() && nRows < maxNRows) {
        size_t rowDataDescId;
        double uInMeters, vInMeters, wInMeters;
        msReader->ReadMeta(uInMeters, vInMeters, wInMeters, rowDataDescId);
        if (rowDataDescId == dataDescId) {
          uvwBuffer[nRows * 3] = uInMeters;
          uvwBuffer[nRows * 3 + 1] = vInMeters;
          uvwBuffer[nRows * 3 + 2] = wInMeters;
          ++nRows;
        }
        msReader->NextInputRow();
      }

      Logger::Info << "Predicting " << nRows << " rows...\n";
      aocommon::UVector<std::complex<float>> visBuffer(maxNRows *
                                                       band.ChannelCount());
      _gridder->PredictVisibilities(nRows, band.ChannelCount(),
                                    uvwBuffer.data(), frequencies.data(),
                                    visBuffer.data());

      Logger::Info << "Writing...\n";
      for (size_t row = 0; row != nRows; ++row) {
        writeVisibilities(*msData.msProvider,
                          &visBuffer[row * band.ChannelCount()]);
      }
      totalNRows += nRows;
    }  // end of chunk
  }    // end of all chunks

  msData.totalRowsProcessed += totalNRows;
}

void WGriddingMSGridder::getTrimmedSize(size_t& trimmedWidth,
                                        size_t& trimmedHeight) const {
  double padding = double(ImageWidth()) / TrimWidth();
  trimmedWidth = std::floor(_actualInversionWidth / padding);
  trimmedHeight = std::floor(_actualInversionHeight / padding);
  if (trimmedWidth & 1) --trimmedWidth;
  if (trimmedHeight & 1) --trimmedHeight;
}

void WGriddingMSGridder::Invert() {
  std::vector<MSData> msDataVector;
  initializeMSDataVector(msDataVector);

  size_t trimmedWidth, trimmedHeight;
  getTrimmedSize(trimmedWidth, trimmedHeight);

  _gridder.reset(new WGriddingGridder_Simple(
      _actualInversionWidth, _actualInversionHeight, trimmedWidth,
      trimmedHeight, _actualPixelSizeX, _actualPixelSizeY, PhaseCentreDL(),
      PhaseCentreDM(), _cpuCount, _accuracy));
  _gridder->InitializeInversion();

  resetVisibilityCounters();

  for (size_t i = 0; i != MeasurementSetCount(); ++i) {
    MSData& msData = msDataVector[i];
    gridMeasurementSet(msData);
  }

  _gridder->FinalizeImage(1.0 / totalWeight());

  Logger::Info << "Gridded visibility count: "
               << double(GriddedVisibilityCount());
  if (Weighting().IsNatural())
    Logger::Info << ", effective count after weighting: "
                 << EffectiveGriddedVisibilityCount();
  Logger::Info << '\n';

  _image = Image(ActualInversionWidth(), ActualInversionHeight());
  {
    std::vector<float> imageFloat = _gridder->RealImage();
    for (size_t i = 0; i < imageFloat.size(); ++i) _image[i] = imageFloat[i];
  }

  if (ImageWidth() != _actualInversionWidth ||
      ImageHeight() != _actualInversionHeight) {
    // Interpolate the image
    // The input is of size _actualInversionWidth x _actualInversionHeight
    FFTResampler resampler(_actualInversionWidth, _actualInversionHeight,
                           ImageWidth(), ImageHeight(), _cpuCount);

    Image resized(ImageWidth(), ImageHeight());
    resampler.Resample(_image.data(), resized.data());
    _image = std::move(resized);
  }

  if (TrimWidth() != ImageWidth() || TrimHeight() != ImageHeight()) {
    Logger::Debug << "Trimming " << ImageWidth() << " x " << ImageHeight()
                  << " -> " << TrimWidth() << " x " << TrimHeight() << '\n';

    Image trimmed(TrimWidth(), TrimHeight());
    Image::Trim(trimmed.data(), TrimWidth(), TrimHeight(), _image.data(),
                ImageWidth(), ImageHeight());
    _image = std::move(trimmed);
  }
}

void WGriddingMSGridder::Predict(Image image) {
  std::vector<MSData> msDataVector;
  initializeMSDataVector(msDataVector);

  size_t trimmedWidth, trimmedHeight;
  getTrimmedSize(trimmedWidth, trimmedHeight);

  _gridder.reset(new WGriddingGridder_Simple(
      _actualInversionWidth, _actualInversionHeight, trimmedWidth,
      trimmedHeight, _actualPixelSizeX, _actualPixelSizeY, PhaseCentreDL(),
      PhaseCentreDM(), _cpuCount, _accuracy));

  if (TrimWidth() != ImageWidth() || TrimHeight() != ImageHeight()) {
    Image untrimmedImage(ImageWidth(), ImageHeight());
    Logger::Debug << "Untrimming " << TrimWidth() << " x " << TrimHeight()
                  << " -> " << ImageWidth() << " x " << ImageHeight() << '\n';
    Image::Untrim(untrimmedImage.data(), ImageWidth(), ImageHeight(),
                  image.data(), TrimWidth(), TrimHeight());
    image = std::move(untrimmedImage);
  }

  if (ImageWidth() != _actualInversionWidth ||
      ImageHeight() != _actualInversionHeight) {
    Image resampledImage(ImageWidth(), ImageHeight());
    FFTResampler resampler(ImageWidth(), ImageHeight(), _actualInversionWidth,
                           _actualInversionHeight, _cpuCount);

    resampler.Resample(image.data(), resampledImage.data());
    image = std::move(resampledImage);
  }

  std::vector<float> imageFloat(_actualInversionWidth * _actualInversionHeight);
  for (size_t i = 0; i < imageFloat.size(); ++i) imageFloat[i] = image[i];
  image.reset();

  _gridder->InitializePrediction(std::move(imageFloat));

  for (size_t i = 0; i != MeasurementSetCount(); ++i)
    predictMeasurementSet(msDataVector[i]);
}
