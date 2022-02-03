
#include "wgriddingmsgridder.h"

#include "wgriddinggridder_simple.h"

#include "../msproviders/msreaders/msreader.h"

#include "../io/logger.h"

#include "../math/fftresampler.h"

#include "../msproviders/msprovider.h"

#include "../system/buffered_lane.h"

#include "../structures/imageweights.h"

#include <aocommon/image.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

using aocommon::Image;

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

template <DDGainMatrix GainEntry>
void WGriddingMSGridder::gridMeasurementSet(MSData& msData) {
  const aocommon::BandData selectedBand(msData.SelectedBand());
  StartMeasurementSet(msData, false);

  aocommon::UVector<std::complex<float>> modelBuffer(
      selectedBand.ChannelCount());
  aocommon::UVector<float> weightBuffer(selectedBand.ChannelCount());
  aocommon::UVector<bool> isSelected(selectedBand.ChannelCount(), true);

  size_t totalNRows = 0;
  aocommon::UVector<double> frequencies(selectedBand.ChannelCount());
  for (size_t i = 0; i != frequencies.size(); ++i)
    frequencies[i] = selectedBand.ChannelFrequency(i);

  size_t maxNRows = calculateMaxNRowsInMemory(selectedBand.ChannelCount());

  aocommon::UVector<std::complex<float>> visBuffer(maxNRows *
                                                   selectedBand.ChannelCount());
  aocommon::UVector<double> uvwBuffer(maxNRows * 3);

  std::unique_ptr<MSReader> msReader = msData.msProvider->MakeReader();
  aocommon::UVector<std::complex<float>> newItemData(
      selectedBand.ChannelCount());
  InversionRow newRowData;
  newRowData.data = newItemData.data();

  // Iterate over chunks until all data has been gridded
  while (msReader->CurrentRowAvailable()) {
    Logger::Debug << "Max " << maxNRows << " rows fit in memory.\n";
    Logger::Info << "Loading data in memory...\n";

    size_t nRows = 0;

    // Read / fill the chunk
    while (msReader->CurrentRowAvailable() && nRows < maxNRows) {
      double uInMeters, vInMeters, wInMeters;
      msReader->ReadMeta(uInMeters, vInMeters, wInMeters);
      newRowData.uvw[0] = uInMeters;
      newRowData.uvw[1] = vInMeters;
      newRowData.uvw[2] = wInMeters;
      readAndWeightVisibilities<1, GainEntry>(
          *msReader, msData.antennaNames, newRowData, selectedBand,
          weightBuffer.data(), modelBuffer.data(), isSelected.data());

      std::copy_n(newRowData.data, selectedBand.ChannelCount(),
                  &visBuffer[nRows * selectedBand.ChannelCount()]);
      std::copy_n(newRowData.uvw, 3, &uvwBuffer[nRows * 3]);

      ++nRows;
      msReader->NextInputRow();
    }

    Logger::Info << "Gridding " << nRows << " rows...\n";
    _gridder->AddInversionData(nRows, selectedBand.ChannelCount(),
                               uvwBuffer.data(), frequencies.data(),
                               visBuffer.data());

    totalNRows += nRows;
  }  // end of chunk

  msData.totalRowsProcessed += totalNRows;
}

template void WGriddingMSGridder::gridMeasurementSet<DDGainMatrix::kXX>(
    MSData& msData);
template void WGriddingMSGridder::gridMeasurementSet<DDGainMatrix::kYY>(
    MSData& msData);
template void WGriddingMSGridder::gridMeasurementSet<DDGainMatrix::kTrace>(
    MSData& msData);

template <DDGainMatrix GainEntry>
void WGriddingMSGridder::predictMeasurementSet(MSData& msData) {
  msData.msProvider->ReopenRW();
  const aocommon::BandData selectedBand(msData.SelectedBand());
  StartMeasurementSet(msData, true);

  size_t totalNRows = 0;

  aocommon::UVector<double> frequencies(selectedBand.ChannelCount());
  for (size_t i = 0; i != frequencies.size(); ++i)
    frequencies[i] = selectedBand.ChannelFrequency(i);

  size_t maxNRows = calculateMaxNRowsInMemory(selectedBand.ChannelCount());

  aocommon::UVector<double> uvwBuffer(maxNRows * 3);
  // Iterate over chunks until all data has been gridded
  msData.msProvider->ResetWritePosition();
  std::unique_ptr<MSReader> msReader = msData.msProvider->MakeReader();
  while (msReader->CurrentRowAvailable()) {
    size_t nRows = 0;
    // Read / fill the chunk
    while (msReader->CurrentRowAvailable() && nRows < maxNRows) {
      double uInMeters, vInMeters, wInMeters;
      msReader->ReadMeta(uInMeters, vInMeters, wInMeters);
      uvwBuffer[nRows * 3] = uInMeters;
      uvwBuffer[nRows * 3 + 1] = vInMeters;
      uvwBuffer[nRows * 3 + 2] = wInMeters;
      ++nRows;
      msReader->NextInputRow();
    }

    Logger::Info << "Predicting " << nRows << " rows...\n";
    aocommon::UVector<std::complex<float>> visBuffer(
        maxNRows * selectedBand.ChannelCount());
    _gridder->PredictVisibilities(nRows, selectedBand.ChannelCount(),
                                  uvwBuffer.data(), frequencies.data(),
                                  visBuffer.data());

    Logger::Info << "Writing...\n";
    for (size_t row = 0; row != nRows; ++row) {
      writeVisibilities<1, GainEntry>(
          *msData.msProvider, msData.antennaNames, selectedBand,
          &visBuffer[row * selectedBand.ChannelCount()]);
    }
    totalNRows += nRows;
  }  // end of chunk

  msData.totalRowsProcessed += totalNRows;
}

template void WGriddingMSGridder::predictMeasurementSet<DDGainMatrix::kXX>(
    MSData& msData);
template void WGriddingMSGridder::predictMeasurementSet<DDGainMatrix::kYY>(
    MSData& msData);
template void WGriddingMSGridder::predictMeasurementSet<DDGainMatrix::kTrace>(
    MSData& msData);

void WGriddingMSGridder::getActualTrimmedSize(size_t& trimmedWidth,
                                              size_t& trimmedHeight) const {
  trimmedWidth = std::ceil(ActualInversionWidth() / ImagePadding());
  trimmedHeight = std::ceil(ActualInversionHeight() / ImagePadding());

  // In facet-based imaging, the alignment is 4, see wsclean.cpp. Also for
  // monolithic imaging - in which just an even number would suffice -
  // the trimmedWidth and trimmedHeight are defined to be divisable by 4.
  const size_t alignment = 4;
  if (trimmedWidth % alignment != 0) {
    trimmedWidth += alignment - (trimmedWidth % alignment);
  }
  if (trimmedHeight % alignment != 0) {
    trimmedHeight += alignment - (trimmedHeight % alignment);
  }
  trimmedWidth = std::min(trimmedWidth, ActualInversionWidth());
  trimmedHeight = std::min(trimmedHeight, ActualInversionHeight());
}

void WGriddingMSGridder::Invert() {
  std::vector<MSData> msDataVector;
  initializeMSDataVector(msDataVector);

  size_t trimmedWidth, trimmedHeight;
  getActualTrimmedSize(trimmedWidth, trimmedHeight);

  _gridder.reset(new WGriddingGridder_Simple(
      ActualInversionWidth(), ActualInversionHeight(), trimmedWidth,
      trimmedHeight, ActualPixelSizeX(), ActualPixelSizeY(), PhaseCentreDL(),
      PhaseCentreDM(), _cpuCount, _accuracy));
  _gridder->InitializeInversion();

  resetVisibilityCounters();

  for (size_t i = 0; i != MeasurementSetCount(); ++i) {
    MSData& msData = msDataVector[i];
    if (Polarization() == aocommon::Polarization::XX) {
      gridMeasurementSet<DDGainMatrix::kXX>(msData);
    } else if (Polarization() == aocommon::Polarization::YY) {
      gridMeasurementSet<DDGainMatrix::kYY>(msData);
    } else {
      gridMeasurementSet<DDGainMatrix::kTrace>(msData);
    }
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

  if (ImageWidth() != ActualInversionWidth() ||
      ImageHeight() != ActualInversionHeight()) {
    // Interpolate the image
    // The input is of size ActualInversionWidth() x ActualInversionHeight()
    FFTResampler resampler(ActualInversionWidth(), ActualInversionHeight(),
                           ImageWidth(), ImageHeight(), _cpuCount);

    Image resized(ImageWidth(), ImageHeight());
    resampler.Resample(_image.Data(), resized.Data());
    _image = std::move(resized);
  }

  if (TrimWidth() != ImageWidth() || TrimHeight() != ImageHeight()) {
    Logger::Debug << "Trimming " << ImageWidth() << " x " << ImageHeight()
                  << " -> " << TrimWidth() << " x " << TrimHeight() << '\n';

    Image trimmed(TrimWidth(), TrimHeight());
    Image::Trim(trimmed.Data(), TrimWidth(), TrimHeight(), _image.Data(),
                ImageWidth(), ImageHeight());
    _image = std::move(trimmed);
  }
}

void WGriddingMSGridder::Predict(std::vector<Image>&& images) {
  std::vector<MSData> msDataVector;
  initializeMSDataVector(msDataVector);

  size_t trimmedWidth, trimmedHeight;
  getActualTrimmedSize(trimmedWidth, trimmedHeight);

  _gridder.reset(new WGriddingGridder_Simple(
      ActualInversionWidth(), ActualInversionHeight(), trimmedWidth,
      trimmedHeight, ActualPixelSizeX(), ActualPixelSizeY(), PhaseCentreDL(),
      PhaseCentreDM(), _cpuCount, _accuracy));

  if (TrimWidth() != ImageWidth() || TrimHeight() != ImageHeight()) {
    Image untrimmedImage(ImageWidth(), ImageHeight());
    Logger::Debug << "Untrimming " << TrimWidth() << " x " << TrimHeight()
                  << " -> " << ImageWidth() << " x " << ImageHeight() << '\n';
    Image::Untrim(untrimmedImage.Data(), ImageWidth(), ImageHeight(),
                  images[0].Data(), TrimWidth(), TrimHeight());
    images[0] = std::move(untrimmedImage);
  }

  if (ImageWidth() != ActualInversionWidth() ||
      ImageHeight() != ActualInversionHeight()) {
    Image resampledImage(ImageWidth(), ImageHeight());
    FFTResampler resampler(ImageWidth(), ImageHeight(), ActualInversionWidth(),
                           ActualInversionHeight(), _cpuCount);

    resampler.Resample(images[0].Data(), resampledImage.Data());
    images[0] = std::move(resampledImage);
  }

  _gridder->InitializePrediction(images[0].Data());
  images[0].Reset();

  for (MSData& msData : msDataVector) {
    if (Polarization() == aocommon::Polarization::XX) {
      predictMeasurementSet<DDGainMatrix::kXX>(msData);
    } else if (Polarization() == aocommon::Polarization::YY) {
      predictMeasurementSet<DDGainMatrix::kYY>(msData);
    } else {
      predictMeasurementSet<DDGainMatrix::kTrace>(msData);
    }
  }
}
