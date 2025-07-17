#include "wsmsgridder.h"

#include "msgriddermanager.h"

#include "../structures/imageweights.h"

#include "../system/buffered_lane.h"

#include "../msproviders/msprovider.h"
#include "../msproviders/msreaders/msreader.h"

#include <aocommon/logger.h>

#include <schaapcommon/math/resampler.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <fftw3.h>

#include <cassert>
#include <queue>
#include <stdexcept>

using aocommon::Image;
using aocommon::Logger;

namespace wsclean {

WSMSGridder::WSMSGridder(const Settings& settings, const Resources& resources,
                         MsProviderCollection& ms_provider_collection)
    : MsGridder(settings, ms_provider_collection),
      _antialiasingKernelSize(settings.antialiasingKernelSize),
      _overSamplingFactor(settings.overSamplingFactor),
      _resources(resources),
      _laneBufferSize(std::max<size_t>(_resources.NCpus() * 2, 1024)) {
  // We do this once here. WStackingGridder does this too, but by default only
  // for the float variant of fftw. schaapcommon::fft::Resampler does double
  // fft's multithreaded, hence this needs to be done here too.
  fftw_make_planner_thread_safe();
}

WSMSGridder::~WSMSGridder() noexcept {
  // In case there is an exception during gridding, the lanes need to
  // be end so that the threads stop waiting:
  for (aocommon::Lane<WSMSGridder::InversionWorkSample>& lane :
       _inversionCPULanes)
    lane.write_end();
  for (std::thread& t : _threadGroup) t.join();
}

void WSMSGridder::countSamplesPerLayer(MsProviderCollection::MsData& msData) {
  aocommon::UVector<size_t> sampleCount(ActualWGridSize(), 0);
  size_t total = 0;
  msData.matching_rows = 0;
  std::unique_ptr<MSReader> msReader = msData.ms_provider->MakeReader();
  const aocommon::MultiBandData& bands = msData.ms_provider->SelectedBands();
  while (msReader->CurrentRowAvailable()) {
    MSProvider::MetaData meta_data;
    msReader->ReadMeta(meta_data);
    const aocommon::BandData& band = bands[meta_data.data_desc_id];
    for (size_t ch = 0; ch != band.ChannelCount(); ++ch) {
      const double w = meta_data.w_in_m / band.ChannelWavelength(ch);
      const size_t wLayerIndex = _gridder->WToLayer(w);
      if (wLayerIndex < ActualWGridSize()) {
        ++sampleCount[wLayerIndex];
        ++total;
      }
    }
    ++msData.matching_rows;
    msReader->NextInputRow();
  }
  Logger::Debug << "Visibility count per layer: ";
  for (size_t& count : sampleCount) {
    Logger::Debug << count << ' ';
  }
  Logger::Debug << "\nTotal nr. of visibilities to be gridded: " << total
                << '\n';
}

size_t WSMSGridder::GetSuggestedWGridSize() const {
  size_t wWidth, wHeight;
  if (HasNWSize()) {
    wWidth = NWWidth();
    wHeight = NWHeight();
  } else {
    wWidth = TrimWidth();
    wHeight = TrimHeight();
  }
  double maxL = wWidth * PixelSizeX() * 0.5 + fabs(LShift()),
         maxM = wHeight * PixelSizeY() * 0.5 + fabs(MShift()),
         lmSq = maxL * maxL + maxM * maxM;
  double cMinW = IsComplex() ? -MaxW() : MinW();
  double radiansForAllLayers;
  if (lmSq < 1.0)
    radiansForAllLayers =
        2 * M_PI * (MaxW() - cMinW) * (1.0 - sqrt(1.0 - lmSq));
  else
    radiansForAllLayers = 2 * M_PI * (MaxW() - cMinW);
  size_t suggestedGridSize = size_t(ceil(radiansForAllLayers * NWFactor()));
  if (suggestedGridSize == 0) suggestedGridSize = 1;
  if (suggestedGridSize < _resources.NCpus()) {
    // When nwlayers is lower than the nr of cores, we cannot parallellize well.
    // However, we don't want extra w-layers if we are low on mem, as that might
    // slow down the process
    double memoryRequired =
        double(_resources.NCpus()) * double(sizeof(GridderType::num_t)) *
        double(ActualInversionWidth() * ActualInversionHeight());
    if (4.0 * memoryRequired < double(_resources.Memory())) {
      Logger::Info << "The theoretically suggested number of w-layers ("
                   << suggestedGridSize
                   << ") is less than the number of availables\n"
                      "cores ("
                   << _resources.NCpus()
                   << "). Changing suggested number of w-layers to "
                   << _resources.NCpus() << ".\n";
      suggestedGridSize = _resources.NCpus();
    } else {
      Logger::Info << "The theoretically suggested number of w-layers ("
                   << suggestedGridSize
                   << ") is less than the number of availables\n"
                      "cores ("
                   << _resources.NCpus()
                   << "), but there is not enough memory available to increase "
                      "the number of w-layers.\n"
                      "Not all cores can be used efficiently.\n";
    }
  }
  if (IsFirstTask())
    Logger::Info << "Suggested number of w-layers: " << ceil(suggestedGridSize)
                 << '\n';
  return suggestedGridSize;
}

size_t WSMSGridder::GridMeasurementSet(
    const MsProviderCollection::MsData& ms_data) {
  const size_t n_vis_polarizations = ms_data.ms_provider->NPolarizations();
  const aocommon::MultiBandData& selected_bands =
      ms_data.ms_provider->SelectedBands();

  const size_t max_data_size =
      selected_bands.MaxBandChannels() * n_vis_polarizations;
  aocommon::UVector<std::complex<float>> model_buffer(max_data_size);
  aocommon::UVector<float> weight_buffer(max_data_size);
  aocommon::UVector<bool> selection_buffer(selected_bands.MaxBandChannels());

  startInversionWorkThreads(selected_bands.MaxBandChannels());
  _gridder->PrepareBands(selected_bands);

  // Samples of the same w-layer are collected in a buffer
  // before they are written into the lane. This is done because writing
  // to a lane is reasonably slow; it requires holding a mutex. Without
  // these buffers, writing the lane was a bottleneck and multithreading
  // did not help. I think.
  std::vector<lane_write_buffer<InversionWorkSample>> buffered_lanes(
      _resources.NCpus());
  size_t lane_buffer_size =
      std::max<size_t>(8u, _inversionCPULanes[0].capacity() / 8);
  lane_buffer_size = std::min<size_t>(
      128, std::min(lane_buffer_size, _inversionCPULanes[0].capacity()));
  for (size_t i = 0; i != _resources.NCpus(); ++i) {
    buffered_lanes[i].reset(&_inversionCPULanes[i], lane_buffer_size);
  }

  InversionRow row_data;
  aocommon::UVector<std::complex<float>> row_visibilities(max_data_size);
  row_data.data = row_visibilities.data();

  const size_t n_parms = NumValuesPerSolution();
  size_t n_total_rows_read = 0;
  try {
    std::unique_ptr<MSReader> ms_reader = ms_data.ms_provider->MakeReader();
    while (ms_reader->CurrentRowAvailable()) {
      MSProvider::MetaData metadata;
      ms_reader->ReadMeta(metadata);
      const double u_in_m = metadata.u_in_m;
      const double v_in_m = metadata.v_in_m;
      const double w_in_m = metadata.w_in_m;

      const aocommon::BandData& band = selected_bands[metadata.data_desc_id];
      const double w1 = w_in_m / band.LongestWavelength();
      const double w2 = w_in_m / band.SmallestWavelength();
      if (_gridder->IsInLayerRange(w1, w2)) {
        row_data.uvw[0] = u_in_m;
        row_data.uvw[1] = v_in_m;
        row_data.uvw[2] = w_in_m;

        // Any visibilities that are not gridded in this pass
        // should not contribute to the weight sum
        for (size_t ch = 0; ch != band.ChannelCount(); ++ch) {
          const double w = row_data.uvw[2] / band.ChannelWavelength(ch);
          selection_buffer[ch] = _gridder->IsInLayerRange(w);
        }

        if (n_parms == 2) {
          GetCollapsedVisibilities<2>(*ms_reader, ms_data.antenna_names.size(),
                                      row_data, weight_buffer.data(),
                                      model_buffer.data(),
                                      selection_buffer.data(), metadata);
        } else {
          GetCollapsedVisibilities<4>(*ms_reader, ms_data.antenna_names.size(),
                                      row_data, weight_buffer.data(),
                                      model_buffer.data(),
                                      selection_buffer.data(), metadata);
        }

        if (HasDenormalPhaseCentre()) {
          const double shiftFactor =
              -2.0 * M_PI *
              (row_data.uvw[0] * LShift() + row_data.uvw[1] * MShift());
          // Because the visibilities have been collapsed, there's only one
          // polarization left:
          RotateVisibilities<1>(band, shiftFactor, row_data.data);
        }

        InversionWorkSample sample_data;
        for (size_t channel = 0; channel != band.ChannelCount(); ++channel) {
          double wavelength = band.ChannelWavelength(channel);
          sample_data.sample = row_data.data[channel];
          sample_data.uInLambda = row_data.uvw[0] / wavelength;
          sample_data.vInLambda = row_data.uvw[1] / wavelength;
          sample_data.wInLambda = row_data.uvw[2] / wavelength;
          size_t cpu =
              _gridder->WToLayer(sample_data.wInLambda) % _resources.NCpus();
          buffered_lanes[cpu].write(sample_data);
        }

        ++n_total_rows_read;
      }

      ms_reader->NextInputRow();
    }

    for (lane_write_buffer<InversionWorkSample>& lane : buffered_lanes)
      lane.write_end();

    if (IsFirstTask())
      Logger::Info << "Rows that were required: " << n_total_rows_read << '/'
                   << ms_data.matching_rows << '\n';
  } catch (...) {
    for (lane_write_buffer<InversionWorkSample>& lane : buffered_lanes)
      lane.write_end();
    throw;
  }

  finishInversionWorkThreads();
  return n_total_rows_read;
}

void WSMSGridder::startInversionWorkThreads(size_t maxChannelCount) {
  _inversionCPULanes.resize(_resources.NCpus());
  _threadGroup.clear();
  for (size_t i = 0; i != _resources.NCpus(); ++i) {
    _inversionCPULanes[i].resize(maxChannelCount * _laneBufferSize);
    set_lane_debug_name(
        _inversionCPULanes[i],
        "Work lane (buffered) containing individual visibility samples");
    _threadGroup.emplace_back(&WSMSGridder::workThreadPerSample, this,
                              &_inversionCPULanes[i]);
  }
}

void WSMSGridder::finishInversionWorkThreads() {
  for (std::thread& thrd : _threadGroup) thrd.join();
  _threadGroup.clear();
  _inversionCPULanes.clear();
}

void WSMSGridder::workThreadPerSample(
    aocommon::Lane<InversionWorkSample>* workLane) {
  size_t bufferSize = std::max<size_t>(8u, workLane->capacity() / 8);
  bufferSize =
      std::min<size_t>(128, std::min(bufferSize, workLane->capacity()));
  lane_read_buffer<InversionWorkSample> buffer(workLane, bufferSize);
  InversionWorkSample sampleData;
  while (buffer.read(sampleData)) {
    _gridder->AddDataSample(sampleData.sample, sampleData.uInLambda,
                            sampleData.vInLambda, sampleData.wInLambda);
  }
}

size_t WSMSGridder::PredictMeasurementSet(
    const MsProviderCollection::MsData& ms_data) {
  ms_data.ms_provider->ReopenRW();
  ms_data.ms_provider->ResetWritePosition();
  const aocommon::MultiBandData& selected_bands =
      ms_data.ms_provider->SelectedBands();
  _gridder->PrepareBands(selected_bands);

  size_t n_total_rows_processed = 0;

  aocommon::Lane<PredictionWorkItem> lane(_laneBufferSize + _resources.NCpus());
  aocommon::Lane<PredictionWorkItem> write_lane(_laneBufferSize);
  set_lane_debug_name(
      lane, "Prediction calculation lane (buffered) containing full row data");
  set_lane_debug_name(write_lane,
                      "Prediction write lane containing full row data");
  lane_write_buffer<PredictionWorkItem> buffered_lane(&lane, _laneBufferSize);
  std::thread writeThread(&WSMSGridder::predictWriteThread, this, &write_lane,
                          &ms_data, SelectGainMode(Polarization(), 1));
  std::vector<std::thread> calcThreads;
  for (size_t i = 0; i != _resources.NCpus(); ++i)
    calcThreads.emplace_back(&WSMSGridder::predictCalcThread, this, &lane,
                             &write_lane, &selected_bands);

  /* Start by reading the u,v,ws in, so we don't need IO access
   * from this thread during further processing */
  std::vector<std::array<double, 3>> uvws;
  std::vector<size_t> row_ids;
  std::vector<size_t> data_desc_ids;
  std::unique_ptr<MSReader> ms_reader = ms_data.ms_provider->MakeReader();
  while (ms_reader->CurrentRowAvailable()) {
    MSProvider::MetaData meta_data;
    ms_reader->ReadMeta(meta_data);
    uvws.emplace_back(std::array<double, 3>{meta_data.u_in_m, meta_data.v_in_m,
                                            meta_data.w_in_m});
    data_desc_ids.emplace_back(meta_data.data_desc_id);
    row_ids.emplace_back(ms_reader->RowId());
    ++n_total_rows_processed;

    ms_reader->NextInputRow();
  }

  for (size_t i = 0; i != uvws.size(); ++i) {
    PredictionWorkItem new_item;
    new_item.uvw = uvws[i];
    new_item.data.reset(
        new std::complex<float>[selected_bands.MaxBandChannels()]);
    new_item.rowId = row_ids[i];
    new_item.data_desc_id = data_desc_ids[i];
    buffered_lane.write(std::move(new_item));
  }
  if (IsFirstTask())
    Logger::Info << "Rows that were required: " << n_total_rows_processed << '/'
                 << ms_data.matching_rows << '\n';

  buffered_lane.write_end();
  for (std::thread& thr : calcThreads) thr.join();
  write_lane.write_end();
  writeThread.join();

  return n_total_rows_processed;
}

void WSMSGridder::predictCalcThread(
    aocommon::Lane<PredictionWorkItem>* inputLane,
    aocommon::Lane<PredictionWorkItem>* outputLane,
    const aocommon::MultiBandData* bands) {
  lane_write_buffer<PredictionWorkItem> writeBuffer(outputLane,
                                                    _laneBufferSize);

  PredictionWorkItem item;
  while (inputLane->read(item)) {
    _gridder->SampleData(item.data.get(), item.data_desc_id, item.uvw[0],
                         item.uvw[1], item.uvw[2]);
    if (HasDenormalPhaseCentre()) {
      const double shiftFactor =
          2.0 * M_PI * (item.uvw[0] * LShift() + item.uvw[1] * MShift());
      RotateVisibilities<1>((*bands)[item.data_desc_id], shiftFactor,
                            item.data.get());
    }

    writeBuffer.write(std::move(item));
  }
}

void WSMSGridder::predictWriteThread(
    aocommon::Lane<PredictionWorkItem>* predictionWorkLane,
    const MsProviderCollection::MsData* msData, GainMode gain_mode) {
  lane_read_buffer<PredictionWorkItem> buffer(
      predictionWorkLane,
      std::min(_laneBufferSize, predictionWorkLane->capacity()));
  PredictionWorkItem workItem;
  auto comparison = [](const PredictionWorkItem& lhs,
                       const PredictionWorkItem& rhs) -> bool {
    return lhs.rowId > rhs.rowId;
  };
  std::priority_queue<PredictionWorkItem, std::vector<PredictionWorkItem>,
                      decltype(comparison)>
      queue(comparison);
  size_t nextRowId = 0;
  while (buffer.read(workItem)) {
    queue.emplace(std::move(workItem));
    while (!queue.empty() && queue.top().rowId == nextRowId) {
      MSProvider::MetaData metadata;
      ReadPredictMetaData(metadata);
      WriteCollapsedVisibilities(
          *msData->ms_provider, msData->antenna_names.size(),
          metadata.data_desc_id, queue.top().data.get(), queue.top().uvw.data(),
          metadata.field_id, metadata.antenna1, metadata.antenna2,
          metadata.time);

      queue.pop();
      ++nextRowId;
    }
  }
  assert(queue.empty());
}

void WSMSGridder::StartInversion() {
  _gridder = std::make_unique<GridderType>(
      ActualInversionWidth(), ActualInversionHeight(), ActualPixelSizeX(),
      ActualPixelSizeY(), _resources.NCpus(), AntialiasingKernelSize(),
      OverSamplingFactor());
  _gridder->SetGridMode(GetGridMode());
  if (HasDenormalPhaseCentre())
    _gridder->SetDenormalPhaseCentre(LShift(), MShift());
  _gridder->SetIsComplex(IsComplex());
  //_imager->SetImageConjugatePart(Polarization() == aocommon::Polarization::YX
  //&& IsComplex());
  _gridder->PrepareWLayers(ActualWGridSize(),
                           double(_resources.Memory()) * (6.0 / 10.0), MinW(),
                           MaxW());
  if (IsFirstTask()) {
    Logger::Info << "Will process "
                 << (_gridder->NWLayers() / _gridder->NPasses()) << "/"
                 << _gridder->NWLayers() << " w-layers per pass.\n";
  }

  if (IsFirstTask() && Logger::IsVerbose()) {
    for (size_t i = 0; i != GetMsCount(); ++i)
      countSamplesPerLayer(GetMsData(i));
  }

  ResetVisibilityCounters();
}

void WSMSGridder::StartInversionPass(size_t pass_index) {
  Logger::Info << "Gridding pass " << pass_index << "... ";
  if (IsFirstTask())
    Logger::Info << '\n';
  else
    Logger::Info.Flush();
  _gridder->StartInversionPass(pass_index);
}

void WSMSGridder::FinishInversionPass(size_t pass_index) {
  Logger::Info << "Fourier transforms...\n";
  _gridder->FinishInversionPass();
}

void WSMSGridder::FinishInversion() {
  if (IsFirstTask()) {
    size_t total_rows_processed = 0;
    size_t total_rows_matching = 0;
    for (size_t i = 0; i != GetMsCount(); ++i) {
      const MsProviderCollection::MsData& ms_data = GetMsData(i);
      total_rows_processed += ms_data.total_rows_processed;
      total_rows_matching += ms_data.matching_rows;
    }

    Logger::Debug << "Total rows read: " << total_rows_processed;
    if (total_rows_matching != 0)
      Logger::Debug << " (overhead: "
                    << std::max(0.0, round(total_rows_processed * 100.0 /
                                               total_rows_matching -
                                           100.0))
                    << "%)";
    Logger::Debug << '\n';
  }

  _gridder->FinalizeImage(1.0 / ImageWeight());
  if (IsFirstTask()) {
    std::string log_message =
        "Gridded visibility count: " + std::to_string(GriddedVisibilityCount());
    if (Weighting().IsNatural()) {
      log_message += ", effective count after weighting: " +
                     std::to_string(EffectiveGriddedVisibilityCount());
    }
    Logger::Info << log_message + '\n';
  }

  _realImage = _gridder->RealImageFloat();
  if (IsComplex())
    _imaginaryImage = _gridder->ImaginaryImageFloat();
  else
    _imaginaryImage = Image();

  if (ImageWidth() != ActualInversionWidth() ||
      ImageHeight() != ActualInversionHeight()) {
    // Interpolate the image
    // The input is of size ActualInversionWidth() x ActualInversionHeight()
    schaapcommon::math::Resampler resampler(
        ActualInversionWidth(), ActualInversionHeight(), ImageWidth(),
        ImageHeight(), _resources.NCpus());

    if (IsComplex()) {
      Image resized_real(ImageWidth(), ImageHeight());
      Image resized_imaginary(ImageWidth(), ImageHeight());
      resampler.Start();
      resampler.AddTask(_realImage.Data(), resized_real.Data());
      resampler.AddTask(_imaginaryImage.Data(), resized_imaginary.Data());
      resampler.Finish();
      _realImage = std::move(resized_real);
      _imaginaryImage = std::move(resized_imaginary);
    } else {
      Image resized(ImageWidth(), ImageHeight());
      resampler.Resample(_realImage.Data(), resized.Data());
      _realImage = std::move(resized);
    }
  }

  if (TrimWidth() != ImageWidth() || TrimHeight() != ImageHeight()) {
    Logger::Debug << "Trimming " << ImageWidth() << " x " << ImageHeight()
                  << " -> " << TrimWidth() << " x " << TrimHeight() << '\n';
    _realImage = _realImage.Trim(TrimWidth(), TrimHeight());

    if (IsComplex()) {
      _imaginaryImage = _imaginaryImage.Trim(TrimWidth(), TrimHeight());
    }
  }
  Logger::Debug << "Inversion finished.\n";
}

void WSMSGridder::StartPredict(std::vector<Image>&& images) {
  if (images.size() != 2 && IsComplex())
    throw std::runtime_error("Missing imaginary in complex prediction");
  if (images.size() != 1 && !IsComplex())
    throw std::runtime_error("Imaginary specified in non-complex prediction");

  _gridder = std::make_unique<GridderType>(
      ActualInversionWidth(), ActualInversionHeight(), ActualPixelSizeX(),
      ActualPixelSizeY(), _resources.NCpus(), AntialiasingKernelSize(),
      OverSamplingFactor());
  _gridder->SetGridMode(GetGridMode());
  if (HasDenormalPhaseCentre())
    _gridder->SetDenormalPhaseCentre(LShift(), MShift());
  _gridder->SetIsComplex(IsComplex());
  //_imager->SetImageConjugatePart(Polarization() == aocommon::Polarization::YX
  //&& IsComplex());
  _gridder->PrepareWLayers(ActualWGridSize(),
                           double(_resources.Memory()) * (6.0 / 10.0), MinW(),
                           MaxW());

  if (IsFirstTask()) {
    for (size_t i = 0; i != GetMsCount(); ++i)
      countSamplesPerLayer(GetMsData(i));
  }

  if (TrimWidth() != ImageWidth() || TrimHeight() != ImageHeight()) {
    Logger::Debug << "Untrimming " << TrimWidth() << " x " << TrimHeight()
                  << " -> " << ImageWidth() << " x " << ImageHeight() << '\n';
    // Undo trimming (i.e., extend with zeros)
    // The input is of size TrimWidth() x TrimHeight()
    // This will make the model image of size ImageWidth() x ImageHeight()
    for (Image& image : images) {
      image = image.Untrim(ImageWidth(), ImageHeight());
    }
  }

  if (ImageWidth() != ActualInversionWidth() ||
      ImageHeight() != ActualInversionHeight()) {
    // Decimate the image
    // Input is ImageWidth() x ImageHeight()
    schaapcommon::math::Resampler resampler(
        ImageWidth(), ImageHeight(), ActualInversionWidth(),
        ActualInversionHeight(), _resources.NCpus());

    if (images.size() == 1) {
      Image resampled(ImageWidth(), ImageHeight());
      resampler.Resample(images[0].Data(), resampled.Data());
      images[0] = std::move(resampled);
    } else {
      std::vector<Image> resampled;
      resampled.reserve(images.size());
      resampler.Start();
      for (Image& image : images) {
        resampled.emplace_back(ImageWidth(), ImageHeight());
        resampler.AddTask(image.Data(), resampled.back().Data());
      }
      resampler.Finish();
      for (size_t i = 0; i != images.size(); ++i)
        images[i] = std::move(resampled[i]);
    }
  }

  if (images.size() == 1)
    _gridder->InitializePrediction(images[0]);
  else
    _gridder->InitializePrediction(images[0], images[1]);
}

void WSMSGridder::StartPredictPass(size_t pass_index) {
  Logger::Info << "Fourier transforms for pass " << pass_index << "... ";
  if (IsFirstTask())
    Logger::Info << '\n';
  else
    Logger::Info.Flush();

  _gridder->StartPredictionPass(pass_index);

  Logger::Info << "Predicting...\n";
}

void WSMSGridder::FinishPredictPass(size_t /*pass_index*/) {}

void WSMSGridder::FinishPredict() {
  size_t total_rows_processed = 0;
  size_t total_rows_matching = 0;
  for (size_t i = 0; i != GetMsCount(); ++i) {
    const MsProviderCollection::MsData& ms_data = GetMsData(i);
    total_rows_processed += ms_data.total_rows_processed;
    total_rows_matching += ms_data.matching_rows;
  }

  Logger::Debug << "Total rows written: " << total_rows_processed;
  if (total_rows_matching != 0)
    Logger::Debug << " (overhead: "
                  << std::max(0.0, round(total_rows_processed * 100.0 /
                                             total_rows_matching -
                                         100.0))
                  << "%)";
  Logger::Debug << '\n';
}

}  // namespace wsclean
