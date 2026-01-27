#include "wgriddingmsgridder.h"

#include "wgridder.h"

#include "../gridding/msgriddermanager.h"
#include "../msproviders/msreaders/msreader.h"

#include "../msproviders/msprovider.h"

#include "../system/buffered_lane.h"

#include "../structures/imageweights.h"

#include <aocommon/image.h>
#include <aocommon/logger.h>

#include <schaapcommon/math/resampler.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

using aocommon::Image;
using aocommon::Logger;

namespace wsclean {
namespace {
constexpr size_t kUvwSize = sizeof(double) * 3;
}  // namespace

WGriddingMSGridder::WGriddingMSGridder(
    const Settings& settings, const Resources& resources,
    MsProviderCollection& ms_provider_collection, bool use_tuned_wgridder)
    : MsGridder(settings, ms_provider_collection),
      resources_(resources),
      accuracy_(GetSettings().gridder_accuracy),
      use_tuned_wgridder_(use_tuned_wgridder) {
  // It may happen that several schaapcommon::fft::Resamplers are created
  // concurrently, so we must make sure that the FFTW planner can deal with
  // this.
  fftwf_make_planner_thread_safe();
}

WGriddingMSGridder::~WGriddingMSGridder() = default;

std::unique_ptr<WGridderBase> WGriddingMSGridder::MakeGridder(
    size_t width, size_t height) const {
  if (accuracy_ <= 1.01e-5) {
    return std::make_unique<WGridder<double>>(
        ActualInversionWidth(), ActualInversionHeight(), width, height,
        ActualPixelSizeX(), ActualPixelSizeY(), LShift(), MShift(),
        resources_.NCpus(), accuracy_, 0, use_tuned_wgridder_);
  } else {
    return std::make_unique<WGridder<float>>(
        ActualInversionWidth(), ActualInversionHeight(), width, height,
        ActualPixelSizeX(), ActualPixelSizeY(), LShift(), MShift(),
        resources_.NCpus(), accuracy_, 0, use_tuned_wgridder_);
  }
}

size_t WGriddingMSGridder::CalculateConstantMemory() const {
  size_t constant_mem = gridder_->ConstantMemoryUsage();
  constant_mem += GetVisibilityModifier().GetCacheParmResponseSize();
  return constant_mem;
}

size_t WGriddingMSGridder::CalculateMaxRowsInMemory(
    int64_t available_memory, size_t constant_memory,
    double additional_per_row_consumption, size_t per_row_uvw_consumption,
    size_t channel_count, size_t num_polarizations_stored) const {
  if (static_cast<int64_t>(constant_memory) >= available_memory) {
    // Assume that half the memory is necessary for the constant parts (like
    // image grid), and the other half remains available for the dynamic buffers
    constant_memory = available_memory / 2;
    Logger::Warn << "Not enough memory available for doing the gridding:\n"
                    "swapping might occur!\n";
  }

  const size_t per_visibility_ducc_overhead =
      gridder_->PerVisibilityMemoryUsage();
  // Internal DUCC size and size we might be storing in memory are different as
  // we collapse the visibilities before passing them in to DUCC
  const size_t per_row_visibility_memory =
      (per_visibility_ducc_overhead * channel_count) +
      (sizeof(std::complex<float>) * channel_count * num_polarizations_stored);
  // Keep computed number floating point to maintain precision.
  // This is because additional_per_row_consumption can be fractional in the
  // case of -shared-facet-writes; where gridders/facets share some per row
  // memory overheads with each other. In this instance we allocate a fractional
  // portion of this shared memory to each gridder to help compute maximum rows
  // more accurately.
  const double memory_per_row =
      additional_per_row_consumption  // external overheads
      + per_row_visibility_memory     // visibilities
      + per_row_uvw_consumption;      // uvw
  const uint64_t memory_for_buffers = available_memory - constant_memory;
  const size_t max_n_rows =
      std::max(uint64_t(memory_for_buffers / memory_per_row), uint64_t(100));
  if (max_n_rows < 1000) {
    Logger::Warn << "Less than 1000 data rows fit in memory: this probably "
                    "means performance is going to be very poor!\n";
  }

  return max_n_rows;
}

size_t WGriddingMSGridder::CalculateMaxVisibilitiesInMemory(
    int64_t available_memory, size_t constant_memory,
    double additional_per_visibility_consumption,
    size_t per_visibility_uvw_consumption,
    size_t num_polarizations_stored) const {
  // Function follows mostly the logic of CalculateMaxRowsInMemory().
  if (static_cast<int64_t>(constant_memory) >= available_memory) {
    constant_memory = available_memory / 2;
    Logger::Warn << "Not enough memory available for doing the gridding:\n"
                    "swapping might occur!\n";
  }

  const size_t per_visibility_ducc_overhead =
      gridder_->PerVisibilityMemoryUsage();
  const size_t per_visibility_memory =
      per_visibility_ducc_overhead +
      (sizeof(std::complex<float>) * num_polarizations_stored);
  const double memory_per_visibility =
      additional_per_visibility_consumption  // external overheads
      + per_visibility_memory                // visibilities
      + per_visibility_uvw_consumption;      // uvw
  const uint64_t memory_for_buffers = available_memory - constant_memory;
  // This value is a bit arbitrary, but gridding less than 10000 vis at a time
  // will be prohabitively slow...
  constexpr uint64_t kMinVisibilities = 10000;
  const size_t max_n_visibilities = std::max(
      uint64_t(memory_for_buffers / memory_per_visibility), kMinVisibilities);

  return max_n_visibilities;
}

void WGriddingMSGridder::GridSharedMeasurementSetChunk(
    bool apply_corrections, size_t n_polarizations, size_t n_rows,
    const double* uvws, const double* frequencies,
    const aocommon::BandData& selected_band, size_t data_desc_id,
    const std::pair<size_t, size_t>* antennas,
    const std::complex<float>* visibilities, const size_t* time_offsets,
    size_t n_antennas, const std::vector<std::complex<float>>& parm_response,
    const BeamResponseCacheChunk& beam_response) {
  // If there are no corrections to apply then we can bypass needing a callback
  // and just use the shared buffer directly
  if (!apply_corrections) {
    gridder_->AddInversionData(n_rows, selected_band.ChannelCount(), uvws,
                               frequencies, visibilities);
  } else {
    VisibilityCallbackData data(selected_band.ChannelCount(), selected_band,
                                data_desc_id, antennas, visibilities, uvws,
                                time_offsets, this, n_antennas,
                                parm_response.data(), beam_response);
    gridder_->AddInversionDataWithCorrectionCallback(
        GetGainMode(), n_polarizations, n_rows, uvws, frequencies, data);
  }
}

size_t WGriddingMSGridder::GridRegularMeasurementSet(
    const MsProviderCollection::MsData& ms_data) {
  // Regular data should always have one band per ms provider...
  const aocommon::MultiBandData& selected_bands(
      ms_data.ms_provider->SelectedBands());
  assert(selected_bands.BandCount() == 1);
  const aocommon::BandData& selected_band = *selected_bands.begin();
  const size_t n_vis_polarizations = ms_data.ms_provider->NPolarizations();

  const size_t data_size = selected_band.ChannelCount() * n_vis_polarizations;
  aocommon::UVector<std::complex<float>> model_buffer(data_size);
  aocommon::UVector<float> weight_buffer(data_size);
  aocommon::UVector<bool> selection_buffer(selected_band.ChannelCount(), true);

  aocommon::UVector<double> frequencies(selected_band.ChannelCount());
  for (size_t i = 0; i != frequencies.size(); ++i)
    frequencies[i] = selected_band.ChannelFrequency(i);

  size_t max_rows_per_chunk =
      CalculateMaxRowsInMemory(resources_.Memory(), CalculateConstantMemory(),
                               0, kUvwSize, selected_band.ChannelCount(), 1);

  aocommon::UVector<std::complex<float>> visibility_buffer(
      max_rows_per_chunk * selected_band.ChannelCount());
  aocommon::UVector<double> uvw_buffer(max_rows_per_chunk * 3);

  std::unique_ptr<MSReader> ms_reader = ms_data.ms_provider->MakeReader();
  aocommon::UVector<std::complex<float>> row_visibilities(data_size);
  InversionRow row_data;
  row_data.data = row_visibilities.data();

  const size_t n_parms = NumValuesPerSolution();

  // Iterate over chunks until all data has been gridded
  size_t n_total_rows_read = 0;
  Logger::Debug << "Max " << max_rows_per_chunk << " rows fit in memory.\n";
  while (ms_reader->CurrentRowAvailable()) {
    Logger::Info << "Loading data in memory...\n";

    size_t n_chunk_rows_read = 0;

    // Read / fill the chunk
    while (ms_reader->CurrentRowAvailable() &&
           n_chunk_rows_read < max_rows_per_chunk) {
      MSProvider::MetaData metadata;
      ms_reader->ReadMeta(metadata);
      row_data.uvw[0] = metadata.u_in_m;
      row_data.uvw[1] = metadata.v_in_m;
      row_data.uvw[2] = metadata.w_in_m;

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

      std::copy_n(
          row_data.data, selected_band.ChannelCount(),
          &visibility_buffer[n_chunk_rows_read * selected_band.ChannelCount()]);
      std::copy_n(row_data.uvw, 3, &uvw_buffer[n_chunk_rows_read * 3]);

      ++n_chunk_rows_read;
      ms_reader->NextInputRow();
    }

    Logger::Info << "Gridding " << n_chunk_rows_read << " rows...\n";
    gridder_->AddInversionData(n_chunk_rows_read, selected_band.ChannelCount(),
                               uvw_buffer.data(), frequencies.data(),
                               visibility_buffer.data());

    n_total_rows_read += n_chunk_rows_read;
  }  // end of chunk
  return n_total_rows_read;
}

size_t WGriddingMSGridder::GridBdaMeasurementSet(
    const MsProviderCollection::MsData& ms_data) {
  const aocommon::MultiBandData& selected_bands(
      ms_data.ms_provider->SelectedBands());
  const size_t n_vis_polarizations = ms_data.ms_provider->NPolarizations();

  const size_t max_data_size =
      ms_data.ms_provider->NMaxChannels() * n_vis_polarizations;
  aocommon::UVector<std::complex<float>> model_buffer(max_data_size);
  aocommon::UVector<float> weight_buffer(max_data_size);
  aocommon::UVector<bool> selection_buffer(max_data_size, true);

  const size_t max_vis_per_chunk = CalculateMaxVisibilitiesInMemory(
      resources_.Memory(), CalculateConstantMemory(), 0, kUvwSize, 1);

  aocommon::UVector<std::complex<float>> visibility_buffer;
  visibility_buffer.reserve(max_vis_per_chunk);
  aocommon::UVector<double> uvw_buffer;
  uvw_buffer.reserve(3 * max_vis_per_chunk);

  std::unique_ptr<MSReader> ms_reader = ms_data.ms_provider->MakeReader();
  aocommon::UVector<std::complex<float>> row_visibilities(max_data_size);
  InversionRow row_data;
  row_data.data = row_visibilities.data();

  const size_t n_parms = NumValuesPerSolution();

  // Iterate over chunks until all data has been gridded
  size_t n_total_rows_read = 0;
  while (ms_reader->CurrentRowAvailable()) {
    Logger::Info << "Loading data in memory...\n";

    size_t n_chunk_rows_read = 0;
    visibility_buffer.clear();
    uvw_buffer.clear();

    // Read / fill the chunk
    while (ms_reader->CurrentRowAvailable()) {
      MSProvider::MetaData metadata;
      ms_reader->ReadMeta(metadata);
      row_data.uvw[0] = metadata.u_in_m;
      row_data.uvw[1] = metadata.v_in_m;
      row_data.uvw[2] = metadata.w_in_m;
      const aocommon::BandData& band = selected_bands[metadata.data_desc_id];
      if (visibility_buffer.size() + band.ChannelCount() >= max_vis_per_chunk)
        break;

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

      for (size_t channel = 0; channel != band.ChannelCount(); ++channel) {
        // Because the gridder doesn't have an option to have different nr of
        // channels per row, the data is "flattened" into a single array and the
        // uvws are scaled so that they become frequency independent.
        visibility_buffer.emplace_back(row_data.data[channel]);
        for (size_t i = 0; i != 3; ++i)
          uvw_buffer.emplace_back(row_data.uvw[i] *
                                  band.ChannelFrequency(channel));
      }

      ++n_chunk_rows_read;
      ms_reader->NextInputRow();
    }

    Logger::Info << "Gridding " << n_chunk_rows_read
                 << " (irregular) rows...\n";

    // The Uvws have been scaled by the frequency already, hence use 1 here.
    constexpr double kDummyFrequency = 1.0;
    // Data has been flattened, so for the gridder there's only one channel:
    constexpr size_t kNGridderChannels = 1;
    gridder_->AddInversionData(visibility_buffer.size(), kNGridderChannels,
                               uvw_buffer.data(), &kDummyFrequency,
                               visibility_buffer.data());

    n_total_rows_read += n_chunk_rows_read;
  }  // end of chunk
  return n_total_rows_read;
}

void WGriddingMSGridder::PredictChunk(size_t n_rows, size_t n_channels,
                                      const double* frequencies,
                                      const double* uvws,
                                      std::complex<float>* visibilities) const {
  Logger::Info << "Predicting " + std::to_string(n_rows) + " rows...\n";
  gridder_->PredictVisibilities(n_rows, n_channels, uvws, frequencies,
                                visibilities);
}

size_t WGriddingMSGridder::PredictRegularMeasurementSet(
    const MsProviderCollection::MsData& ms_data) {
  const aocommon::MultiBandData& selected_bands(
      ms_data.ms_provider->SelectedBands());
  // Regular data should always have one band...
  assert(selected_bands.BandCount() == 1);
  const aocommon::BandData& selected_band = *selected_bands.begin();

  size_t n_total_rows_written = 0;

  aocommon::UVector<double> frequencies(selected_band.ChannelCount());
  for (size_t i = 0; i != frequencies.size(); ++i)
    frequencies[i] = selected_band.ChannelFrequency(i);

  const size_t max_rows_per_chunk =
      CalculateMaxRowsInMemory(resources_.Memory(), CalculateConstantMemory(),
                               0, kUvwSize, selected_band.ChannelCount(), 1);

  aocommon::UVector<double> uvw_buffer(max_rows_per_chunk * 3);
  // Iterate over chunks until all data has been gridded
  ms_data.ms_provider->ResetWritePosition();
  std::unique_ptr<MSReader> ms_reader = ms_data.ms_provider->MakeReader();
  while (ms_reader->CurrentRowAvailable()) {
    size_t n_chunk_rows_read = 0;

    // Read / fill the chunk
    Logger::Info << "Loading metadata...\n";
    // Read from metadata buffer
    std::vector<MSProvider::MetaData> metadata_buffer;
    while (ms_reader->CurrentRowAvailable() &&
           n_chunk_rows_read < max_rows_per_chunk) {
      MSProvider::MetaData metadata;
      ReadPredictMetaData(metadata);
      uvw_buffer[n_chunk_rows_read * 3] = metadata.u_in_m;
      uvw_buffer[n_chunk_rows_read * 3 + 1] = metadata.v_in_m;
      uvw_buffer[n_chunk_rows_read * 3 + 2] = metadata.w_in_m;
      metadata_buffer.emplace_back(std::move(metadata));
      n_chunk_rows_read++;

      ms_reader->NextInputRow();
    }

    Logger::Info << "Predicting " << n_chunk_rows_read << " rows...\n";
    aocommon::UVector<std::complex<float>> visibility_buffer(
        max_rows_per_chunk * selected_band.ChannelCount());
    gridder_->PredictVisibilities(
        n_chunk_rows_read, selected_band.ChannelCount(), uvw_buffer.data(),
        frequencies.data(), visibility_buffer.data());

    Logger::Info << "Writing...\n";
    for (size_t row = 0; row != n_chunk_rows_read; ++row) {
      WriteCollapsedVisibilities(
          *ms_data.ms_provider, ms_data.antenna_names.size(),
          metadata_buffer[row].data_desc_id,
          &visibility_buffer[row * selected_band.ChannelCount()],
          &uvw_buffer[row * 3], metadata_buffer[row].field_id,
          metadata_buffer[row].antenna1, metadata_buffer[row].antenna2,
          metadata_buffer[row].time);
    }
    n_total_rows_written += n_chunk_rows_read;
  }  // end of chunk
  return n_total_rows_written;
}

size_t WGriddingMSGridder::PredictBdaMeasurementSet(
    const MsProviderCollection::MsData& ms_data) {
  const aocommon::MultiBandData& selected_bands(
      ms_data.ms_provider->SelectedBands());

  size_t n_total_rows_written = 0;

  const size_t max_vis_per_chunk = CalculateMaxVisibilitiesInMemory(
      resources_.Memory(), CalculateConstantMemory(), 0, kUvwSize, 1);

  aocommon::UVector<double> uvw_buffer;
  uvw_buffer.reserve(3 * max_vis_per_chunk);

  // Iterate over chunks until all data has been gridded
  ms_data.ms_provider->ResetWritePosition();
  std::unique_ptr<MSReader> ms_reader = ms_data.ms_provider->MakeReader();
  std::vector<MSProvider::MetaData> metadata_buffer;
  while (ms_reader->CurrentRowAvailable()) {
    size_t n_chunk_rows_read = 0;

    // Read / fill the chunk
    Logger::Info << "Loading metadata...\n";
    metadata_buffer.clear();
    uvw_buffer.clear();
    size_t visibility_count = 0;
    while (ms_reader->CurrentRowAvailable()) {
      MSProvider::MetaData metadata;
      ReadPredictMetaData(metadata);
      const aocommon::BandData& band = selected_bands[metadata.data_desc_id];
      if (visibility_count + band.ChannelCount() >= max_vis_per_chunk) break;

      visibility_count += band.ChannelCount();
      for (size_t channel = 0; channel != band.ChannelCount(); ++channel) {
        const double frequency = band.ChannelFrequency(channel);
        uvw_buffer.emplace_back(metadata.u_in_m * frequency);
        uvw_buffer.emplace_back(metadata.v_in_m * frequency);
        uvw_buffer.emplace_back(metadata.w_in_m * frequency);
      }
      metadata_buffer.emplace_back(std::move(metadata));
      n_chunk_rows_read++;

      ms_reader->NextInputRow();
    }

    Logger::Info << "Predicting " << n_chunk_rows_read
                 << " (irregular) rows...\n";
    aocommon::UVector<std::complex<float>> visibility_buffer(visibility_count);
    // The Uvws have been scaled by the frequency already, hence use 1 here.
    constexpr double kDummyFrequency = 1.0;
    // Data has been flattened, so for the gridder there's only one channel:
    constexpr size_t kNGridderChannels = 1;
    gridder_->PredictVisibilities(visibility_count, kNGridderChannels,
                                  uvw_buffer.data(), &kDummyFrequency,
                                  visibility_buffer.data());

    Logger::Info << "Writing...\n";
    std::complex<float>* visibility_ptr = visibility_buffer.data();
    for (size_t row = 0; row != n_chunk_rows_read; ++row) {
      // The uvw_buffer has scaled uvws, so reload them from the metadata
      const double uvw[3] = {metadata_buffer[row].u_in_m,
                             metadata_buffer[row].v_in_m,
                             metadata_buffer[row].w_in_m};
      WriteCollapsedVisibilities(
          *ms_data.ms_provider, ms_data.antenna_names.size(),
          metadata_buffer[row].data_desc_id, visibility_ptr, uvw,
          metadata_buffer[row].field_id, metadata_buffer[row].antenna1,
          metadata_buffer[row].antenna2, metadata_buffer[row].time);
      const aocommon::BandData& band =
          selected_bands[metadata_buffer[row].data_desc_id];
      visibility_ptr += band.ChannelCount();
    }
    n_total_rows_written += n_chunk_rows_read;
  }  // end of chunk
  return n_total_rows_written;
}

size_t WGriddingMSGridder::GridMeasurementSet(
    const MsProviderCollection::MsData& ms_data) {
  // If the data isn't regular, the data is flattened before calling the
  // gridder. This costs more memory and may potentially be a bit slower or
  // cause more gridder calls, so flattening is done only when it is necessary.
  if (ms_data.ms_provider->IsRegular()) {
    return GridRegularMeasurementSet(ms_data);
  } else {
    return GridBdaMeasurementSet(ms_data);
  }
}

size_t WGriddingMSGridder::PredictMeasurementSet(
    const MsProviderCollection::MsData& ms_data) {
  ms_data.ms_provider->ReopenRW();

  // See comment in GridMeasurementSet().
  if (ms_data.ms_provider->IsRegular()) {
    return PredictRegularMeasurementSet(ms_data);
  } else {
    return PredictBdaMeasurementSet(ms_data);
  }
}

void WGriddingMSGridder::GetActualTrimmedSize(size_t& trimmedWidth,
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

void WGriddingMSGridder::StartInversion() {
  size_t trimmed_width;
  size_t trimmed_height;
  GetActualTrimmedSize(trimmed_width, trimmed_height);

  gridder_ = MakeGridder(trimmed_width, trimmed_height);
  gridder_->InitializeInversion();

  ResetVisibilityCounters();
}

void WGriddingMSGridder::FinishInversion() {
  gridder_->FinalizeImage(1.0 / ImageWeight());

  std::string log_message =
      "Gridded visibility count: " + std::to_string(GriddedVisibilityCount());
  if (Weighting().IsNatural()) {
    log_message += ", effective count after weighting: " +
                   std::to_string(EffectiveGriddedVisibilityCount());
  }
  Logger::Info << log_message + '\n';

  image_ = Image(ActualInversionWidth(), ActualInversionHeight());
  {
    std::vector<float> image_float = gridder_->RealImage();
    std::copy(image_float.begin(), image_float.end(), image_.begin());
  }

  if (ImageWidth() != ActualInversionWidth() ||
      ImageHeight() != ActualInversionHeight()) {
    // Interpolate the image
    // The input is of size ActualInversionWidth() x ActualInversionHeight()
    schaapcommon::math::Resampler resampler(
        ActualInversionWidth(), ActualInversionHeight(), ImageWidth(),
        ImageHeight(), resources_.NCpus());

    Image resized(ImageWidth(), ImageHeight());
    resampler.Resample(image_.Data(), resized.Data());
    image_ = std::move(resized);
  }

  if (TrimWidth() != ImageWidth() || TrimHeight() != ImageHeight()) {
    Logger::Debug << "Trimming " + std::to_string(ImageWidth()) + " x " +
                         std::to_string(ImageHeight()) + " -> " +
                         std::to_string(TrimWidth()) + " x " +
                         std::to_string(TrimHeight()) + "\n";

    image_ = image_.Trim(TrimWidth(), TrimHeight());
  }
}

void WGriddingMSGridder::StartPredict(std::vector<Image>&& images) {
  size_t trimmed_width;
  size_t trimmed_height;
  GetActualTrimmedSize(trimmed_width, trimmed_height);

  gridder_ = MakeGridder(trimmed_width, trimmed_height);

  if (TrimWidth() != ImageWidth() || TrimHeight() != ImageHeight()) {
    Image untrimmed_image(ImageWidth(), ImageHeight());
    Logger::Debug << "Untrimming " + std::to_string(TrimWidth()) + " x " +
                         std::to_string(TrimHeight()) + " -> " +
                         std::to_string(ImageWidth()) + " x " +
                         std::to_string(ImageHeight()) + "\n";
    Image::Untrim(untrimmed_image.Data(), ImageWidth(), ImageHeight(),
                  images[0].Data(), TrimWidth(), TrimHeight());
    images[0] = std::move(untrimmed_image);
  }

  if (ImageWidth() != ActualInversionWidth() ||
      ImageHeight() != ActualInversionHeight()) {
    Image resampled_image(ActualInversionWidth(), ActualInversionHeight());
    schaapcommon::math::Resampler resampler(
        ImageWidth(), ImageHeight(), ActualInversionWidth(),
        ActualInversionHeight(), resources_.NCpus());

    resampler.Resample(images[0].Data(), resampled_image.Data());
    images[0] = std::move(resampled_image);
  }

  gridder_->InitializePrediction(images[0].Data());
  images[0].Reset();
}

void WGriddingMSGridder::FinishPredict() {}

}  // namespace wsclean
