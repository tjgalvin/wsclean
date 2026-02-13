#include "wtowersmsgridder.h"

#include <aocommon/image.h>
#include <aocommon/logger.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <schaapcommon/math/resampler.h>

#include "wtowers_gridder.h"
#include "wtowers_gridder_implementation.h"
#include "../gridding/msgriddermanager.h"
#include "../msproviders/msreaders/msreader.h"
#include "../msproviders/msprovider.h"
#include "../structures/imageweights.h"

template class wsclean::WTowersGridder<double>;
template class wsclean::WTowersGridder<float>;

using aocommon::Image;
using aocommon::Logger;

namespace wsclean {

WTowersMsGridder::WTowersMsGridder(const Settings& settings,
                                   const Resources& resources,
                                   MsProviderCollection& ms_provider_collection)
    : MsGridder(settings, ms_provider_collection),
      resources_(resources),
      wtowers_subgrid_size_(GetSettings().wtowers_subgrid_size),
      wtowers_support_(GetSettings().wtowers_support),
      wtowers_w_support_(GetSettings().wtowers_w_support),
      wtowers_padding_(GetSettings().wtowers_padding),
      wtowers_w_padding_(GetSettings().wtowers_w_padding),
      accuracy_(GetSettings().gridder_accuracy) {
  // It may happen that several schaapcommon::fft::Resamplers are created
  // concurrently, so we must make sure that the FFTW planner can deal with
  // this.
  fftwf_make_planner_thread_safe();
}

WTowersMsGridder::~WTowersMsGridder() = default;

std::unique_ptr<WTowersGridderBase> WTowersMsGridder::MakeGridder(
    size_t width, size_t height) const {
  if (accuracy_ <= 1.01e-5) {
    return std::make_unique<WTowersGridder<double>>(
        ActualInversionWidth(), ActualInversionHeight(), width, height,
        ActualPixelSizeX(), ActualPixelSizeY(), LShift(), MShift(),
        wtowers_subgrid_size_, wtowers_support_, wtowers_w_support_,
        wtowers_padding_, wtowers_w_padding_, resources_.NCpus(), accuracy_,
        max_abs_w_);
  } else {
    return std::make_unique<WTowersGridder<float>>(
        ActualInversionWidth(), ActualInversionHeight(), width, height,
        ActualPixelSizeX(), ActualPixelSizeY(), LShift(), MShift(),
        wtowers_subgrid_size_, wtowers_support_, wtowers_w_support_,
        wtowers_padding_, wtowers_w_padding_, resources_.NCpus(), accuracy_,
        max_abs_w_);
  }
}

size_t WTowersMsGridder::CalculateConstantMemory() const {
  size_t constant_mem = gridder_->ConstantMemoryUsage();
  constant_mem += GetVisibilityModifier().GetCacheParmResponseSize();
  return constant_mem;
}

size_t WTowersMsGridder::CalculateMaxRowsInMemory(
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

size_t WTowersMsGridder::GridMeasurementSet(
    const MsProviderCollection::MsData& ms_data) {
  const size_t n_vis_polarizations = ms_data.ms_provider->NPolarizations();

  // TODO For now we do not allow multiple bands in one msprovider
  const aocommon::MultiBandData& selected_bands(
      ms_data.ms_provider->SelectedBands());
  if (!ms_data.ms_provider->IsRegular())
    throw std::runtime_error(
        "w-towers implementation does not support irregular data yet");
  // Regular data should always have one band...
  assert(selected_bands.BandCount() == 1);
  const aocommon::BandData& selected_band = *selected_bands.begin();

  const size_t data_size = selected_band.ChannelCount() * n_vis_polarizations;
  aocommon::UVector<std::complex<float>> model_buffer(data_size);
  aocommon::UVector<float> weight_buffer(data_size);
  aocommon::UVector<bool> selection_buffer(selected_band.ChannelCount(), true);

  aocommon::UVector<double> frequencies(selected_band.ChannelCount());
  for (size_t i = 0; i != frequencies.size(); ++i)
    frequencies[i] = selected_band.ChannelFrequency(i);

  const size_t per_row_uvw_memory_consumption = sizeof(double) * 3;
  size_t max_rows_per_chunk = CalculateMaxRowsInMemory(
      resources_.Memory(), CalculateConstantMemory(), 0,
      per_row_uvw_memory_consumption, selected_band.ChannelCount(), 1);

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
  while (ms_reader->CurrentRowAvailable()) {
    Logger::Debug << "Max " << max_rows_per_chunk << " rows fit in memory.\n";
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

      // Negate v; otherwise image is flipped compared to other gridders like
      // DUCC, that also do this. We must do this after GetCollapsedVisibilities
      // call otherwise we corrupt the image.
      uvw_buffer[(n_chunk_rows_read * 3) + 1] *= -1;
      // Negate w; otherwise we get distortion as we move away from the centre
      uvw_buffer[(n_chunk_rows_read * 3) + 2] *= -1;

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

size_t WTowersMsGridder::PredictMeasurementSet(
    const MsProviderCollection::MsData& ms_data) {
  ms_data.ms_provider->ReopenRW();

  // TODO For now we do not allow multiple bands in one msprovider
  const aocommon::MultiBandData& selected_bands(
      ms_data.ms_provider->SelectedBands());
  if (!ms_data.ms_provider->IsRegular())
    throw std::runtime_error(
        "w-towers implementation does not support irregular data yet");
  // Regular data should always have one band...
  assert(selected_bands.BandCount() == 1);
  const aocommon::BandData& selected_band = *selected_bands.begin();
  const size_t data_desc_id = *selected_bands.DataDescIds().begin();

  size_t n_total_rows_read = 0;

  aocommon::UVector<double> frequencies(selected_band.ChannelCount());
  for (size_t i = 0; i != frequencies.size(); ++i)
    frequencies[i] = selected_band.ChannelFrequency(i);

  const size_t per_row_uvw_memory_consumption = sizeof(double) * 3;
  size_t max_rows_per_chunk = CalculateMaxRowsInMemory(
      resources_.Memory(), CalculateConstantMemory(), 0,
      per_row_uvw_memory_consumption, selected_band.ChannelCount(), 1);

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
      // Negate v; otherwise image is flipped compared to other gridders like
      // DUCC that also do this
      uvw_buffer[n_chunk_rows_read * 3 + 1] = -metadata.v_in_m;
      // Negate w; otherwise we get distortion as we move away from the centre
      uvw_buffer[n_chunk_rows_read * 3 + 2] = -metadata.w_in_m;
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
          *ms_data.ms_provider, ms_data.antenna_names.size(), data_desc_id,
          &visibility_buffer[row * selected_band.ChannelCount()],
          &uvw_buffer[row * 3], metadata_buffer[row].field_id,
          metadata_buffer[row].antenna1, metadata_buffer[row].antenna2,
          metadata_buffer[row].time);
    }
    n_total_rows_read += n_chunk_rows_read;
  }  // end of chunk
  return n_total_rows_read;
}

void WTowersMsGridder::GetActualTrimmedSize(size_t& trimmed_width,
                                            size_t& trimmed_height) const {
  trimmed_width = std::ceil(ActualInversionWidth() / ImagePadding());
  trimmed_height = std::ceil(ActualInversionHeight() / ImagePadding());

  // In facet-based imaging, the alignment is 4, see wsclean.cpp. Also for
  // monolithic imaging - in which just an even number would suffice -
  // the trimmed_width and trimmed_height are defined to be divisable by 4.
  const size_t alignment = 4;
  if (trimmed_width % alignment != 0) {
    trimmed_width += alignment - (trimmed_width % alignment);
  }
  if (trimmed_height % alignment != 0) {
    trimmed_height += alignment - (trimmed_height % alignment);
  }
  trimmed_width = std::min(trimmed_width, ActualInversionWidth());
  trimmed_height = std::min(trimmed_height, ActualInversionHeight());
}

void WTowersMsGridder::CalculateGridderMetaData() {
  max_abs_w_ = 0;
  for (size_t ms_index = 0; ms_index < GetMsCount(); ++ms_index) {
    MsProviderCollection::MsData& ms_data = GetMsData(ms_index);
    max_abs_w_ = std::max(max_abs_w_, -ms_data.min_w);
    max_abs_w_ = std::max(max_abs_w_, ms_data.max_w);
  }
}

void WTowersMsGridder::StartInversion() {
  size_t trimmed_width;
  size_t trimmed_height;
  GetActualTrimmedSize(trimmed_width, trimmed_height);

  gridder_ = MakeGridder(trimmed_width, trimmed_height);
  gridder_->InitializeInversion();

  ResetVisibilityCounters();
}

void WTowersMsGridder::FinishInversion() {
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
    Logger::Debug << "Trimming " << ImageWidth() << " x " << ImageHeight()
                  << " -> " << TrimWidth() << " x " << TrimHeight() << '\n';

    image_ = image_.Trim(TrimWidth(), TrimHeight());
  }
}

void WTowersMsGridder::StartPredict(std::vector<Image>&& images) {
  size_t trimmed_width;
  size_t trimmed_height;
  GetActualTrimmedSize(trimmed_width, trimmed_height);

  gridder_ = MakeGridder(trimmed_width, trimmed_height);

  if (TrimWidth() != ImageWidth() || TrimHeight() != ImageHeight()) {
    Image untrimmed_image(ImageWidth(), ImageHeight());
    Logger::Debug << "Untrimming " << TrimWidth() << " x " << TrimHeight()
                  << " -> " << ImageWidth() << " x " << ImageHeight() << '\n';
    Image::Untrim(untrimmed_image.Data(), ImageWidth(), ImageHeight(),
                  images[0].Data(), TrimWidth(), TrimHeight());
    images[0] = std::move(untrimmed_image);
  }

  if (ImageWidth() != ActualInversionWidth() ||
      ImageHeight() != ActualInversionHeight()) {
    Image resampled_image(ImageWidth(), ImageHeight());
    schaapcommon::math::Resampler resampler(
        ImageWidth(), ImageHeight(), ActualInversionWidth(),
        ActualInversionHeight(), resources_.NCpus());

    resampler.Resample(images[0].Data(), resampled_image.Data());
    images[0] = std::move(resampled_image);
  }

  gridder_->InitializePrediction(images[0].Data());
  images[0].Reset();
}

void WTowersMsGridder::FinishPredict() {}

}  // namespace wsclean
