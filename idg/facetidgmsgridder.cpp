#include "facetidgmsgridder.h"

#include "../msproviders/msreaders/timestepbufferreader.h"

#include <cmath>
#include <thread>

#include <idg-api.h>

#include <aocommon/coordinatesystem.h>
#include <aocommon/fits/fitsreader.h>
#include <aocommon/logger.h>

#include "../msproviders/msprovider.h"
#include "../msproviders/timestepbuffer.h"

#include "../io/findmwacoefffile.h"
#include "../io/imagefilename.h"
#include "../io/parsetreader.h"

#include "../structures/imagingtable.h"

#include "../main/settings.h"

#include "idgconfiguration.h"

using aocommon::CoordinateSystem;
using aocommon::Image;
using aocommon::Logger;

namespace wsclean {

namespace {
constexpr const size_t kGridderIndex = 0;
}

FacetIdgMsGridder::FacetIdgMsGridder(
    const Settings& settings, const Resources& resources,
    MsProviderCollection& ms_provider_collection)
    : MsGridder(settings, ms_provider_collection),
      output_provider_(nullptr),
      proxy_type_(idg::api::Type::CPU_OPTIMIZED),
      buffer_size_(0),
      resources_(resources) {
  IdgConfiguration::Read(proxy_type_, buffer_size_, options_);

  proxy_type_ = GetIdgType(GetSettings());

  buffer_set_ = std::unique_ptr<idg::api::BufferSet>(
      idg::api::BufferSet::create(proxy_type_));
  options_["max_threads"] = int(resources.NCpus());

  if (settings.gridMode == GriddingKernelMode::BlackmanHarris) {
    options_["taper"] = std::string("blackman-harris");
  }

  // This MsGridder is only meant to be used in combination with faceting,
  // hence, if multiple polarizations are requested, they will be gridded
  // separately. Thus, as far as IDG is concerned, we're only gridding Stokes I.
  options_["stokes_I_only"] = true;
}

void FacetIdgMsGridder::StartInversion() {
  const size_t untrimmed_width = ImageWidth();
  assert(TrimWidth() == TrimHeight());
  assert(untrimmed_width == ImageHeight());

  options_["padded_size"] = untrimmed_width;

  double max_w = 0;
  for (size_t i = 0; i != GetMsCount(); ++i) {
    max_w = std::max(max_w, GetMsData(i).max_w_with_flags);
  }

  const double shift_l = LShift();
  const double shift_m = MShift();
  buffer_set_->init(TrimWidth(), ActualPixelSizeX(), max_w + 1.0, shift_l,
                    shift_m, options_);
  Logger::Debug << "IDG subgrid size: " << buffer_set_->get_subgridsize()
                << '\n';
}

void FacetIdgMsGridder::FinishInversion() {
  // GridMeasurementSet calls have added the gridding result to image_ member
  image_.assign(TrimWidth() * TrimHeight(), 0.0);
  buffer_set_->get_image(image_.data());
  if (GetPsfMode() != PsfMode::kNone) {
    Logger::Debug << "Total weight: " << ImageWeight() << '\n';
  }
  // result is now in image_ member
  // Can be accessed by subsequent calls to ResultImages()
}

size_t FacetIdgMsGridder::GridMeasurementSet(
    const MsProviderCollection::MsData& ms_data) {
  if (!PrepareForMeasurementSet(ms_data, idg::api::BufferSetType::gridding)) {
    return 0;
  }

  const aocommon::BandData& selected_band = *selected_bands_.begin();

  // Even though facet-idg will only grid a single polarization, the data and
  // weights buffers are sent to IDG, which always expects four polarizations.
  constexpr size_t n_idg_polarizations = 4;
  const size_t data_size = selected_band.ChannelCount() * n_idg_polarizations;
  aocommon::UVector<float> weight_buffer(data_size);
  aocommon::UVector<bool> selection_buffer(selected_band.ChannelCount(), true);

  // Since the model data buffer is only used internal to WSClean, it can be
  // smaller.
  const size_t model_data_size =
      selected_band.ChannelCount() * ms_data.ms_provider->NPolarizations();
  aocommon::UVector<std::complex<float>> model_buffer(model_data_size);

  // The gridder doesn't need to know the absolute time index; this value
  // indexes relatively to where we start in the measurement set, and only
  // increases when the time changes.
  int time_index = -1;
  double current_time = -1.0;
  aocommon::UVector<double> uvws(ms_data.ms_provider->NAntennas() * 3, 0.0);

  TimestepBuffer timestep_buffer(ms_data.ms_provider, DoSubtractModel());
  std::unique_ptr<MSReader> ms_reader = timestep_buffer.MakeReader();
  TimestepBufferReader& timestep_reader =
      static_cast<TimestepBufferReader&>(*ms_reader);
  aocommon::UVector<std::complex<float>> row_visibilities(data_size);
  IDGInversionRow row_data;
  row_data.data = row_visibilities.data();

  while (ms_reader->CurrentRowAvailable()) {
    MSProvider::MetaData metadata;
    timestep_reader.ReadMeta(metadata);

    if (current_time != metadata.time) {
      current_time = metadata.time;
      time_index++;
    }

    row_data.uvw[0] = metadata.u_in_m;
    row_data.uvw[1] = metadata.v_in_m;
    row_data.uvw[2] = metadata.w_in_m;

    row_data.antenna1 = metadata.antenna1;
    row_data.antenna2 = metadata.antenna2;
    row_data.time_index = time_index;

    const size_t n_parms = NumValuesPerSolution();
    if (n_parms == 2) {
      GetCollapsedVisibilities<2>(*ms_reader, ms_data.antenna_names.size(),
                                  row_data, weight_buffer.data(),
                                  model_buffer.data(), selection_buffer.data(),
                                  metadata);
    } else {
      GetCollapsedVisibilities<4>(*ms_reader, ms_data.antenna_names.size(),
                                  row_data, weight_buffer.data(),
                                  model_buffer.data(), selection_buffer.data(),
                                  metadata);
    }

    size_t source_index = row_visibilities.size() / 4;
    for (size_t i = row_visibilities.size(); i != 0; i -= 4) {
      row_visibilities[i - 1] = row_visibilities[source_index - 1];
      row_visibilities[i - 2] = 0.0;
      row_visibilities[i - 3] = 0.0;
      row_visibilities[i - 4] = row_visibilities[source_index - 1];
      weight_buffer[i - 1] = weight_buffer[source_index - 1];
      weight_buffer[i - 2] = weight_buffer[source_index - 1];
      weight_buffer[i - 3] = weight_buffer[source_index - 1];
      weight_buffer[i - 4] = weight_buffer[source_index - 1];
      source_index--;
    }

    row_data.uvw[1] = -metadata.v_in_m;  // DEBUG vdtol, flip axis
    row_data.uvw[2] = -metadata.w_in_m;  //

    buffer_set_->get_gridder(kGridderIndex)
        ->grid_visibilities(time_index, metadata.antenna1, metadata.antenna2,
                            row_data.uvw, row_data.data, weight_buffer.data());

    ms_reader->NextInputRow();
  }
  buffer_set_->finished();

  return 0;
}

void FacetIdgMsGridder::StartPredict(std::vector<Image>&& images) {
  if (images.size() == 2)
    throw std::runtime_error("IDG gridder cannot make complex images");
  const size_t untrimmed_width = ImageWidth();
  const size_t width = TrimWidth();
  const size_t height = TrimHeight();

  assert(width == height);
  assert(untrimmed_width == ImageHeight());

  options_["padded_size"] = untrimmed_width;

  // Since we only (de)grid a single polarization at a time in when using the
  // faceting engine, the number of image polarizations will always be one.
  constexpr size_t n_image_polarizations = 1;
  image_.assign(n_image_polarizations * width * height, 0.0);

  assert(images.size() == n_image_polarizations);
  std::copy_n(images[0].Data(), width * height, image_.data());

  double max_w = 0;
  for (size_t i = 0; i != GetMsCount(); ++i) {
    max_w = std::max(max_w, GetMsData(i).max_w_with_flags);
  }

  const double shift_l = LShift();
  const double shift_m = MShift();
  buffer_set_->init(width, ActualPixelSizeX(), max_w + 1.0, shift_l, shift_m,
                    options_);

  // FacetIdgMsGridder doesn't require scaling since we don't apply the average
  // beam.
  constexpr bool do_scale = false;
  buffer_set_->set_image(image_.data(), do_scale);
}

size_t FacetIdgMsGridder::PredictMeasurementSet(
    const MsProviderCollection::MsData& ms_data) {
  if (!PrepareForMeasurementSet(ms_data, idg::api::BufferSetType::degridding)) {
    return 0;
  }

  ms_data.ms_provider->ReopenRW();

  output_provider_ = ms_data.ms_provider;

  const aocommon::BandData& selected_band = *selected_bands_.begin();

  constexpr size_t n_idg_polarizations = 4;
  aocommon::UVector<std::complex<float>> buffer(selected_band.ChannelCount() *
                                                n_idg_polarizations);

  int time_index = -1;
  double current_time = -1.0;
  aocommon::UVector<double> uvws(ms_data.ms_provider->NAntennas() * 3, 0.0);

  TimestepBuffer timestep_buffer(ms_data.ms_provider, false);
  timestep_buffer.ResetWritePosition();
  for (std::unique_ptr<MSReader> ms_reader = timestep_buffer.MakeReader();
       ms_reader->CurrentRowAvailable(); ms_reader->NextInputRow()) {
    TimestepBufferReader& timestep_reader =
        static_cast<TimestepBufferReader&>(*ms_reader);

    MSProvider::MetaData metadata;
    timestep_reader.ReadMeta(metadata);

    const size_t provider_row_id = timestep_reader.RowId();
    if (current_time != metadata.time) {
      current_time = metadata.time;
      time_index++;
    }

    IDGPredictionRow row;
    row.uvw[0] = metadata.u_in_m;
    row.uvw[1] = -metadata.v_in_m;
    row.uvw[2] = -metadata.w_in_m;
    row.antenna1 = metadata.antenna1;
    row.antenna2 = metadata.antenna2;
    row.time_index = time_index;
    row.row_id = provider_row_id;
    PredictRow(row, ms_data.antenna_names);
  }

  ComputePredictionBuffer(ms_data.antenna_names);
  return 0;
}

void FacetIdgMsGridder::PredictRow(
    IDGPredictionRow& row, const std::vector<std::string>& antenna_names) {
  while (buffer_set_->get_degridder(kGridderIndex)
             ->request_visibilities(row.row_id, row.time_index, row.antenna1,
                                    row.antenna2, row.uvw)) {
    ComputePredictionBuffer(antenna_names);
  }
}

void FacetIdgMsGridder::ComputePredictionBuffer(
    const std::vector<std::string>& antenna_names) {
  auto available_row_ids = buffer_set_->get_degridder(kGridderIndex)->compute();
  Logger::Debug << "Computed " << available_row_ids.size() << " rows.\n";
  for (std::pair<long unsigned, std::complex<float>*>& row :
       available_row_ids) {
    MSProvider::MetaData metadata;
    ReadPredictMetaData(metadata);
    const aocommon::BandData& band = selected_bands_[metadata.data_desc_id];
    // Place the single polarization in the first quarter of the array
    for (size_t i = 0; i != band.ChannelCount(); ++i) {
      row.second[i] = (row.second[i * 4] + row.second[i * 4 + 3]) / 2.0f;
    }
    const double* uvw = nullptr;
    WriteCollapsedVisibilities(*output_provider_, antenna_names.size(),
                               metadata.data_desc_id, row.second, uvw,
                               metadata.field_id, metadata.antenna1,
                               metadata.antenna2, metadata.time);
  }
  buffer_set_->get_degridder(kGridderIndex)->finished_reading();
}

std::vector<Image> FacetIdgMsGridder::ResultImages() {
  const size_t width = TrimWidth();
  const size_t height = TrimHeight();
  std::vector<Image> images;

  images.emplace_back(width, height);
  std::copy_n(image_.data(), width * height, images[0].Data());

  return images;
}

bool FacetIdgMsGridder::PrepareForMeasurementSet(
    const MsProviderCollection::MsData& ms_data,
    idg::api::BufferSetType bufferSetType) {
  const float max_baseline = ms_data.max_baseline_in_m;
  // Skip this ms if there is no data in it
  if (!max_baseline) return false;

  selected_bands_ = ms_data.ms_provider->SelectedBands();

  // TODO For now we do not allow multiple bands in one msprovider
  if (!ms_data.ms_provider->IsRegular())
    throw std::runtime_error(
        "IDG implementation does not support irregular data yet");
  // Regular data should always have one band...
  assert(selected_bands_.BandCount() == 1);
  const aocommon::BandData& selected_band = *selected_bands_.begin();

  // TODO for now we map the ms antennas directly to the gridder's antenna,
  // including non-selected antennas. Later this can be made more efficient.
  const size_t nStations = ms_data.ms_provider->NAntennas();

  std::vector<std::vector<double>> bands;
  bands.emplace_back(selected_band.begin(), selected_band.end());
  const size_t nChannels = selected_band.ChannelCount();

  // Only one-third of the mem is allocated to the buffers, so that memory
  // remains available for the images and other things done by IDG.
  // Never use more than 16 GB
  const size_t memSize = std::min<uint64_t>(16ul * 1024ul * 1024ul * 1024ul,
                                            resources_.Memory() / 3);
  uint64_t memPerTimestep =
      idg::api::BufferSet::get_memory_per_timestep(nStations, nChannels);

  // IDG can allocate two visibility buffers: (for parallel processing)
  memPerTimestep *= 2;

  buffer_size_ = std::max<size_t>(1, memSize / memPerTimestep);

  Logger::Debug << "Allocatable timesteps (" << nStations << " stations, "
                << nChannels << " channels, " << memSize / (1024 * 1024 * 1024)
                << " GB mem): " << buffer_size_ << '\n';
  buffer_set_->init_buffers(buffer_size_, bands, nStations, max_baseline,
                            options_, bufferSetType);

  return true;
}

}  // namespace wsclean
