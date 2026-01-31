#include "reorderedmsprovider.h"

#include "operations.h"

#include "msreaders/reorderedmsreader.h"

#include "rowproviders/averagingmsrowprovider.h"
#include "rowproviders/directmsrowprovider.h"
#include "rowproviders/msrowprovider.h"
#include "rowproviders/noisemsrowprovider.h"

#include "../main/progressbar.h"
#include "../main/settings.h"

#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <memory>
#include <vector>

#include <aocommon/logger.h>

#include <schaapcommon/reordering/msselection.h>

#include <casacore/measures/Measures/MEpoch.h>

#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/ms/MeasurementSets/MSDataDescription.h>

using aocommon::Logger;

using casacore::ArrayColumn;
using casacore::MSMainEnums;
using casacore::ScalarColumn;

using schaapcommon::reordering::ChannelRange;
using schaapcommon::reordering::FileWriter;
using schaapcommon::reordering::GetMaxChannels;
using schaapcommon::reordering::GetMetaFilename;
using schaapcommon::reordering::GetPartPrefix;
using schaapcommon::reordering::HandleData;
using schaapcommon::reordering::MetaHeader;
using schaapcommon::reordering::MSSelection;
using schaapcommon::reordering::PartHeader;
using schaapcommon::reordering::StorageManagerType;

/**
 * MAP_NORESERVE is unsuported AND not defined on hurd-i386, so
 * assign it to zero in this case.
 */
#ifndef MAP_NORESERVE
#define MAP_NORESERVE 0
#endif

namespace wsclean {

namespace {

std::map<size_t, std::set<aocommon::PolarizationEnum>>
GetMSPolarizationsPerDataDescId(casacore::MeasurementSet& ms) {
  casacore::MSDataDescription data_description_table = ms.dataDescription();
  std::map<size_t, std::set<aocommon::PolarizationEnum>>
      ms_polarizations_per_data_desc_id;
  for (size_t data_desc_id = 0; data_desc_id != data_description_table.nrow();
       ++data_desc_id) {
    ms_polarizations_per_data_desc_id.emplace(
        data_desc_id, GetMSPolarizations(data_desc_id, ms));
  }
  return ms_polarizations_per_data_desc_id;
}

/**
 * Returns a map that maps data_desc_id to spw index.
 */
std::map<size_t, size_t> GetSpwMap(const aocommon::MultiBandData& all_bands) {
  std::map<size_t, size_t> result;
  for (size_t data_desc_id : all_bands.DataDescIds()) {
    const size_t spw_index = all_bands.GetBandIndex(data_desc_id);
    result.emplace(data_desc_id, spw_index);
  }
  return result;
}

}  // namespace

ReorderedMsProvider::ReorderedMsProvider(
    const ReorderedHandle& handle, size_t part_index,
    aocommon::PolarizationEnum polarization)
    : handle_(handle),
      part_index_(part_index),
      current_output_row_(0),
      polarization_(polarization),
      polarization_count_in_file_(
          aocommon::Polarization::GetVisibilityCount(polarization_)) {
  const size_t meta_file_index = handle_.data_->metadata_indices_[part_index_];
  std::ifstream meta_file(GetMetaFilename(handle.data_->ms_path_,
                                          handle.data_->temporary_directory_,
                                          meta_file_index));
  if (!meta_file) {
    throw std::runtime_error("Error opening meta file for ms " +
                             handle.data_->ms_path_);
  }

  meta_header_.Read(meta_file);
  std::vector<char> ms_path(meta_header_.filename_length + 1, char(0));
  meta_file.read(ms_path.data(), meta_header_.filename_length);
  Logger::Info << "Opening reordered part " << part_index << " for "
               << ms_path.data() << '\n';
  std::string part_prefix =
      GetPartPrefix(ms_path.data(), part_index, polarization,
                    handle.data_->temporary_directory_);

  std::ifstream data_file(part_prefix + ".tmp", std::ios::in);
  if (!data_file.good())
    throw std::runtime_error("Error opening temporary data file '" +
                             part_prefix + ".tmp'");
  part_header_.Read(data_file);
  if (!data_file.good())
    throw std::runtime_error("Error reading header from file '" + part_prefix +
                             ".tmp'");

  if (part_header_.has_model) {
    model_file_path_ = part_prefix + "-m.tmp";
    model_file_ = MappedFile(model_file_path_);
  }
  meta_file.close();
  data_file.close();

  if (!meta_header_.data_desc_id.HasValue()) reader_.emplace(this);
}

ReorderedMsProvider::~ReorderedMsProvider() = default;

std::unique_ptr<MSReader> ReorderedMsProvider::MakeReader() {
  return std::make_unique<ReorderedMsReader>(this);
}

std::string ReorderedMsProvider::PartDescription() const {
  return handle_.data_->ms_path_ + ":" + std::to_string(part_index_);
}

void ReorderedMsProvider::NextOutputRow() {
  ++current_output_row_;
  if (reader_) {
    current_output_position_ +=
        SelectedBands()[reader_->CurrentDataDescId()].ChannelCount() *
        polarization_count_in_file_;
    reader_->NextInputRow();
  }
}

void ReorderedMsProvider::ResetWritePosition() {
  current_output_row_ = 0;
  current_output_position_ = 0;
  if (!meta_header_.data_desc_id.HasValue()) reader_.emplace(this);
}

void ReorderedMsProvider::WriteModel(const std::complex<float>* buffer,
                                     bool add_to_ms) {
#ifndef NDEBUG
  if (!part_header_.has_model)
    throw std::runtime_error("Reordered MS initialized without model");
#endif
  size_t n_channels;
  uint64_t position;
  if (reader_) {
    n_channels = SelectedBands()[reader_->CurrentDataDescId()].ChannelCount();
    position = current_output_position_;
  } else {
    n_channels = part_header_.max_channel_count;
    position = n_channels * polarization_count_in_file_ * current_output_row_;
  }

  std::complex<float>* model_write_ptr =
      reinterpret_cast<std::complex<float>*>(model_file_.Data()) + position;

  // In case the value was not sampled in this pass, it has been set to infinite
  // and should not overwrite the current value in the set.
  if (add_to_ms) {
    for (size_t i = 0; i != n_channels * polarization_count_in_file_; ++i) {
      if (std::isfinite(buffer[i].real())) model_write_ptr[i] += buffer[i];
    }
  } else {
    for (size_t i = 0; i != n_channels * polarization_count_in_file_; ++i) {
      if (std::isfinite(buffer[i].real())) model_write_ptr[i] = buffer[i];
    }
  }
}

void ReorderedMsProvider::ResetModelColumn() {
  model_file_.Close();
  std::remove(model_file_path_.c_str());
  schaapcommon::reordering::AllocateFile(model_file_path_,
                                         NMaxChannels() * NPolarizations() *
                                             meta_header_.selected_row_count *
                                             sizeof(std::complex<float>));
  model_file_ = MappedFile(model_file_path_);
}

/*
 * When reordered:
 * One global file stores:
 * - Metadata:
 *   * Number of selected rows
 *   * Filename length + string
 *   * [ UVW, data_desc_id ]
 * The binary parts store the following information:
 * - Number of channels
 * - Start channel in MS
 * - Total weight in part
 * - Data    (single polarization, as requested)
 * - Weights (single)
 * - Model, optionally
 */
ReorderedHandle ReorderMS(
    const std::string& ms_path,
    const std::vector<aocommon::VectorMap<ChannelRange>>& channels,
    const MSSelection& selection, const std::string& data_column_name,
    const std::string& model_column_name,
    StorageManagerType model_storage_manager, bool include_model,
    bool initial_model_required, const Settings& settings) {
  const bool model_update_required = settings.modelUpdateRequired;
  std::set<aocommon::PolarizationEnum> pols_out;
  for (aocommon::PolarizationEnum p : settings.polarizations)
    pols_out.insert(settings.GetProviderPolarization(p));
  const std::string& temporary_directory = settings.temporaryDirectory;

  const size_t channel_parts = channels.size();

  // This maps data_desc_id to spw index.
  // TODO use VectorMap instead
  std::map<size_t, size_t> selected_data_desc_ids;
  bool has_frequency_bda;
  {  // close ms asap
    casacore::MeasurementSet ms(ms_path);
    selected_data_desc_ids = GetSpwMap(aocommon::MultiBandData(ms));
    has_frequency_bda = wsclean::HasFrequencyBda(ms);
  }

  std::unique_ptr<MsRowProviderBase> row_provider;
  if (settings.baselineDependentAveragingInWavelengths == 0.0) {
    if (settings.simulateNoise) {
      std::unique_ptr<NoiseMSRowProvider> noise_row_provider(
          new NoiseMSRowProvider(ms_path, selection, selected_data_desc_ids,
                                 data_column_name, model_column_name,
                                 initial_model_required));
      if (settings.simulatedBaselineNoiseFilename.empty())
        noise_row_provider->SetNoiseLevel(settings.simulatedNoiseStdDev);
      else
        noise_row_provider->SetNoiseBaselineFile(
            settings.simulatedBaselineNoiseFilename);
      row_provider = std::move(noise_row_provider);
    } else
      row_provider = MakeMsRowProvider(
          ms_path, selection, selected_data_desc_ids, data_column_name,
          model_column_name, initial_model_required);
  } else {
    if (initial_model_required)
      throw std::runtime_error(
          "Baseline-dependent averaging is enabled together with a mode that "
          "requires the model data (e.g. -continue or -subtract-model). This "
          "is not possible.");
    row_provider = std::make_unique<AveragingMSRowProvider>(
        settings.baselineDependentAveragingInWavelengths, ms_path, selection,
        selected_data_desc_ids, settings.fieldIds[0], data_column_name,
        model_column_name, initial_model_required);
  }

  const std::map<size_t, std::set<aocommon::PolarizationEnum>>
      ms_polarizations_per_data_desc_id =
          GetMSPolarizationsPerDataDescId(row_provider->Ms());
  const size_t nAntennas = row_provider->Ms().antenna().nrow();
  const aocommon::MultiBandData original_bands(row_provider->Ms());

  // This handle is just for the writer
  const std::vector<aocommon::MultiBandData> bands_per_part =
      MakeSelectedBands(original_bands, channels);
  auto handle_data = std::make_unique<HandleData>(
      ms_path, data_column_name, model_column_name, model_storage_manager,
      temporary_directory, channels, initial_model_required,
      model_update_required, pols_out, selection, bands_per_part,
      has_frequency_bda, nAntennas, settings.saveReorder,
      ReorderedMsProvider::StoreReorderedInMS);

  std::vector<aocommon::OptionalNumber<size_t>> data_desc_ids;
  std::tie(handle_data->metadata_indices_, data_desc_ids) =
      schaapcommon::reordering::MakeMetaFilesMap(handle_data->channels_);

  Logger::Debug << "Using " << data_desc_ids.size()
                << " temporary metadata files.\n";

  FileWriter reordered_file_writer(*handle_data,
                                   ms_polarizations_per_data_desc_id,
                                   data_desc_ids, row_provider->StartTime());

  if (settings.parallelReordering == 1)
    Logger::Info << "Reordering " << ms_path << " into " << channel_parts
                 << " x " << pols_out.size() << " parts.\n";

  casacore::Array<std::complex<float>> data_array;
  casacore::Array<std::complex<float>> model_array;
  casacore::Array<float> weight_spectrum_array;
  casacore::Array<bool> flag_array;

  std::unique_ptr<ProgressBar> progress1;
  if (settings.parallelReordering == 1)
    progress1.reset(new ProgressBar("Reordering"));

  size_t selected_rows_total = 0;
  aocommon::UVector<size_t> selected_row_count_per_spw_index(
      selected_data_desc_ids.size(), 0);
  while (!row_provider->AtEnd()) {
    if (progress1)
      progress1->SetProgress(row_provider->CurrentProgress(),
                             row_provider->TotalProgress());

    double time, u, v, w;
    uint32_t data_desc_id, antenna1, antenna2, field_id;
    row_provider->ReadData(data_array, flag_array, weight_spectrum_array, u, v,
                           w, data_desc_id, antenna1, antenna2, field_id, time);

    if (initial_model_required) row_provider->ReadModel(model_array);

    reordered_file_writer.WriteMetaRow(u, v, w, time, data_desc_id, antenna1,
                                       antenna2, field_id);
    reordered_file_writer.WriteDataRow(data_array.data(), model_array.data(),
                                       weight_spectrum_array.data(),
                                       flag_array.data(), data_desc_id);

    row_provider->NextRow();
    ++selected_rows_total;
  }
  progress1.reset();
  Logger::Debug << "Total selected rows: " << selected_rows_total << '\n';
  row_provider->OutputStatistics();

  reordered_file_writer.UpdateMetaHeaders();
  reordered_file_writer.UpdatePartHeaders(include_model);

  std::unique_ptr<ProgressBar> progress2;
  if (include_model && !initial_model_required) {
    if (settings.parallelReordering) {
      progress2 =
          std::make_unique<ProgressBar>("Initializing model visibilities");
    }
    auto update_progress = [progress2 = std::move(progress2)](size_t progress,
                                                              size_t total) {
      if (progress2) {
        progress2->SetProgress(progress, total);
      }
    };
    constexpr bool zero_initialize_model = true;
    reordered_file_writer.PopulateModel(zero_initialize_model,
                                        std::cref(update_progress));
  }
  progress2.reset();

  return ReorderedHandle(std::move(handle_data));
}  // namespace wsclean

void ReorderedMsProvider::StoreReorderedInMS(const HandleData& handle) {
  const std::set<aocommon::PolarizationEnum> pols = handle.polarizations_;

  std::ifstream first_data_file(GetPartPrefix(handle.ms_path_, 0, *pols.begin(),
                                              handle.temporary_directory_) +
                                    ".tmp",
                                std::ios::in);
  if (!first_data_file.good())
    throw std::runtime_error("Error opening temporary data file");
  PartHeader firstpart_header_;
  firstpart_header_.Read(first_data_file);
  if (!first_data_file.good())
    throw std::runtime_error("Error reading from temporary data file");

  if (firstpart_header_.has_model) {
    const size_t channel_parts = handle.channels_.size();

    // Open the temporary files
    std::vector<std::unique_ptr<std::ifstream>> model_files(channel_parts *
                                                            pols.size());
    size_t file_index = 0;
    for (size_t part = 0; part != channel_parts; ++part) {
      for (aocommon::PolarizationEnum p : pols) {
        std::string part_prefix = GetPartPrefix(handle.ms_path_, part, p,
                                                handle.temporary_directory_);
        model_files[file_index] =
            std::make_unique<std::ifstream>(part_prefix + "-m.tmp");
        if (!*model_files[file_index])
          throw std::runtime_error("Error opening temporary model data file '" +
                                   part_prefix + "-m.tmp' for reading");
        ++file_index;
      }
    }

    casacore::MeasurementSet ms(handle.ms_path_, casacore::Table::Update);
    const aocommon::VectorMap<std::set<aocommon::PolarizationEnum>>
        ms_polarizations_per_data_desc_id(GetMSPolarizationsPerDataDescId(ms));
    InitializeModelColumn(ms, handle.data_column_name_,
                          handle.model_column_name_,
                          handle.model_storage_manager_);
    ScalarColumn<int> antenna1_column(ms, ms.columnName(MSMainEnums::ANTENNA1));
    ScalarColumn<int> antenna2_column(ms, ms.columnName(MSMainEnums::ANTENNA2));
    ScalarColumn<int> field_id_column(ms, ms.columnName(MSMainEnums::FIELD_ID));
    ScalarColumn<double> time_column(ms, ms.columnName(MSMainEnums::TIME));
    ScalarColumn<int> data_desc_id_column(
        ms, ms.columnName(MSMainEnums::DATA_DESC_ID));
    ArrayColumn<casacore::Complex> data_column(ms, handle.data_column_name_);
    ArrayColumn<casacore::Complex> model_column(ms, handle.model_column_name_);
    ArrayColumn<double> uvw_column(ms, ms.columnName(MSMainEnums::UVW));

    casacore::IPosition shape(data_column.shape(0));
    const size_t maxchannels_ = GetMaxChannels(handle.channels_);

    const size_t polarizations_per_file =
        aocommon::Polarization::GetVisibilityCount(*pols.begin());
    std::vector<std::complex<float>> model_data_buffer(maxchannels_ *
                                                       polarizations_per_file);
    casacore::Array<std::complex<float>> model_data_array(shape);

    ProgressBar progress(std::string("Writing changed model back to ") +
                         handle.ms_path_);
    size_t start_row, end_row;
    GetRowRange(ms, handle.selection_, start_row, end_row);
    size_t timestep =
        handle.selection_.HasInterval() ? handle.selection_.IntervalStart() : 0;
    double time = time_column(start_row);
    size_t selected_row_count_for_debug = 0;
    for (size_t row = start_row; row != end_row; ++row) {
      progress.SetProgress(row - start_row, end_row - start_row);
      const int antenna1 = antenna1_column(row);
      const int antenna2 = antenna2_column(row);
      const int field_id = field_id_column(row);
      const size_t data_desc_id = data_desc_id_column(row);
      casacore::Vector<double> uvw = uvw_column(row);

      if (time != time_column(row)) {
        ++timestep;
        time = time_column(row);
      }
      if (handle.selection_.IsSelected(field_id, timestep, antenna1, antenna2,
                                       uvw.data())) {
        if (handle.model_storage_manager_ == StorageManagerType::Sisco) {
          shape = data_column.shape(row);
          model_data_array.resize(shape);
        } else {
          model_column.get(row, model_data_array, true);
        }
        size_t file_index = 0;
        for (size_t part = 0; part != channel_parts; ++part) {
          const ChannelRange& range = handle.channels_[part][data_desc_id];
          if (!range.Empty()) {
            const size_t part_start_ch = range.start;
            const size_t part_end_ch = range.end;
            const std::set<aocommon::PolarizationEnum>& ms_polarizations =
                ms_polarizations_per_data_desc_id[data_desc_id];
            for (aocommon::PolarizationEnum p : pols) {
              model_files[file_index]->read(
                  reinterpret_cast<char*>(model_data_buffer.data()),
                  (part_end_ch - part_start_ch) * polarizations_per_file *
                      sizeof(std::complex<float>));
              if (!model_files[file_index]->good())
                throw std::runtime_error(
                    "Error reading from temporary model data file");
              schaapcommon::reordering::StoreData<false>(
                  model_data_array.data(), part_start_ch, part_end_ch,
                  ms_polarizations, model_data_buffer.data(), p);

              ++file_index;
            }
          } else {
            file_index += pols.size();
          }
        }
        model_column.put(row, model_data_array);
        selected_row_count_for_debug++;
      } else if (handle.model_storage_manager_ == StorageManagerType::Sisco) {
        // When Sisco is used, the shape of every cell should be set, even
        // unselected ones.
        shape = data_column.shape(row);
        model_data_array.resize(shape);
        model_data_array = 0.0f;
        model_column.put(row, model_data_array);
      }
    }
    progress.SetProgress(end_row - start_row, end_row - start_row);

    Logger::Debug
        << "Row count for writing reordered data back to the measurement set: "
        << selected_row_count_for_debug << '\n';
  }
}

double ReorderedMsProvider::Interval() { return GetMsInterval(*MS()); }

ObservationInfo ReorderedMsProvider::GetObservationInfo() {
  return ReadObservationInfo(*MS(),
                             handle_.data_->selection_.FieldIds().front());
}

std::vector<std::string> ReorderedMsProvider::GetAntennaNames() {
  return wsclean::GetAntennaNames(MS()->antenna());
}

}  // namespace wsclean
