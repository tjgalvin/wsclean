#include "reorderedmsprovider.h"
#include "msreaders/reorderedmsreader.h"

#include "averagingmsrowprovider.h"
#include "directmsrowprovider.h"
#include "msrowprovider.h"
#include "noisemsrowprovider.h"

#include "../main/progressbar.h"
#include "../main/settings.h"

#include <fstream>
#include <sstream>
#include <memory>
#include <vector>

#include <aocommon/logger.h>

#include <casacore/measures/Measures/MEpoch.h>

#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>

using aocommon::Logger;

using schaapcommon::reordering::ChannelRange;
using schaapcommon::reordering::GetDataDescIdMap;
using schaapcommon::reordering::GetMaxChannels;
using schaapcommon::reordering::GetMetaFilename;
using schaapcommon::reordering::GetPartPrefix;
using schaapcommon::reordering::MetaHeader;
using schaapcommon::reordering::MSSelection;
using schaapcommon::reordering::PartHeader;
using schaapcommon::reordering::ReorderedFileWriter;
using schaapcommon::reordering::ReorderedHandleData;
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
GetMSPolarizationsPerDataDescId(const std::vector<ChannelRange>& ranges,
                                casacore::MeasurementSet& ms) {
  std::map<size_t, std::set<aocommon::PolarizationEnum>>
      ms_polarizations_per_data_desc_id;
  for (const ChannelRange& range : ranges) {
    ms_polarizations_per_data_desc_id.emplace(
        range.data_desc_id,
        ReorderedMsProvider::GetMSPolarizations(range.data_desc_id, ms));
  }
  return ms_polarizations_per_data_desc_id;
}

}  // namespace

ReorderedMsProvider::ReorderedMsProvider(
    const ReorderedHandle& handle, size_t part_index,
    aocommon::PolarizationEnum polarization, size_t data_desc_id)
    : handle_(handle),
      part_index_(part_index),
      data_desc_id_(data_desc_id),
      current_output_row_(0),
      polarization_(polarization),
      polarization_count_in_file_(
          aocommon::Polarization::GetVisibilityCount(polarization_)) {
  std::ifstream meta_file(GetMetaFilename(handle.data_->ms_path_,
                                          handle.data_->temporary_directory_,
                                          data_desc_id));
  if (!meta_file) {
    throw std::runtime_error("Error opening meta file for ms " +
                             handle.data_->ms_path_ + ", data_desc_id " +
                             std::to_string(data_desc_id));
  }

  meta_header_.Read(meta_file);
  std::vector<char> ms_path(meta_header_.filename_length + 1, char(0));
  meta_file.read(ms_path.data(), meta_header_.filename_length);
  Logger::Info << "Opening reordered part " << part_index << " spw "
               << data_desc_id << " for " << ms_path.data() << '\n';
  std::string part_prefix =
      GetPartPrefix(ms_path.data(), part_index, polarization, data_desc_id,
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
    const size_t length =
        part_header_.channel_count * meta_header_.selected_row_count *
        polarization_count_in_file_ * sizeof(std::complex<float>);
    model_file_ = MappedFile(part_prefix + "-m.tmp", length);
  }
  meta_file.close();
  data_file.close();
}

ReorderedMsProvider::~ReorderedMsProvider() {}

std::unique_ptr<MSReader> ReorderedMsProvider::MakeReader() {
  std::unique_ptr<MSReader> reader(new ReorderedMsReader(this));
  return reader;
}

void ReorderedMsProvider::NextOutputRow() { ++current_output_row_; }

void ReorderedMsProvider::WriteModel(const std::complex<float>* buffer,
                                     bool add_to_ms) {
#ifndef NDEBUG
  if (!part_header_.has_model)
    throw std::runtime_error("Reordered MS initialized without model");
#endif
  size_t row_length = part_header_.channel_count * polarization_count_in_file_ *
                      sizeof(std::complex<float>);
  std::complex<float>* model_write_ptr = reinterpret_cast<std::complex<float>*>(
      model_file_.Data() + row_length * current_output_row_);

  // In case the value was not sampled in this pass, it has been set to infinite
  // and should not overwrite the current value in the set.
  if (add_to_ms) {
    for (size_t i = 0;
         i != part_header_.channel_count * polarization_count_in_file_; ++i) {
      if (std::isfinite(buffer[i].real())) model_write_ptr[i] += buffer[i];
    }
  } else {
    for (size_t i = 0;
         i != part_header_.channel_count * polarization_count_in_file_; ++i) {
      if (std::isfinite(buffer[i].real())) model_write_ptr[i] = buffer[i];
    }
  }
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
ReorderedMsProvider::ReorderedHandle ReorderMS(
    const std::string& ms_path, const std::vector<ChannelRange>& channels,
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
  const std::map<size_t, size_t> selected_data_desc_ids =
      GetDataDescIdMap(channels);

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
          GetMSPolarizationsPerDataDescId(channels, row_provider->Ms());
  const size_t nAntennas = row_provider->Ms().antenna().nrow();
  const aocommon::MultiBandData bands(row_provider->Ms());

  // This handle is just for the writer
  ReorderedHandleData handle_data(
      ms_path, data_column_name, model_column_name, model_storage_manager,
      temporary_directory, channels, initial_model_required,
      model_update_required, pols_out, selection, bands, nAntennas, true,
      ReorderedMsProvider::StoreReorderedInMS);

  ReorderedFileWriter reordered_file_writer(handle_data,
                                            ms_polarizations_per_data_desc_id,
                                            row_provider->StartTime());

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
    reordered_file_writer.PopulateModelWithZeros(std::cref(update_progress));
  }
  progress2.reset();

  return ReorderedMsProvider::ReorderedHandle(
      ms_path, data_column_name, model_column_name, model_storage_manager,
      temporary_directory, channels, initial_model_required,
      model_update_required, pols_out, selection, bands, nAntennas,
      settings.saveReorder, ReorderedMsProvider::StoreReorderedInMS);
}  // namespace wsclean

void ReorderedMsProvider::StoreReorderedInMS(
    const ReorderedHandleData& handle) {
  const std::set<aocommon::PolarizationEnum> pols = handle.polarizations_;

  const std::map<size_t, size_t> data_desc_ids =
      GetDataDescIdMap(handle.channels_);

  std::vector<MetaHeader> meta_headers(data_desc_ids.size());
  for (const std::pair<const size_t, size_t>& data_desc_id : data_desc_ids) {
    std::ifstream meta_file(GetMetaFilename(
        handle.ms_path_, handle.temporary_directory_, data_desc_id.first));
    MetaHeader& meta_header = meta_headers[data_desc_id.second];
    meta_header.Read(meta_file);
    std::vector<char> ms_path(meta_header.filename_length + 1, char(0));
    meta_file.read(ms_path.data(), meta_header.filename_length);
  }

  ChannelRange firstRange = handle.channels_[0];
  std::ifstream first_data_file(
      GetPartPrefix(handle.ms_path_, 0, *pols.begin(), firstRange.data_desc_id,
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
      size_t data_desc_id = handle.channels_[part].data_desc_id;
      for (aocommon::PolarizationEnum p : pols) {
        std::string part_prefix =
            GetPartPrefix(handle.ms_path_, part, p, data_desc_id,
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
    const std::map<size_t, std::set<aocommon::PolarizationEnum>>
        ms_polarizations_per_data_desc_id =
            GetMSPolarizationsPerDataDescId(handle.channels_, ms);
    InitializeModelColumn(ms, handle.model_column_name_,
                          handle.model_storage_manager_);
    casacore::ScalarColumn<int> antenna1_column(
        ms, ms.columnName(casacore::MSMainEnums::ANTENNA1));
    casacore::ScalarColumn<int> antenna2_column(
        ms, ms.columnName(casacore::MSMainEnums::ANTENNA2));
    casacore::ScalarColumn<int> field_id_column(
        ms, ms.columnName(casacore::MSMainEnums::FIELD_ID));
    casacore::ScalarColumn<double> time_column(
        ms, ms.columnName(casacore::MSMainEnums::TIME));
    casacore::ScalarColumn<int> data_desc_id_column(
        ms, ms.columnName(casacore::MSMainEnums::DATA_DESC_ID));
    casacore::ArrayColumn<casacore::Complex> data_column(
        ms, handle.data_column_name_);
    casacore::ArrayColumn<casacore::Complex> model_column(
        ms, handle.model_column_name_);
    casacore::ArrayColumn<double> uvw_column(
        ms, ms.columnName(casacore::MSMainEnums::UVW));

    const casacore::IPosition shape(data_column.shape(0));
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
        std::map<size_t, size_t>::const_iterator data_desc_id_iter =
            data_desc_ids.find(data_desc_id);
        if (data_desc_id_iter != data_desc_ids.end()) {
          model_column.get(row, model_data_array, true);
          size_t file_index = 0;
          for (size_t part = 0; part != channel_parts; ++part) {
            const size_t partDataDescId = handle.channels_[part].data_desc_id;
            if (data_desc_id == partDataDescId) {
              const size_t part_start_ch = handle.channels_[part].start;
              const size_t part_end_ch = handle.channels_[part].end;
              const std::set<aocommon::PolarizationEnum>& ms_polarizations =
                  ms_polarizations_per_data_desc_id.find(data_desc_id)->second;
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
        }
      }
    }
    progress.SetProgress(end_row - start_row, end_row - start_row);

    Logger::Debug
        << "Row count for writing reordered data back to the measurement set: "
        << selected_row_count_for_debug << '\n';
  }
}

void ReorderedMsProvider::MakeIdToMSRowMapping(
    std::vector<size_t>& id_to_ms_row) {
  const MSSelection& selection = handle_.data_->selection_;
  const std::map<size_t, size_t> data_desc_ids =
      GetDataDescIdMap(handle_.data_->channels_);
  std::set<size_t> data_desc_idSet;
  for (std::map<size_t, size_t>::const_iterator i = data_desc_ids.begin();
       i != data_desc_ids.end(); ++i)
    data_desc_idSet.insert(i->first);
  size_t start_row, end_row;
  SynchronizedMS ms = MS();
  GetRowRangeAndIDMap(*ms, selection, start_row, end_row, data_desc_idSet,
                      id_to_ms_row);
}

void ReorderedMsProvider::ReorderedHandle::Serialize(
    aocommon::SerialOStream& stream) const {
  stream.Ptr(data_);
}

void ReorderedMsProvider::ReorderedHandle::Unserialize(
    aocommon::SerialIStream& stream) {
  stream.Ptr(data_);
}

}  // namespace wsclean
