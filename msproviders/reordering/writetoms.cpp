#include "writetoms.h"

#include <fstream>
#include <set>

#include <aocommon/logger.h>
#include <aocommon/polarization.h>
#include <aocommon/vectorset.h>

#include <schaapcommon/reordering/handledata.h>
#include <schaapcommon/reordering/reordering.h>

#include "bdarowwriter.h"
#include "regularrowwriter.h"
#include "../rowproviders/msrowprovider.h"
#include "../operations.h"
#include "../../main/progressbar.h"

namespace wsclean {

void WriteReorderedToMS(const schaapcommon::reordering::HandleData& handle) {
  const std::set<aocommon::PolarizationEnum> pols = handle.polarizations_;

  const std::string prefix = schaapcommon::reordering::GetPartPrefix(
      handle.ms_path_, 0, *pols.begin(), handle.temporary_directory_);
  std::ifstream first_data_file(prefix + ".tmp", std::ios::in);
  if (!first_data_file.good()) throw std::runtime_error("Error opening temporary data file");
  schaapcommon::reordering::PartHeader firstpart_header_;
  firstpart_header_.Read(first_data_file);
  if (!first_data_file.good()) throw std::runtime_error("Error reading from temporary data file");

  if (firstpart_header_.has_model) {
    const size_t channel_parts = handle.channels_.size();

    // Open the temporary files
    std::vector<std::unique_ptr<std::ifstream>> model_files(channel_parts * pols.size());
    size_t file_index = 0;
    for (size_t part = 0; part != channel_parts; ++part) {
      for (aocommon::PolarizationEnum p : pols) {
        const std::string part_prefix = schaapcommon::reordering::GetPartPrefix(
            handle.ms_path_, part, p, handle.temporary_directory_);
        model_files[file_index] = std::make_unique<std::ifstream>(part_prefix + "-m.tmp");
        if (!*model_files[file_index])
          throw std::runtime_error("Error opening temporary model data file '" + part_prefix +
                                   "-m.tmp' for reading");
        ++file_index;
      }
    }

    aocommon::VectorSet<size_t> selected_data_desc_ids;
    {  // close MS as soon as possible
      casacore::MeasurementSet ms(handle.ms_path_, casacore::Table::Update);
      selected_data_desc_ids = aocommon::MultiBandData(ms).GetDataDescIdSet();
    }
    std::unique_ptr<MsRowProviderBase> row_provider = MakeMsRowProvider(
        handle.ms_path_, handle.selection_, selected_data_desc_ids, handle.data_column_name_,
        handle.model_column_name_, false, handle.sort_in_time_);
    row_provider->Ms().reopenRW();
    const aocommon::VectorMap<std::set<aocommon::PolarizationEnum>>
        ms_polarizations_per_data_desc_id = GetMsPolarizationsPerDataDescId(row_provider->Ms());
    InitializeModelColumn(row_provider->Ms(), handle.data_column_name_, handle.model_column_name_,
                          handle.model_storage_manager_);
    casacore::ArrayColumn<casacore::Complex> model_column(row_provider->Ms(),
                                                          handle.model_column_name_);
    casacore::ArrayColumn<casacore::Complex> data_column(row_provider->Ms(),
                                                         handle.data_column_name_);
    const size_t maxchannels_ = GetMaxChannels(handle.channels_);

    const bool has_time_bda = HasTimeBda(row_provider->Ms());
    std::unique_ptr<RowWriter> row_writer;
    if (has_time_bda) {
      row_writer = std::make_unique<BdaRowWriter>(row_provider->Ms(), data_column, model_column);
    } else {
      row_writer = std::make_unique<RegularRowWriter>(data_column, model_column);
    }

    const size_t polarizations_per_file = aocommon::Polarization::GetVisibilityCount(*pols.begin());
    std::vector<std::complex<float>> model_data_buffer(maxchannels_ * polarizations_per_file);
    casacore::Array<std::complex<float>> model_data_array;
    casacore::Array<double> uvw(casacore::IPosition{3});

    ProgressBar progress(std::string("Writing changed model back to ") + handle.ms_path_);

    size_t selected_row_count_for_debug = 0;
    const size_t total_progress = row_provider->TotalProgress();
    while (!row_provider->AtEnd()) {
      const MSRowProvider::WriteInfo write_info = row_provider->GetWriteInfo();
      progress.SetProgress(write_info.current_progress, total_progress);

      const casacore::IPosition shape = data_column.shape(write_info.row_index);
      model_data_array.resize(shape);

      size_t file_index = 0;
      for (size_t part = 0; part != channel_parts; ++part) {
        const schaapcommon::reordering::ChannelRange& range =
            handle.channels_[part][write_info.data_desc_id];
        if (!range.Empty()) {
          const size_t part_start_ch = range.start;
          const size_t part_end_ch = range.end;
          const std::set<aocommon::PolarizationEnum>& ms_polarizations =
              ms_polarizations_per_data_desc_id[write_info.data_desc_id];
          for (aocommon::PolarizationEnum p : pols) {
            model_files[file_index]->read(reinterpret_cast<char*>(model_data_buffer.data()),
                                          (part_end_ch - part_start_ch) * polarizations_per_file *
                                              sizeof(std::complex<float>));
            if (!model_files[file_index]->good())
              throw std::runtime_error("Error reading from temporary model data file");
            schaapcommon::reordering::StoreData<false>(model_data_array.data(), part_start_ch,
                                                       part_end_ch, ms_polarizations,
                                                       model_data_buffer.data(), p);

            ++file_index;
          }
        } else {
          file_index += pols.size();
        }
      }
      row_writer->Write(write_info.row_index, model_data_array);

      ++selected_row_count_for_debug;
      row_provider->NextRow();
    }
    progress.SetProgress(total_progress, total_progress);

    aocommon::Logger::Debug << "Row count for writing reordered data back to the measurement set: "
                            << selected_row_count_for_debug << '\n';
  }
}

}  // namespace wsclean
