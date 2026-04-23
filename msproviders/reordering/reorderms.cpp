#include "reorderms.h"

#include <aocommon/logger.h>
#include <aocommon/polarization.h>
#include <aocommon/vectorset.h>

#include <schaapcommon/reordering/channelrange.h>
#include <schaapcommon/reordering/filewriter.h>
#include <schaapcommon/reordering/handledata.h>
#include <schaapcommon/reordering/msselection.h>

#include "writetoms.h"

#include "../../main/progressbar.h"
#include "../../main/settings.h"
#include "../rowproviders/averagingmsrowprovider.h"
#include "../rowproviders/noisemsrowprovider.h"
#include "../operations.h"

using aocommon::Logger;

namespace wsclean {

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
    const std::vector<aocommon::VectorMap<schaapcommon::reordering::ChannelRange>>& channels,
    const schaapcommon::reordering::MSSelection& selection, const std::string& data_column_name,
    const std::string& model_column_name,
    schaapcommon::reordering::StorageManagerType model_storage_manager, bool include_model,
    bool initial_model_required, const Settings& settings) {
  const bool model_update_required = settings.modelUpdateRequired;
  std::set<aocommon::PolarizationEnum> pols_out;
  for (aocommon::PolarizationEnum p : settings.polarizations)
    pols_out.insert(settings.GetProviderPolarization(p));
  const std::string& temporary_directory = settings.temporaryDirectory;

  const size_t channel_parts = channels.size();

  aocommon::VectorSet<size_t> selected_data_desc_ids;
  bool has_frequency_bda;
  {  // close ms asap
    casacore::MeasurementSet ms(ms_path);
    selected_data_desc_ids = aocommon::MultiBandData(ms).GetDataDescIdSet();
    has_frequency_bda = wsclean::HasFrequencyBda(ms);
  }

  std::unique_ptr<MsRowProviderBase> row_provider;
  if (settings.baselineDependentAveragingInWavelengths == 0.0) {
    if (settings.simulateNoise) {
      std::unique_ptr<NoiseMSRowProvider> noise_row_provider(
          new NoiseMSRowProvider(ms_path, selection, selected_data_desc_ids, data_column_name,
                                 model_column_name, initial_model_required));
      if (settings.simulatedBaselineNoiseFilename.empty())
        noise_row_provider->SetNoiseLevel(settings.simulatedNoiseStdDev);
      else
        noise_row_provider->SetNoiseBaselineFile(settings.simulatedBaselineNoiseFilename);
      row_provider = std::move(noise_row_provider);
    } else
      row_provider =
          MakeMsRowProvider(ms_path, selection, selected_data_desc_ids, data_column_name,
                            model_column_name, initial_model_required, settings.sortMsInTime);
  } else {
    if (initial_model_required) {
      throw std::runtime_error(
          "Baseline-dependent averaging is enabled together with a mode that "
          "requires the model data (e.g. -continue or -subtract-model). This "
          "is not possible.");
    }
    row_provider = std::make_unique<AveragingMSRowProvider>(
        settings.baselineDependentAveragingInWavelengths, ms_path, selection,
        selected_data_desc_ids, settings.fieldIds[0], data_column_name, model_column_name,
        initial_model_required);
  }

  const std::map<size_t, std::set<aocommon::PolarizationEnum>> ms_polarizations_per_data_desc_id =
      GetMsPolarizationsPerDataDescId(row_provider->Ms());
  const size_t nAntennas = row_provider->Ms().antenna().nrow();
  const aocommon::MultiBandData original_bands(row_provider->Ms());

  // This handle is just for the writer
  const std::vector<aocommon::MultiBandData> bands_per_part =
      MakeSelectedBands(original_bands, channels);
  auto handle_data = std::make_unique<schaapcommon::reordering::HandleData>(
      ms_path, data_column_name, model_column_name, model_storage_manager, temporary_directory,
      channels, initial_model_required, model_update_required, pols_out, selection,
      settings.sortMsInTime, bands_per_part, has_frequency_bda, nAntennas, settings.saveReorder,
      WriteReorderedToMS);

  std::vector<aocommon::OptionalNumber<size_t>> data_desc_ids;
  std::tie(handle_data->metadata_indices_, data_desc_ids) =
      schaapcommon::reordering::MakeMetaFilesMap(handle_data->channels_);

  Logger::Debug << "Using " << data_desc_ids.size() << " temporary metadata files.\n";

  schaapcommon::reordering::FileWriter reordered_file_writer(
      *handle_data, ms_polarizations_per_data_desc_id, data_desc_ids, row_provider->StartTime());

  if (settings.parallelReordering == 1) {
    Logger::Info << "Reordering " << ms_path << " into " << channel_parts << " x "
                 << pols_out.size() << " parts.\n";
  }

  casacore::Array<std::complex<float>> data_array;
  casacore::Array<std::complex<float>> model_array;
  casacore::Array<float> weight_spectrum_array;
  casacore::Array<bool> flag_array;

  std::unique_ptr<ProgressBar> reordering_progress;
  if (settings.parallelReordering == 1) reordering_progress.reset(new ProgressBar("Reordering"));

  size_t selected_rows_total = 0;
  while (!row_provider->AtEnd()) {
    if (reordering_progress)
      reordering_progress->SetProgress(row_provider->CurrentProgress(),
                                       row_provider->TotalProgress());

    double time, u, v, w;
    uint32_t data_desc_id, antenna1, antenna2, field_id;
    row_provider->ReadData(data_array, flag_array, weight_spectrum_array, u, v, w, data_desc_id,
                           antenna1, antenna2, field_id, time);

    if (initial_model_required) row_provider->ReadModel(model_array);

    reordered_file_writer.WriteMetaRow(u, v, w, time, data_desc_id, antenna1, antenna2, field_id);
    reordered_file_writer.WriteDataRow(data_array.data(), model_array.data(),
                                       weight_spectrum_array.data(), flag_array.data(),
                                       data_desc_id);

    row_provider->NextRow();
    ++selected_rows_total;
  }
  reordering_progress.reset();
  Logger::Debug << "Total selected rows: " << selected_rows_total << '\n';
  row_provider->OutputStatistics();

  reordered_file_writer.UpdateMetaHeaders();
  reordered_file_writer.UpdatePartHeaders(include_model);

  std::unique_ptr<ProgressBar> model_column_progress;
  if (include_model && !initial_model_required) {
    if (settings.parallelReordering) {
      model_column_progress = std::make_unique<ProgressBar>("Initializing model visibilities");
    }
    auto update_progress = [progress_bar = std::move(model_column_progress)](size_t progress,
                                                                             size_t total) {
      if (progress_bar) {
        progress_bar->SetProgress(progress, total);
      }
    };
    constexpr bool zero_initialize_model = true;
    reordered_file_writer.PopulateModel(zero_initialize_model, std::cref(update_progress));
  }
  model_column_progress.reset();

  return ReorderedHandle(std::move(handle_data));
}

}  // namespace wsclean
