#include "mshelper.h"

#include <mutex>

#include <aocommon/counting_semaphore.h>
#include <aocommon/dynamicfor.h>
#include <aocommon/logger.h>

using aocommon::Logger;

using schaapcommon::reordering::ChannelRange;
using schaapcommon::reordering::MSSelection;

namespace wsclean {

aocommon::VectorMap<ChannelRange> MsHelper::GenerateChannelPartInfo(
    const ImagingTableEntry& entry, size_t ms_index) const {
  const aocommon::MultiBandData& band_data = ms_bands_[ms_index];
  aocommon::VectorMap<ChannelRange> result;
  for (const size_t data_desc_id : band_data.DataDescIds()) {
    const size_t band_index = band_data.GetBandIndex(data_desc_id);
    MSSelection selection(global_selection_);
    if (settings_.IsBandSelected(band_index)) {
      const ChannelRange range =
          SelectMsChannels(band_data, data_desc_id, entry.lowestFrequency,
                           entry.highestFrequency);
      if (!range.Empty()) result.AlwaysEmplace(data_desc_id, range);
    }
  }
  return result;
}

std::vector<aocommon::VectorMap<ChannelRange>> MsHelper::GenerateChannelInfo(
    const ImagingTable& imaging_table, size_t ms_index) const {
  const std::string ms_name(settings_.filenames[ms_index]);
  std::vector<aocommon::VectorMap<ChannelRange>> channels;
  std::map<aocommon::PolarizationEnum, size_t> next_index;
  for (size_t sq_index = 0; sq_index != imaging_table.SquaredGroupCount();
       ++sq_index) {
    // One facet_group contains all the entries for one output channel and
    // one polarization.
    // The band information is determined from the first facet in the group of
    // the first polarization. After this, all facet entries inside the
    // facet_groups are updated.
    const ImagingTable::Groups facet_groups =
        imaging_table.FacetGroups([sq_index](const ImagingTableEntry& e) {
          return e.squaredDeconvolutionIndex == sq_index;
        });
    const ImagingTableEntry& entry = *facet_groups.front().front();
    aocommon::VectorMap<ChannelRange> range =
        GenerateChannelPartInfo(entry, ms_index);
    if (!range.Empty()) {
      const size_t part_index = channels.size();
      channels.emplace_back(std::move(range));
      for (const ImagingTable::Group& facet_group : facet_groups) {
        for (const std::shared_ptr<ImagingTableEntry>& entry : facet_group) {
          entry->part_index_per_ms[ms_index] = part_index;
        }
      }
    }
  }

  return channels;
}

void MsHelper::ReuseReorderedFiles(const ImagingTable& imaging_table) {
  assert(reordered_ms_handles_.empty());

  reordered_ms_handles_.resize(settings_.filenames.size());
  bool initial_model_required =
      settings_.subtractModel || settings_.continuedRun;

  Logger::Info << "Reading pre-generated reordered data...\n";

  std::set<aocommon::PolarizationEnum> polarization_types;
  for (aocommon::PolarizationEnum p : settings_.polarizations)
    polarization_types.insert(settings_.GetProviderPolarization(p));

  for (size_t ms_index = 0; ms_index < settings_.filenames.size(); ++ms_index) {
    const std::vector<aocommon::VectorMap<ChannelRange>> channels =
        GenerateChannelInfo(imaging_table, ms_index);
    casacore::MeasurementSet ms_data_obj(settings_.filenames[ms_index]);
    const bool has_frequency_bda = wsclean::HasFrequencyBda(ms_data_obj);
    const size_t n_antennas = ms_data_obj.antenna().nrow();
    const aocommon::MultiBandData original_bands(ms_data_obj);
    const std::vector<aocommon::MultiBandData> bands_per_part =
        MakeSelectedBands(original_bands, channels);

    auto data = std::make_unique<schaapcommon::reordering::HandleData>(
        settings_.filenames[ms_index], settings_.dataColumnName,
        settings_.modelColumnName, settings_.modelStorageManager,
        settings_.temporaryDirectory, channels, initial_model_required,
        settings_.modelUpdateRequired, polarization_types, global_selection_,
        bands_per_part, has_frequency_bda, n_antennas, settings_.saveReorder,
        ReorderedMsProvider::StoreReorderedInMS);

    std::vector<aocommon::OptionalNumber<size_t>> data_desc_ids;
    std::tie(data->metadata_indices_, data_desc_ids) =
        schaapcommon::reordering::MakeMetaFilesMap(data->channels_);
    ReorderedHandle part_ms = ReorderedHandle(std::move(data));

    reordered_ms_handles_[ms_index] = std::move(part_ms);
  }
}

void MsHelper::PerformReordering(const ImagingTable& imaging_table,
                                 bool is_predict_mode) {
  std::mutex mutex;

  assert(reordered_ms_handles_.empty());

  reordered_ms_handles_.resize(settings_.filenames.size());
  bool use_model = settings_.deconvolutionMGain != 1.0 || is_predict_mode ||
                   settings_.subtractModel || settings_.continuedRun;
  bool initial_model_required =
      settings_.subtractModel || settings_.continuedRun;

  if (settings_.parallelReordering != 1) Logger::Info << "Reordering...\n";

  aocommon::CountingSemaphore semaphore(settings_.parallelReordering);
  aocommon::DynamicFor<size_t> loop;
  loop.Run(0, settings_.filenames.size(), [&](size_t ms_index) {
    aocommon::ScopedCountingSemaphoreLock semaphore_lock(semaphore);
    std::vector<aocommon::VectorMap<ChannelRange>> channels =
        GenerateChannelInfo(imaging_table, ms_index);
    ReorderedHandle part_ms =
        ReorderMS(settings_.filenames[ms_index], channels, global_selection_,
                  settings_.dataColumnName, settings_.modelColumnName,
                  settings_.modelStorageManager, use_model,
                  initial_model_required, settings_);
    std::lock_guard<std::mutex> lock(mutex);
    reordered_ms_handles_[ms_index] = std::move(part_ms);
    if (settings_.parallelReordering != 1)
      Logger::Info << "Finished reordering " << settings_.filenames[ms_index]
                   << " [" << ms_index << "]\n";
  });
}

std::vector<MsListItem> MsHelper::InitializeMsList(
    const ImagingTableEntry& entry) const {
  const aocommon::PolarizationEnum polarization =
      settings_.GetProviderPolarization(entry.polarization);

  std::vector<MsListItem> ms_list;

  for (size_t ms_index = 0; ms_index != settings_.filenames.size();
       ++ms_index) {
    MsListItem item;
    aocommon::VectorMap<ChannelRange> channels =
        GenerateChannelPartInfo(entry, ms_index);
    if (!channels.Empty()) {
      if (settings_.doReorder) {
        item.ms_description = MSDataDescription::ForReordered(
            reordered_ms_handles_[ms_index], entry.part_index_per_ms[ms_index],
            polarization, settings_.UseMpi());
      } else {
        item.ms_description = MSDataDescription::ForContiguous(
            settings_.filenames[ms_index], settings_.dataColumnName,
            settings_.modelColumnName, settings_.modelStorageManager,
            global_selection_, channels, polarization, settings_.UseMpi());
      }
      item.ms_index = ms_index;
      ms_list.emplace_back(std::move(item));
    }
  }

  return ms_list;
}

}  // namespace wsclean
