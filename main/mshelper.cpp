#include "mshelper.h"

#include <mutex>

#include <aocommon/counting_semaphore.h>
#include <aocommon/dynamicfor.h>
#include <aocommon/logger.h>

using aocommon::Logger;

using schaapcommon::reordering::ChannelRange;
using schaapcommon::reordering::MSSelection;

namespace wsclean {

const std::vector<ChannelRange> MsHelper::GenerateChannelInfo(
    const ImagingTable& imaging_table, size_t ms_index) const {
  const aocommon::MultiBandData& band_data = ms_bands_[ms_index];
  std::vector<ChannelRange> channels;
  // The partIndex needs to increase per data desc ids and channel ranges
  std::map<aocommon::PolarizationEnum, size_t> next_index;
  for (size_t sq_index = 0; sq_index != imaging_table.SquaredGroupCount();
       ++sq_index) {
    const ImagingTable::Groups facet_groups =
        imaging_table.FacetGroups([sq_index](const ImagingTableEntry& e) {
          return e.squaredDeconvolutionIndex == sq_index;
        });
    for (const ImagingTable::Group& facet_group : facet_groups) {
      // The band information is determined from the first facet in the group.
      // After this, all facet entries inside the group are updated.
      const ImagingTableEntry& entry = *facet_group.front();
      for (size_t d = 0; d != band_data.DataDescCount(); ++d) {
        MSSelection selection{global_selection_};
        const size_t band_index = band_data.GetBandIndex(d);

        if (settings_.IsBandSelected(band_index) &&
            SelectMsChannels(selection, band_data, d, entry)) {
          if (entry.polarization == *settings_.polarizations.begin()) {
            ChannelRange r;
            r.data_desc_id = d;
            r.start = selection.ChannelRangeStart();
            r.end = selection.ChannelRangeEnd();
            channels.push_back(r);
          }
          for (const std::shared_ptr<ImagingTableEntry>& facet_entry :
               facet_group) {
            facet_entry->msData[ms_index].bands[d].partIndex =
                next_index[entry.polarization];
          }
          ++next_index[entry.polarization];
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
    const std::vector<ChannelRange> channels =
        GenerateChannelInfo(imaging_table, ms_index);
    casacore::MeasurementSet ms_data_obj(settings_.filenames[ms_index]);
    const size_t n_antennas = ms_data_obj.antenna().nrow();
    const aocommon::MultiBandData bands(ms_data_obj);
    ReorderedMsProvider::ReorderedHandle part_ms =
        ReorderedMsProvider::ReorderedHandle(
            settings_.filenames[ms_index], settings_.dataColumnName,
            settings_.modelColumnName, settings_.modelStorageManager,
            settings_.temporaryDirectory, channels, initial_model_required,
            settings_.modelUpdateRequired, polarization_types,
            global_selection_, bands, n_antennas, settings_.saveReorder,
            ReorderedMsProvider::StoreReorderedInMS);

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
    std::vector<ChannelRange> channels =
        GenerateChannelInfo(imaging_table, ms_index);
    ReorderedMsProvider::ReorderedHandle part_ms =
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
    const aocommon::MultiBandData& band_data = ms_bands_[ms_index];

    for (size_t data_description_id = 0;
         data_description_id != band_data.DataDescCount();
         ++data_description_id) {
      MSSelection selection{global_selection_};
      const size_t band_index = band_data.GetBandIndex(data_description_id);

      if (settings_.IsBandSelected(band_index) &&
          SelectMsChannels(selection, band_data, data_description_id, entry)) {
        MsListItem item;
        if (settings_.doReorder)
          item.ms_description = MSDataDescription::ForReordered(
              reordered_ms_handles_[ms_index], selection,
              entry.msData[ms_index].bands[data_description_id].partIndex,
              polarization, data_description_id, settings_.UseMpi());
        else
          item.ms_description = MSDataDescription::ForContiguous(
              settings_.filenames[ms_index], settings_.dataColumnName,
              settings_.modelColumnName, settings_.modelStorageManager,
              selection, polarization, data_description_id, settings_.UseMpi());
        item.ms_index = ms_index;
        ms_list.emplace_back(std::move(item));
      }
    }
  }

  return ms_list;
}

}  // namespace wsclean
