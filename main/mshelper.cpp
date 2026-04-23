#include "mshelper.h"

#include <mutex>
#include <semaphore>

#include <aocommon/counting_semaphore.h>
#include <aocommon/dynamicfor.h>
#include <aocommon/logger.h>

#include "../io/readinmemory.h"
#include "../msproviders/operations.h"
#include "../msproviders/reordering/reorderms.h"
#include "../msproviders/reordering/writetoms.h"

using aocommon::Logger;

using schaapcommon::reordering::ChannelRange;
using schaapcommon::reordering::MSSelection;

namespace wsclean {
namespace {

/**
 * Sets the part index for all the entries in the image table that are supposed to have
 * the same part index. Different facets and different polarizations (either independent or joined)
 * all have the same part index, and are all set.
 */
void SetPartIndex(const ImagingTable& imaging_table, size_t output_channel, size_t sq_index,
                  size_t ms_index, size_t part_index) {
  for (size_t i = 0; i != imaging_table.IndependentGroupCount(); ++i) {
    const ImagingTable independent_group = imaging_table.GetIndependentGroup(i);
    const ImagingTable squared_group = independent_group.GetSquaredGroup(sq_index);
    for (ImagingTableEntry& entry : squared_group) {
      if (entry.outputChannelIndex == output_channel) {
        entry.part_index_per_ms[ms_index] = part_index;
      }
    }
  }
}
}  // namespace

aocommon::VectorMap<ChannelRange> MsHelper::GenerateChannelPartInfo(const ImagingTableEntry& entry,
                                                                    size_t ms_index) const {
  const aocommon::MultiBandData& band_data = ms_bands_[ms_index];
  aocommon::VectorMap<ChannelRange> result;
  for (const size_t data_desc_id : band_data.DataDescIds()) {
    const size_t band_index = band_data.GetBandIndex(data_desc_id);
    MSSelection selection(global_selection_);
    if (settings_.IsBandSelected(band_index)) {
      const ChannelRange range =
          SelectMsChannels(band_data, data_desc_id, entry.lowestFrequency, entry.highestFrequency);
      if (!range.Empty()) result.AlwaysEmplace(data_desc_id, range);
    }
  }
  return result;
}

std::vector<aocommon::VectorMap<ChannelRange>> MsHelper::GenerateChannelInfo(
    const ImagingTable& imaging_table, size_t ms_index) const {
  assert(imaging_table.IndependentGroupCount() > 0);
  std::vector<aocommon::VectorMap<ChannelRange>> channels;
  std::map<aocommon::PolarizationEnum, size_t> next_index;
  // Independent groups may be independent because they are a non-joined polarization or a
  // non-joined output channel index. Polarizations have the same channel info, but output channels
  // do not. The first polarization is used as template.
  const aocommon::PolarizationEnum first_polarization = imaging_table.Front().polarization;
  for (size_t independent_index = 0; independent_index != imaging_table.IndependentGroupCount();
       ++independent_index) {
    const ImagingTable independent_group = imaging_table.GetIndependentGroup(independent_index);
    for (size_t sq_index = 0; sq_index != independent_group.SquaredGroupCount(); ++sq_index) {
      const ImagingTable squared_group = independent_group.GetSquaredGroup(sq_index);
      // One facet_group contains all the entries for one output channel and
      // one polarization. Like for polarizations, the band information is determined from the first
      // facet in the group of the first polarization. After this, all facet entries and
      // polarizations inside the facet_groups for the same output channel are updated using
      // SetPartIndex().
      const ImagingTable::Groups facet_groups =
          squared_group.FacetGroups([first_polarization](const ImagingTableEntry& entry) {
            return entry.polarization == first_polarization;
          });
      if (!facet_groups.empty()) {
        const ImagingTableEntry& entry = *facet_groups.front().front();
        aocommon::VectorMap<ChannelRange> range = GenerateChannelPartInfo(entry, ms_index);
        if (!range.Empty()) {
          const size_t part_index = channels.size();
          channels.emplace_back(std::move(range));
          SetPartIndex(imaging_table, entry.outputChannelIndex, sq_index, ms_index, part_index);
        }
      }
    }
  }

  return channels;
}

void MsHelper::ReuseReorderedFiles(const ImagingTable& imaging_table) {
  assert(reordered_ms_handles_.empty());

  reordered_ms_handles_.resize(settings_.filenames.size());
  bool initial_model_required = settings_.subtractModel || settings_.continuedRun;

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
        settings_.filenames[ms_index], settings_.dataColumnName, settings_.modelColumnName,
        settings_.modelStorageManager, settings_.temporaryDirectory, channels,
        initial_model_required, settings_.modelUpdateRequired, polarization_types,
        global_selection_, settings_.sortMsInTime, bands_per_part, has_frequency_bda, n_antennas,
        settings_.saveReorder, WriteReorderedToMS);

    std::vector<aocommon::OptionalNumber<size_t>> data_desc_ids;
    std::tie(data->metadata_indices_, data_desc_ids) =
        schaapcommon::reordering::MakeMetaFilesMap(data->channels_);
    ReorderedHandle part_ms = ReorderedHandle(std::move(data));

    reordered_ms_handles_[ms_index] = std::move(part_ms);
  }
}

void MsHelper::PerformReordering(const ImagingTable& imaging_table, bool is_predict_mode) {
  assert(reordered_ms_handles_.empty());
  assert(in_memory_handles_.empty());

  std::mutex mutex;

  reordered_ms_handles_.resize(settings_.filenames.size());
  bool use_model = settings_.deconvolutionMGain != 1.0 || is_predict_mode ||
                   settings_.subtractModel || settings_.continuedRun;
  bool initial_model_required = settings_.subtractModel || settings_.continuedRun;

  if (settings_.parallelReordering != 1) Logger::Info << "Reordering...\n";

  std::counting_semaphore semaphore(settings_.parallelReordering);
  aocommon::DynamicFor<size_t> loop;
  loop.Run(0, settings_.filenames.size(), [&](size_t ms_index) {
    aocommon::ScopedCountingSemaphoreLock semaphore_lock(semaphore);
    std::vector<aocommon::VectorMap<ChannelRange>> channels =
        GenerateChannelInfo(imaging_table, ms_index);
    ReorderedHandle part_ms =
        ReorderMS(settings_.filenames[ms_index], channels, global_selection_,
                  settings_.dataColumnName, settings_.modelColumnName,
                  settings_.modelStorageManager, use_model, initial_model_required, settings_);
    std::lock_guard<std::mutex> lock(mutex);
    reordered_ms_handles_[ms_index] = std::move(part_ms);
    if (settings_.parallelReordering != 1)
      Logger::Info << "Finished reordering " << settings_.filenames[ms_index] << " [" << ms_index
                   << "]\n";
  });
}

void MsHelper::PrepareInMemoryData(const ImagingTable& imaging_table, bool is_predict_mode) {
  assert(reordered_ms_handles_.empty());
  assert(in_memory_handles_.empty());
  // This function handles two cases:
  // - If the user asked for in-memory data via the cmdline, the data must still
  // be read
  //   from disk and then reordered.
  // - If the API was used and data was supplied, it only needs to be reordered.
  // If the following variable is true, we are in the second case:
  const bool is_in_memory_data_supplied = !in_memory_data_before_reordering_.empty();
  assert((!is_in_memory_data_supplied && !settings_.filenames.empty()) ||
         (is_in_memory_data_supplied && settings_.filenames.empty()));

  const size_t ms_count = is_in_memory_data_supplied ? in_memory_data_before_reordering_.size()
                                                     : settings_.filenames.size();
  in_memory_handles_.resize(ms_count);

  const size_t parallel_reordering = std::clamp<size_t>(settings_.parallelReordering, 0, ms_count);
  if (parallel_reordering != 1) {
    if (is_in_memory_data_supplied)
      Logger::Info << "In-memory reordering of the data...\n";
    else
      Logger::Info << "Reading into memory...\n";
  }

  std::mutex mutex;
  std::counting_semaphore semaphore(parallel_reordering);
  aocommon::DynamicFor<size_t> loop;
  loop.Run(0, ms_count, [&](size_t ms_index) {
    aocommon::ScopedCountingSemaphoreLock semaphore_lock(semaphore);

    InMemoryMs in_memory_ms;

    if (is_in_memory_data_supplied) {
      in_memory_ms = std::move(*in_memory_data_before_reordering_[ms_index]);
      ApplyRowSelectionToMs(in_memory_ms, global_selection_);
    } else {
      in_memory_ms = ReadMsInMemory(settings_.filenames[ms_index], settings_.dataColumnName,
                                    global_selection_, settings_.sortMsInTime);
    }

    const std::vector<aocommon::VectorMap<ChannelRange>> channels =
        GenerateChannelInfo(imaging_table, ms_index);
    // The model data is always overwritten in the first iteration, so it does
    // not need to be initialized.
    constexpr bool initial_model_required = false;
    InMemoryHandle handle =
        ReorderInMemory(std::move(in_memory_ms), channels, settings_.dataColumnName,
                        settings_.modelColumnName, initial_model_required, settings_);

    std::lock_guard<std::mutex> lock(mutex);
    in_memory_handles_[ms_index] = std::move(handle);
    if (settings_.parallelReordering != 1) {
      if (is_in_memory_data_supplied) {
        Logger::Info << "Finished reordering ms " << ms_index << "\n";
      } else {
        Logger::Info << "Finished reading " << settings_.filenames[ms_index] << " [" << ms_index
                     << "]\n";
      }
    }
  });
}

std::vector<MsListItem> MsHelper::InitializeMsList(const ImagingTableEntry& entry) {
  const aocommon::PolarizationEnum polarization =
      settings_.GetProviderPolarization(entry.polarization);

  std::vector<MsListItem> ms_list;

  const size_t ms_count =
      settings_.inMemoryData ? in_memory_handles_.size() : settings_.filenames.size();
  for (size_t ms_index = 0; ms_index != ms_count; ++ms_index) {
    MsListItem item;
    aocommon::VectorMap<ChannelRange> channels = GenerateChannelPartInfo(entry, ms_index);
    if (!channels.Empty()) {
      if (settings_.inMemoryData) {
        item.ms_description = MSDataDescription::ForInMemory(in_memory_handles_[ms_index],
                                                             entry.part_index_per_ms[ms_index],
                                                             polarization, settings_.UseMpi());
      } else if (settings_.doReorder) {
        item.ms_description = MSDataDescription::ForReordered(reordered_ms_handles_[ms_index],
                                                              entry.part_index_per_ms[ms_index],
                                                              polarization, settings_.UseMpi());
      } else {
        item.ms_description = MSDataDescription::ForContiguous(
            settings_.filenames[ms_index], settings_.dataColumnName, settings_.modelColumnName,
            settings_.sortMsInTime, settings_.modelStorageManager, global_selection_, channels,
            polarization, settings_.UseMpi());
      }
      item.ms_index = ms_index;
      ms_list.emplace_back(std::move(item));
    }
  }

  return ms_list;
}

}  // namespace wsclean
