#ifndef WSCLEAN_MSHELPER_H_
#define WSCLEAN_MSHELPER_H_

#include <memory>
#include <vector>

#include <aocommon/multibanddata.h>
#include <aocommon/vectormap.h>

#include <schaapcommon/reordering/handledata.h>

#include "../msproviders/reorderedmsprovider.h"

#include "../structures/imagingtable.h"
#include "../structures/inmemoryreordering.h"
#include "../structures/mslistitem.h"
#include "../structures/msselection.h"

#include "settings.h"

namespace wsclean {

struct InMemoryMs;

/**
 * Class with helper routines for managing measurement sets.
 */
class MsHelper {
 public:
  /**
   * @param ms_bands List such that element ms_bands[i] holds the bands for
   * settings.filenames[i].
   */
  explicit MsHelper(
      const Settings& settings,
      const schaapcommon::reordering::MSSelection& global_selection,
      const std::vector<aocommon::MultiBandData>& ms_bands)
      : settings_{settings},
        global_selection_{global_selection},
        ms_bands_{ms_bands},
        reordered_ms_handles_{} {}

  const std::vector<ReorderedHandle>& GetReorderedMsHandles() const {
    return reordered_ms_handles_;
  }

  std::vector<aocommon::VectorMap<schaapcommon::reordering::ChannelRange>>
  GenerateChannelInfo(const ImagingTable& imaging_table, size_t ms_index) const;

  aocommon::VectorMap<schaapcommon::reordering::ChannelRange>
  GenerateChannelPartInfo(const ImagingTableEntry& entry,
                          size_t ms_index) const;

  void SetInMemoryData(
      std::vector<std::unique_ptr<InMemoryMs>>&& in_memory_data) {
    in_memory_data_before_reordering_ = std::move(in_memory_data);
  }

  void ReuseReorderedFiles(const ImagingTable& imaging_table);

  void PrepareInMemoryData(const ImagingTable& imaging_table,
                           bool is_predict_mode);

  void PerformReordering(const ImagingTable& imaging_table,
                         bool is_predict_mode);

  std::vector<MsListItem> InitializeMsList(const ImagingTableEntry& entry);

 private:
  const Settings& settings_;
  const schaapcommon::reordering::MSSelection& global_selection_;
  const std::vector<aocommon::MultiBandData>& ms_bands_;
  std::vector<ReorderedHandle> reordered_ms_handles_;
  std::vector<InMemoryHandle> in_memory_handles_;
  std::vector<std::unique_ptr<InMemoryMs>> in_memory_data_before_reordering_;
};

}  // namespace wsclean

#endif
