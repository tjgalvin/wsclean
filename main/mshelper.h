#ifndef WSCLEAN_MSHELPER_H_
#define WSCLEAN_MSHELPER_H_

#include <memory>
#include <vector>

#include <aocommon/multibanddata.h>

#include <schaapcommon/reordering/reorderedhandle.h>

#include "../msproviders/reorderedmsprovider.h"
#include "../structures/imagingtable.h"
#include "../structures/mslistitem.h"
#include "../structures/msselection.h"

#include "settings.h"

namespace wsclean {

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

  const std::vector<ReorderedMsProvider::ReorderedHandle>&
  GetReorderedMsHandles() const {
    return reordered_ms_handles_;
  }

  const std::vector<schaapcommon::reordering::ChannelRange> GenerateChannelInfo(
      const ImagingTable& imaging_table, size_t ms_index) const;

  void ReuseReorderedFiles(const ImagingTable& imaging_table);

  void PerformReordering(const ImagingTable& imaging_table,
                         bool is_predict_mode);

  std::vector<MsListItem> InitializeMsList(
      const ImagingTableEntry& entry) const;

 private:
  const Settings& settings_;
  const schaapcommon::reordering::MSSelection& global_selection_;
  const std::vector<aocommon::MultiBandData>& ms_bands_;
  std::vector<ReorderedMsProvider::ReorderedHandle> reordered_ms_handles_;
};

}  // namespace wsclean

#endif
