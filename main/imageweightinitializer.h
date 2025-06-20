#ifndef WSCLEAN_IMAGE_WEIGHT_INITIALIZER_H_
#define WSCLEAN_IMAGE_WEIGHT_INITIALIZER_H_

#include <memory>
#include <vector>

#include <aocommon/multibanddata.h>
#include <aocommon/polarization.h>

#include "../io/imageweightcache.h"
#include "../msproviders/reorderedmsprovider.h"
#include "../structures/imageweights.h"
#include "../structures/imagingtable.h"
#include "../structures/mslistitem.h"
#include "../structures/msselection.h"

#include "settings.h"

namespace wsclean {

/**
 * Helper class with various initialization routines for WSClean.
 */
class ImageWeightInitializer {
 public:
  explicit ImageWeightInitializer(
      const Settings& settings,
      const schaapcommon::reordering::MSSelection& global_selection,
      const std::vector<aocommon::MultiBandData>& ms_bands,
      const std::vector<ReorderedMsProvider::ReorderedHandle>&
          reordered_ms_handles)
      : settings_(settings),
        global_selection_(global_selection),
        ms_bands_(ms_bands),
        reordered_ms_handles_(reordered_ms_handles) {}

  const Settings& GetSettings() const { return settings_; }

  std::shared_ptr<ImageWeights> Initialize(
      const ImagingTableEntry& entry, const std::vector<MsListItem>& ms_list,
      ImageWeightCache& cache) const;

  void InitializeMf(const ImagingTable& imaging_table, ImageWeightCache& cache);

 private:
  void GridMfReorderedBand(size_t data_desc_id,
                           const aocommon::MultiBandData& bands,
                           ImageWeights& weights, size_t ms_index,
                           const ImagingTableEntry& entry) const;
  void GridMfContiguousBand(size_t filename_index, size_t data_desc_id,
                            ImageWeights& weights) const;

  const Settings& settings_;
  const schaapcommon::reordering::MSSelection& global_selection_;
  const std::vector<aocommon::MultiBandData>& ms_bands_;
  const std::vector<ReorderedMsProvider::ReorderedHandle>&
      reordered_ms_handles_;
};

}  // namespace wsclean

#endif
