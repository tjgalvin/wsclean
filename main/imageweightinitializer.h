#ifndef WSCLEAN_IMAGE_WEIGHT_INITIALIZER_H_
#define WSCLEAN_IMAGE_WEIGHT_INITIALIZER_H_

#include <memory>
#include <vector>

namespace wsclean {

class ImageWeights;
class ImageWeightCache;
class ImagingTable;
class ImagingTableEntry;
class MsHelper;
class MsListItem;
class Settings;

std::shared_ptr<ImageWeights> InitializeWeights(
    const Settings& settings, const ImagingTableEntry& entry,
    const std::vector<MsListItem>& ms_list, ImageWeightCache& cache);

void InitializeMfWeights(const Settings& settings,
                         const ImagingTable& imaging_table,
                         ImageWeightCache& cache, MsHelper& ms_helper);

}  // namespace wsclean

#endif
