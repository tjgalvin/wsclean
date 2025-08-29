#include "imageweightinitializer.h"

#include <aocommon/logger.h>

#include "../io/imagefilename.h"
#include "../io/imageweightcache.h"

#include "../msproviders/msprovider.h"

#include "../structures/imageweights.h"
#include "../structures/imagingtable.h"
#include "../structures/mslistitem.h"

#include "mshelper.h"
#include "settings.h"

using aocommon::Logger;
using schaapcommon::reordering::MSSelection;

namespace wsclean {

std::shared_ptr<ImageWeights> InitializeWeights(
    const Settings& settings, const ImagingTableEntry& entry,
    const std::vector<MsListItem>& ms_list, ImageWeightCache& cache) {
  if (settings.mfWeighting) {
    return cache.GetMFWeights();
  } else {
    std::shared_ptr<ImageWeights> weights =
        cache.Get(ms_list, entry.outputChannelIndex, entry.outputIntervalIndex);
    if (settings.isWeightImageSaved) {
      const std::string prefix = ImageFilename::GetPSFPrefix(
          settings, entry.outputChannelIndex, entry.outputIntervalIndex);
      weights->Save(prefix + "-weights.fits");
    }
    return weights;
  }
}

void InitializeMfWeights(const Settings& settings,
                         const ImagingTable& imaging_table,
                         ImageWeightCache& cache, MsHelper& ms_helper) {
  Logger::Info << "Precalculating MF weights for "
               << settings.weightMode.ToString() << " weighting...\n";
  std::unique_ptr<ImageWeights> weights = cache.MakeEmptyWeights();
  for (const ImagingTable::Group& group : imaging_table.SquaredGroups()) {
    const ImagingTableEntry& entry = *group.front();
    std::vector<MsListItem> ms_list = ms_helper.InitializeMsList(entry);
    for (MsListItem& item : ms_list) {
      weights->Grid(*item.ms_description->GetProvider());
    }
  }
  weights->FinishGridding();
  cache.SetMFWeights(std::move(weights));
  if (settings.isWeightImageSaved)
    cache.GetMFWeights()->Save(settings.prefixName + "-weights.fits");
}

}  // namespace wsclean
