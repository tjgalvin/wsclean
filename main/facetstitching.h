#ifndef FACET_STITCHING_H
#define FACET_STITCHING_H

#include "../io/cachedimageset.h"
#include "../io/wscfitswriter.h"
#include "../structures/imagingtable.h"
#include "../structures/outputchannelinfo.h"

#include "settings.h"
#include "stitchingmemory.h"

#include <functional>
#include <vector>

namespace wsclean {

struct StitchFacetsConfig {
  size_t facet_count = 0;
  const Settings* settings = nullptr;
  const std::vector<OutputChannelInfo>* info_per_channel = nullptr;
  std::function<void(const ImagingTableEntry&, aocommon::PolarizationEnum, size_t)>
      initialize_model_images;
  std::function<void(aocommon::Image&, const ImagingTableEntry&)> process_full_psf;
  std::function<WSCFitsWriter(const ImagingTableEntry&, bool, bool)> create_writer;
  std::function<WSCFitsWriter(const ImagingTableEntry&, aocommon::PolarizationEnum, bool, bool)>
      create_writer_for_pol;
};

void StitchFacets(const ImagingTable& table, CachedImageSet& image_cache, StitchImageType type,
                  const StitchFacetsConfig& config);

}  // namespace wsclean

#endif
