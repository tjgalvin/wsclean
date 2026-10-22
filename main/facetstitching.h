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

#include <schaapcommon/facets/facet.h>
#include <schaapcommon/facets/facetimage.h>

namespace wsclean {

/**
 * Create the feathered mask of a facet: the mask of the facet convolved with a
 * tophat kernel of the given radius.
 *
 * The convolution is performed with FFTs, which leaves rounding noise of either
 * sign in the mask. That noise is removed, because the masks of the facets are
 * summed to normalize the stitched image, and noise of opposite signs can
 * cancel out the real coverage of a pixel. See
 * @ref wsclean::tophat_convolution::RemoveMaskNoise().
 *
 * @param facetImage The facet image. The facet must have been set with
 *        SetFacet() before calling this function.
 * @param facet The facet to create the mask for.
 * @param feather_size The radius of the tophat kernel in pixels.
 * @return The feathered mask, with the size of the facet's trimmed bounding
 *         box.
 */
aocommon::Image CreateFeatheredFacetMask(const schaapcommon::facets::FacetImage& facetImage,
                                         const schaapcommon::facets::Facet& facet,
                                         size_t feather_size);

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
