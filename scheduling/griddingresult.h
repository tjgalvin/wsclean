#ifndef GRIDDING_RESULT_H
#define GRIDDING_RESULT_H

#include "../scheduling/metadatacache.h"
#include "../structures/image.h"

#include <aocommon/io/serialstreamfwd.h>

#include <string>

struct GriddingResult {
  GriddingResult()
      : startTime(0.0),
        beamSize(0.0),
        imageWeight(0.0),
        normalizationFactor(0.0),
        actualWGridSize(0),
        griddedVisibilityCount(0),
        effectiveGriddedVisibilityCount(0),
        visibilityWeightSum(0),
        cache() {}

  /**
   * List of produced images. When performing complex images, images[0] will be
   * the real part and images[1] will be the imaginary part. When performing
   * full polarization imaging (indicated with Polarization::FullStokes) with
   * IDG, the list will contain all four images ordered IQUV. In all other
   * cases, this list will only hold one image.
   */
  std::vector<Image> images;
  double startTime;
  double beamSize;
  double imageWeight;
  double normalizationFactor;
  size_t actualWGridSize;
  size_t griddedVisibilityCount;
  double effectiveGriddedVisibilityCount;
  double visibilityWeightSum;
  std::unique_ptr<MetaDataCache> cache;

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);
};

#endif
