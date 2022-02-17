#ifndef GRIDDING_RESULT_H
#define GRIDDING_RESULT_H

#include "../scheduling/metadatacache.h"

#include <aocommon/image.h>
#include <aocommon/io/serialstreamfwd.h>

#include <string>

class AverageBeam;

struct GriddingResult {
  GriddingResult();
  GriddingResult(const GriddingResult& source) = delete;
  GriddingResult(GriddingResult&& source) noexcept;
  ~GriddingResult();
  GriddingResult& operator=(const GriddingResult& rhs) = delete;
  GriddingResult& operator=(GriddingResult&& rhs) noexcept;

  /**
   * List of produced images. When performing complex images, images[0] will be
   * the real part and images[1] will be the imaginary part. When performing
   * full polarization imaging (indicated with Polarization::FullStokes) with
   * IDG, the list will contain all four images ordered IQUV. In all other
   * cases, this list will only hold one image.
   */
  std::vector<aocommon::Image> images;
  double startTime;
  double beamSize;
  double imageWeight;
  double normalizationFactor;
  size_t actualWGridSize;
  size_t griddedVisibilityCount;
  double effectiveGriddedVisibilityCount;
  double visibilityWeightSum;
  std::unique_ptr<MetaDataCache> cache;
  std::unique_ptr<AverageBeam> averageBeam;

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);
};

#endif
