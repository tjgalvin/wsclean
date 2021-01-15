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
        actualInversionWidth(0),
        actualInversionHeight(0) {}

  ImageF imageRealResult;
  ImageF imageImaginaryResult;
  double startTime;
  double beamSize;
  double imageWeight;
  double normalizationFactor;
  size_t actualWGridSize;
  size_t griddedVisibilityCount;
  double effectiveGriddedVisibilityCount;
  double visibilityWeightSum;
  size_t actualInversionWidth, actualInversionHeight;
  std::unique_ptr<MetaDataCache> cache;

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);
};

#endif
