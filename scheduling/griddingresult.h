#ifndef GRIDDING_RESULT_H
#define GRIDDING_RESULT_H

#include "../wsclean/observationinfo.h"

#include <string>

#include "../scheduling/metadatacache.h"

#include "../image.h"

struct GriddingResult {
  GriddingResult()
      : beamSize(0.0),
        imageWeight(0.0),
        normalizationFactor(0.0),
        actualWGridSize(0),
        griddedVisibilityCount(0),
        effectiveGriddedVisibilityCount(0),
        visibilityWeightSum(0),
        actualInversionWidth(0),
        actualInversionHeight(0) {}

  Image imageRealResult;
  Image imageImaginaryResult;
  ObservationInfo observationInfo;
  double beamSize;
  double imageWeight;
  double normalizationFactor;
  size_t actualWGridSize;
  size_t griddedVisibilityCount;
  double effectiveGriddedVisibilityCount;
  double visibilityWeightSum;
  size_t actualInversionWidth, actualInversionHeight;
  std::unique_ptr<MetaDataCache> cache;

  void Serialize(class SerialOStream& stream) const;
  void Unserialize(class SerialIStream& stream);
};

#endif
