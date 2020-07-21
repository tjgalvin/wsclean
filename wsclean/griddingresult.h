#ifndef GRIDDING_RESULT_H
#define GRIDDING_RESULT_H

#include "imagebufferallocator.h"
#include "observationinfo.h"

#include <string>

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

  ImageBufferAllocator::Ptr imageRealResult;
  ImageBufferAllocator::Ptr imageImaginaryResult;
  ObservationInfo observationInfo;
  double beamSize;
  double imageWeight;
  double normalizationFactor;
  size_t actualWGridSize;
  size_t griddedVisibilityCount;
  double effectiveGriddedVisibilityCount;
  double visibilityWeightSum;
  size_t actualInversionWidth, actualInversionHeight;
};

#endif
