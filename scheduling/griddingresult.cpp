#include "griddingresult.h"

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

#include "../idg/averagebeam.h"

GriddingResult::GriddingResult()
    : startTime(0.0),
      beamSize(0.0),
      imageWeight(0.0),
      normalizationFactor(0.0),
      actualWGridSize(0),
      griddedVisibilityCount(0),
      effectiveGriddedVisibilityCount(0),
      visibilityWeightSum(0),
      cache() {}

GriddingResult::GriddingResult(GriddingResult&& source) noexcept = default;
GriddingResult::~GriddingResult() = default;
GriddingResult& GriddingResult::operator=(GriddingResult&& rhs) noexcept =
    default;

void GriddingResult::Serialize(aocommon::SerialOStream& stream) const {
  stream.ObjectVector(images)
      .Double(beamSize)
      .Double(imageWeight)
      .Double(normalizationFactor)
      .UInt64(actualWGridSize)
      .UInt64(griddedVisibilityCount)
      .Double(effectiveGriddedVisibilityCount)
      .Double(visibilityWeightSum)
      .Ptr(cache)
      .Ptr(averageBeam);
}

void GriddingResult::Unserialize(aocommon::SerialIStream& stream) {
  stream.ObjectVector(images)
      .Double(beamSize)
      .Double(imageWeight)
      .Double(normalizationFactor)
      .UInt64(actualWGridSize)
      .UInt64(griddedVisibilityCount)
      .Double(effectiveGriddedVisibilityCount)
      .Double(visibilityWeightSum)
      .Ptr(cache)
      .Ptr(averageBeam);
}
