#include "griddingresult.h"

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

void GriddingResult::Serialize(aocommon::SerialOStream& stream) const {
  stream.ObjectVector(images)
      .Double(beamSize)
      .Double(imageWeight)
      .Double(normalizationFactor)
      .UInt64(actualWGridSize)
      .UInt64(griddedVisibilityCount)
      .Double(effectiveGriddedVisibilityCount)
      .Double(visibilityWeightSum)
      .UInt64(actualInversionWidth)
      .UInt64(actualInversionHeight)
      .Ptr(cache);
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
      .UInt64(actualInversionWidth)
      .UInt64(actualInversionHeight)
      .Ptr(cache);
}
