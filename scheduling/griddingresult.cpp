#include "griddingresult.h"

#include "../serialostream.h"
#include "../serialistream.h"

void GriddingResult::Serialize(SerialOStream& stream) const {
  imageRealResult.Serialize(stream);
  imageImaginaryResult.Serialize(stream);
  observationInfo.Serialize(stream);
  stream.Double(beamSize)
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

void GriddingResult::Unserialize(SerialIStream& stream) {
  imageRealResult.Unserialize(stream);
  imageImaginaryResult.Unserialize(stream);
  observationInfo.Unserialize(stream);
  stream.Double(beamSize)
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
