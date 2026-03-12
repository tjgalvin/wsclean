#include "griddingresult.h"

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

#include "../idg/averagebeam.h"

namespace wsclean {

GriddingResult::GriddingResult(GriddingResult&& source) noexcept = default;
GriddingResult::~GriddingResult() = default;
GriddingResult& GriddingResult::operator=(GriddingResult&& rhs) noexcept = default;

void GriddingResult::Serialize(aocommon::SerialOStream& stream) const {
  stream.UInt32(unique_id)
      .Double(startTime)
      .Double(beamSize)
      .UInt64(griddedVisibilityCount)
      .Double(visibilityWeightSum)
      .ObjectVector(facets);
}

void GriddingResult::Unserialize(aocommon::SerialIStream& stream) {
  stream.UInt32(unique_id)
      .Double(startTime)
      .Double(beamSize)
      .UInt64(griddedVisibilityCount)
      .Double(visibilityWeightSum)
      .ObjectVector(facets);
}

void GriddingResult::FacetData::Serialize(aocommon::SerialOStream& stream) const {
  stream.ObjectVector(images)
      .Double(imageWeight)
      .Double(normalizationFactor)
      .UInt64(actualWGridSize)
      .Double(effectiveGriddedVisibilityCount)
      .Object(averageCorrection)
      .Object(averageBeamCorrection)
      .Ptr(cache)
      .Ptr(averageBeam);
}

void GriddingResult::FacetData::Unserialize(aocommon::SerialIStream& stream) {
  stream.ObjectVector(images)
      .Double(imageWeight)
      .Double(normalizationFactor)
      .UInt64(actualWGridSize)
      .Double(effectiveGriddedVisibilityCount)
      .Object(averageCorrection)
      .Object(averageBeamCorrection)
      .Ptr(cache)
      .Ptr(averageBeam);
}

}  // namespace wsclean
