#include "griddingresult.h"

#include "../serializable.h"

void GriddingResult::Serialize(std::ostream& stream) const
{
	imageRealResult.Serialize(stream);
	imageImaginaryResult.Serialize(stream);
	ObservationInfo observationInfo;
	Serializable::SerializeToDouble(stream, beamSize);
	Serializable::SerializeToDouble(stream, imageWeight);
	Serializable::SerializeToDouble(stream, normalizationFactor);
	Serializable::SerializeToUInt64(stream, actualWGridSize);
	Serializable::SerializeToUInt64(stream, griddedVisibilityCount);
	Serializable::SerializeToDouble(stream, effectiveGriddedVisibilityCount);
	Serializable::SerializeToDouble(stream, visibilityWeightSum);
	Serializable::SerializeToUInt64(stream, actualInversionWidth);
	Serializable::SerializeToUInt64(stream, actualInversionHeight);
	Serializable::SerializePtr(stream, cache);
}

void GriddingResult::Unserialize(std::istream& stream)
{
	imageRealResult.Unserialize(stream);
	imageImaginaryResult.Unserialize(stream);
	ObservationInfo observationInfo;
	beamSize = Serializable::UnserializeDouble(stream);
	imageWeight = Serializable::UnserializeDouble(stream);
	normalizationFactor = Serializable::UnserializeDouble(stream);
	actualWGridSize = Serializable::UnserializeUInt64(stream);
	griddedVisibilityCount = Serializable::UnserializeUInt64(stream);
	effectiveGriddedVisibilityCount = Serializable::UnserializeDouble(stream);
	visibilityWeightSum = Serializable::UnserializeDouble(stream);
	actualInversionWidth = Serializable::UnserializeUInt64(stream);
	actualInversionHeight = Serializable::UnserializeUInt64(stream);
	Serializable::UnserializePtr(stream, cache);
}
