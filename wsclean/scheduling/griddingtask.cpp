#include "griddingtask.h"

#include "../serializable.h"

void GriddingTask::Serialize(std::ostream& stream) const
{
	Serializable::SerializeToUInt32(stream, operation);
	Serializable::SerializeToBool(stream, imagePSF);
	Serializable::SerializeToBool(stream, subtractModel);
	Serializable::SerializeToUInt32(stream, polarization);
	Serializable::SerializeToBool(stream, verbose);
	Serializable::SerializePtr(stream, cache);
	Serializable::SerializeToBool(stream, storeImagingWeights);
	
	Serializable::SerializePtr(stream, imageWeights);
	
	// msList
	Serializable::SerializeToUInt64(stream, msList.size());
	for(const std::unique_ptr<MSDataDescription>& dataDesc : msList)
		dataDesc->Serialize(stream);
	
	Serializable::SerializeToBool(stream, addToModel);
	modelImageReal.Serialize(stream);
	modelImageImaginary.Serialize(stream);
}

void GriddingTask::Unserialize(std::istream& stream)
{
	operation = (Operation) Serializable::UnserializeUInt32(stream);
	imagePSF = Serializable::UnserializeBool(stream);
	subtractModel = Serializable::UnserializeBool(stream);
	polarization = (PolarizationEnum) Serializable::UnserializeUInt32(stream);
	verbose = Serializable::UnserializeBool(stream);
	Serializable::UnserializePtr(stream, cache);
	storeImagingWeights = Serializable::UnserializeBool(stream);
	
	Serializable::UnserializePtr(stream, imageWeights);
	
	// msList
	msList.resize(Serializable::UnserializeUInt64(stream));
	for(std::unique_ptr<MSDataDescription>& dataDesc : msList)
		dataDesc = MSDataDescription::Unserialize(stream);
	
	addToModel = Serializable::UnserializeBool(stream);
	modelImageReal.Unserialize(stream);
	modelImageImaginary.Unserialize(stream);
}

