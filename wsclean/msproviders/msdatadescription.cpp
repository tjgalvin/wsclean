#include "msdatadescription.h"

#include "contiguousms.h"

#include "../wsclean/wscleansettings.h"

#include "../serializable.h"

std::unique_ptr<MSProvider> MSDataDescription::GetProvider() const
{
	if(_isPartitioned)
		return std::unique_ptr<MSProvider>(new PartitionedMS(_partitionHandle, _partIndex, _polarization, _dataDescId));
	else
		return std::unique_ptr<MSProvider>(new ContiguousMS(_filename, _dataColumnName, _selection, _polarization, _dataDescId));
}

void MSDataDescription::Serialize(std::ostream& stream) const
{
	Serializable::SerializeToBool(stream, _isPartitioned);
	Serializable::SerializeToUInt16(stream, _polarization);
	Serializable::SerializeToUInt32(stream, _dataDescId);
	_selection.Serialize(stream);
	Serializable::SerializeToString(stream, _filename);
	Serializable::SerializeToString(stream, _dataColumnName);
	_partitionHandle.Serialize(stream);
	Serializable::SerializeToUInt64(stream, _partIndex);
}

std::unique_ptr<MSDataDescription> MSDataDescription::Unserialize(std::istream& stream)
{
	std::unique_ptr<MSDataDescription> mdd(new MSDataDescription());
	mdd->_isPartitioned = Serializable::UnserializeBool(stream);
	mdd->_polarization = (PolarizationEnum) Serializable::UnserializeUInt16(stream);
	mdd->_dataDescId = Serializable::UnserializeUInt32(stream);
	mdd->_selection.Unserialize(stream);
	mdd->_filename = Serializable::UnserializeString(stream);
	mdd->_dataColumnName = Serializable::UnserializeString(stream);
	mdd->_partitionHandle.Unserialize(stream);
	mdd->_partIndex = Serializable::UnserializeUInt64(stream);
	return mdd;
}
