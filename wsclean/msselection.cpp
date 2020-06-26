#include "msselection.h"

#include "serializable.h"

void MSSelection::Serialize(std::ostream& stream) const
{
	Serializable::SerializeVectorUInt64(stream, _fieldIds);
	Serializable::SerializeToUInt64(stream, _bandId);
	Serializable::SerializeToUInt64(stream, _startChannel);
	Serializable::SerializeToUInt64(stream, _endChannel);
	Serializable::SerializeToUInt64(stream, _startTimestep);
	Serializable::SerializeToUInt64(stream, _endTimestep);
	Serializable::SerializeToDouble(stream, _minUVWInM);
	Serializable::SerializeToDouble(stream, _maxUVWInM);
	Serializable::SerializeToBool(stream, _autoCorrelations);
	Serializable::SerializeToUInt32(stream, _evenOddSelection);
}

void MSSelection::Unserialize(std::istream& stream)
{
	Serializable::UnserializeVectorUInt64(stream, _fieldIds);
	_bandId = Serializable::UnserializeUInt64(stream);
	_startChannel = Serializable::UnserializeUInt64(stream);
	_endChannel = Serializable::UnserializeUInt64(stream);
	_startTimestep = Serializable::UnserializeUInt64(stream);
	_endTimestep = Serializable::UnserializeUInt64(stream);
	_minUVWInM = Serializable::UnserializeDouble(stream);
	_maxUVWInM = Serializable::UnserializeDouble(stream);
	_autoCorrelations = Serializable::UnserializeBool(stream);
	_evenOddSelection = (EvenOddSelection) Serializable::UnserializeUInt32(stream);
}
