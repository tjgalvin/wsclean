#include "averagebeam.h"

#include "../serializable.h"

void AverageBeam::Serialize ( std::ostream& stream ) const
{
	if(_scalarBeam)
	{
		Serializable::SerializeToBool(stream, true);
		Serializable::SerializeVectorFloat(stream, *_scalarBeam);
	}
	else {
		Serializable::SerializeToBool(stream, false);
	}
	
	if(_matrixInverseBeam)
	{
		Serializable::SerializeToBool(stream, true);
		Serializable::SerializeVectorFloatC(stream, *_matrixInverseBeam);
	}
	else {
		Serializable::SerializeToBool(stream, false);
	}
}

void AverageBeam::Unserialize ( std::istream& stream )
{
	bool hasScalar = Serializable::UnserializeBool(stream);
	if(hasScalar)
		_scalarBeam.reset(new std::vector<float>(Serializable::UnserializeVectorFloat(stream)));
	else
		_scalarBeam.reset();
	
	bool hasMatrixInverse = Serializable::UnserializeBool(stream);
	if(hasMatrixInverse)
		_matrixInverseBeam.reset(new std::vector<std::complex<float>>(Serializable::UnserializeVectorFloatC(stream)));
	else
		_matrixInverseBeam.reset();
}
