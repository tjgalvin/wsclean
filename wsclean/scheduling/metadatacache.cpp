#include "metadatacache.h"

#include "../serializable.h"

void MetaDataCache::Serialize ( std::ostream& stream ) const
{
	Serializable::SerializeToUInt64(stream, msDataVector.size());
	for(const Entry& entry : msDataVector)
	{
		Serializable::SerializeToDouble(stream, entry.minW);
		Serializable::SerializeToDouble(stream, entry.maxW);
		Serializable::SerializeToDouble(stream, entry.maxWWithFlags);
		Serializable::SerializeToDouble(stream, entry.maxBaselineUVW);
		Serializable::SerializeToDouble(stream, entry.maxBaselineInM);
		Serializable::SerializeToDouble(stream, entry.integrationTime);
	}
	
	Serializable::SerializePtr(stream, averageBeam);
}

void MetaDataCache::Unserialize ( std::istream& stream )
{
	msDataVector.resize(Serializable::UnserializeUInt64(stream));
	for(Entry& entry : msDataVector)
	{
		entry.minW = Serializable::UnserializeDouble(stream);
		entry.maxW = Serializable::UnserializeDouble(stream);
		entry.maxWWithFlags = Serializable::UnserializeDouble(stream);
		entry.maxBaselineUVW = Serializable::UnserializeDouble(stream);
		entry.maxBaselineInM = Serializable::UnserializeDouble(stream);
		entry.integrationTime = Serializable::UnserializeDouble(stream);
	}
	
	Serializable::UnserializePtr(stream, averageBeam);
}
