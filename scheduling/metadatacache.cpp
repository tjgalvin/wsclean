#include "metadatacache.h"

#include "../serialostream.h"
#include "../serialistream.h"

void MetaDataCache::Serialize(class SerialOStream& stream) const {
  stream.UInt64(msDataVector.size());
  for (const Entry& entry : msDataVector) {
    stream.Double(entry.minW)
        .Double(entry.maxW)
        .Double(entry.maxWWithFlags)
        .Double(entry.maxBaselineUVW)
        .Double(entry.maxBaselineInM)
        .Double(entry.integrationTime);
  }

  stream.Ptr(averageBeam);
}

void MetaDataCache::Unserialize(class SerialIStream& stream) {
  msDataVector.resize(stream.UInt64());
  for (Entry& entry : msDataVector) {
    stream.Double(entry.minW)
        .Double(entry.maxW)
        .Double(entry.maxWWithFlags)
        .Double(entry.maxBaselineUVW)
        .Double(entry.maxBaselineInM)
        .Double(entry.integrationTime);
  }

  stream.Ptr(averageBeam);
}
