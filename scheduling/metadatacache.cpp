#include "metadatacache.h"

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

void MetaDataCache::Serialize(aocommon::SerialOStream& stream) const {
  stream.UInt64(msDataVector.size());
  for (const Entry& entry : msDataVector) {
    stream.Double(entry.minW)
        .Double(entry.maxW)
        .Double(entry.maxWWithFlags)
        .Double(entry.maxBaselineUVW)
        .Double(entry.maxBaselineInM)
        .Double(entry.integrationTime);
  }

  stream.Ptr(averageBeam).Float(h5Sum).Float(correctionSum);
}

void MetaDataCache::Unserialize(aocommon::SerialIStream& stream) {
  msDataVector.resize(stream.UInt64());
  for (Entry& entry : msDataVector) {
    stream.Double(entry.minW)
        .Double(entry.maxW)
        .Double(entry.maxWWithFlags)
        .Double(entry.maxBaselineUVW)
        .Double(entry.maxBaselineInM)
        .Double(entry.integrationTime);
  }

  stream.Ptr(averageBeam).Float(h5Sum).Float(correctionSum);
}
