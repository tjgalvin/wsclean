#ifndef META_DATA_CACHE_H
#define META_DATA_CACHE_H

#include <aocommon/io/serialstreamfwd.h>

#include <memory>
#include <vector>

struct MetaDataCache {
  struct Entry {
    double minW, maxW, maxWWithFlags, maxBaselineUVW, maxBaselineInM,
        integrationTime;
  };
  std::vector<Entry> msDataVector;
  float h5Sum = 0.0;
  float correctionSum = 0.0;

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);
};

#endif
