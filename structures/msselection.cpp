#include "msselection.h"

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

#include <limits>

const size_t MSSelection::ALL_FIELDS = std::numeric_limits<size_t>::max();

void MSSelection::Serialize(aocommon::SerialOStream& stream) const {
  stream.VectorUInt64(_fieldIds)
      .UInt64(_bandId)
      .UInt64(_startChannel)
      .UInt64(_endChannel)
      .UInt64(_startTimestep)
      .UInt64(_endTimestep)
      .Double(_minUVWInM)
      .Double(_maxUVWInM)
      .Bool(_autoCorrelations)
      .UInt32(_evenOddSelection);
}

void MSSelection::Unserialize(aocommon::SerialIStream& stream) {
  stream.VectorUInt64(_fieldIds)
      .UInt64(_bandId)
      .UInt64(_startChannel)
      .UInt64(_endChannel)
      .UInt64(_startTimestep)
      .UInt64(_endTimestep)
      .Double(_minUVWInM)
      .Double(_maxUVWInM)
      .Bool(_autoCorrelations)
      .UInt32(_evenOddSelection);
}
