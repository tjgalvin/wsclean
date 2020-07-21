#include "msselection.h"

#include "serialostream.h"
#include "serialistream.h"

void MSSelection::Serialize(SerialOStream& stream) const {
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

void MSSelection::Unserialize(SerialIStream& stream) {
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
