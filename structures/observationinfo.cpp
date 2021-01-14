#include "observationinfo.h"

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

void ObservationInfo::Serialize(aocommon::SerialOStream& stream) const {
  stream.Double(phaseCentreRA)
      .Double(phaseCentreDec)
      .Double(startTime)
      .Bool(hasDenormalPhaseCentre)
      .Double(phaseCentreDL)
      .Double(phaseCentreDM)
      .String(telescopeName)
      .String(observer)
      .String(fieldName);
}

void ObservationInfo::Unserialize(aocommon::SerialIStream& stream) {
  stream.Double(phaseCentreRA)
      .Double(phaseCentreDec)
      .Double(startTime)
      .Bool(hasDenormalPhaseCentre)
      .Double(phaseCentreDL)
      .Double(phaseCentreDM)
      .String(telescopeName)
      .String(observer)
      .String(fieldName);
}
