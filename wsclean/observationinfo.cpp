#include "observationinfo.h"

#include "../serialostream.h"
#include "../serialistream.h"

void ObservationInfo::Serialize(SerialOStream& stream) const {
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

void ObservationInfo::Unserialize(SerialIStream& stream) {
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
