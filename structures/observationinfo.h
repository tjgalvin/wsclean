#ifndef OBSERVATION_INFO_H
#define OBSERVATION_INFO_H

#include <aocommon/io/serialstreamfwd.h>

#include <string>

struct ObservationInfo {
  double phaseCentreRA = 0.0, phaseCentreDec = 0.0;
  double startTime = 0.0;
  bool hasDenormalPhaseCentre = false;
  double phaseCentreDL = 0.0, phaseCentreDM = 0.0;
  std::string telescopeName;
  std::string observer;
  std::string fieldName;

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);
};

#endif
