#ifndef OBSERVATION_INFO_H
#define OBSERVATION_INFO_H

#include <aocommon/io/serialstreamfwd.h>

#include <string>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

namespace wsclean {

struct ObservationInfo {
  double phaseCentreRA = 0.0;
  double phaseCentreDec = 0.0;
  std::string telescopeName;
  std::string observer;
  std::string fieldName;

  /// Comparison operator, for use in tests.
  bool operator==(const ObservationInfo& other) const {
    return phaseCentreRA == other.phaseCentreRA && phaseCentreDec == other.phaseCentreDec &&
           telescopeName == other.telescopeName && observer == other.observer &&
           fieldName == other.fieldName;
  }

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);
};

/// Generates observation info from the measurement set tables.
struct ObservationInfo ReadObservationInfo(casacore::MeasurementSet& ms, size_t fieldId);

}  // namespace wsclean

#endif
