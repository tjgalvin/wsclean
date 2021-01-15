#ifndef OBSERVATION_INFO_H
#define OBSERVATION_INFO_H

#include <aocommon/io/serialstreamfwd.h>

#include <string>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/tables/Tables/TableRecord.h>

#include "../structures/msselection.h"

struct ObservationInfo {
  double phaseCentreRA = 0.0, phaseCentreDec = 0.0;
  bool hasDenormalPhaseCentre = false;
  double phaseCentreDL = 0.0, phaseCentreDM = 0.0;
  std::string telescopeName;
  std::string observer;
  std::string fieldName;

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);
};

/// Generates observation info from the measurement set tables.
struct ObservationInfo ReadObservationInfo(casacore::MeasurementSet& ms,
                                           size_t fieldId);

#endif
