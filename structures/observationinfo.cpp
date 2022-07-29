#include "observationinfo.h"

#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/tables/Tables/TableRecord.h>

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

#include "../structures/msselection.h"

void ObservationInfo::Serialize(aocommon::SerialOStream& stream) const {
  stream.Double(phaseCentreRA)
      .Double(phaseCentreDec)
      .String(telescopeName)
      .String(observer)
      .String(fieldName);
}

void ObservationInfo::Unserialize(aocommon::SerialIStream& stream) {
  stream.Double(phaseCentreRA)
      .Double(phaseCentreDec)
      .String(telescopeName)
      .String(observer)
      .String(fieldName);
}

ObservationInfo ReadObservationInfo(casacore::MeasurementSet& ms,
                                    size_t fieldId) {
  ObservationInfo obsInfo;

  casacore::MSAntenna aTable = ms.antenna();
  size_t antennaCount = aTable.nrow();
  if (antennaCount == 0) throw std::runtime_error("No antennae in set");
  casacore::MPosition::ScalarColumn antPosColumn(
      aTable, aTable.columnName(casacore::MSAntennaEnums::POSITION));
  casacore::MPosition ant1Pos = antPosColumn(0);

  size_t fieldRow = (fieldId == MSSelection::ALL_FIELDS) ? 0 : fieldId;
  casacore::MEpoch::ScalarColumn timeColumn(
      ms, ms.columnName(casacore::MSMainEnums::TIME));
  casacore::MSField fTable(ms.field());
  casacore::MDirection::ScalarColumn phaseDirColumn(
      fTable, fTable.columnName(casacore::MSFieldEnums::PHASE_DIR));
  casacore::MDirection phaseDir = phaseDirColumn(fieldRow);
  casacore::MEpoch curtime = timeColumn(0);
  casacore::MeasFrame frame(ant1Pos, curtime);
  casacore::MDirection::Ref j2000Ref(casacore::MDirection::J2000, frame);
  casacore::MDirection j2000 =
      casacore::MDirection::Convert(phaseDir, j2000Ref)();
  casacore::Vector<casacore::Double> j2000Val = j2000.getValue().get();

  obsInfo.phaseCentreRA = j2000Val[0];
  obsInfo.phaseCentreDec = j2000Val[1];

  if (fTable.keywordSet().isDefined("WSCLEAN_DL") ||
      fTable.keywordSet().isDefined("WSCLEAN_DM")) {
    throw std::runtime_error(
        "This measurement was processed with chgcentre to have a shifted phase "
        "centre. \n"
        "This kind of shifting has been directly implemented in wsclean. "
        "Because of this,\n"
        "the chgcentre approach is no longer supported. Use chgcentre to undo "
        "the shift\n"
        "and use wsclean's -shift parameter.\n");
  }

  casacore::MSObservation oTable = ms.observation();
  size_t obsCount = oTable.nrow();
  if (obsCount == 0) throw std::runtime_error("No observations in set");
  casacore::ScalarColumn<casacore::String> telescopeNameColumn(
      oTable, oTable.columnName(casacore::MSObservation::TELESCOPE_NAME));
  casacore::ScalarColumn<casacore::String> observerColumn(
      oTable, oTable.columnName(casacore::MSObservation::OBSERVER));
  obsInfo.telescopeName = telescopeNameColumn(0);
  obsInfo.observer = observerColumn(0);

  casacore::ScalarColumn<casacore::String> fieldNameColumn(
      fTable, fTable.columnName(casacore::MSField::NAME));
  obsInfo.fieldName = fieldNameColumn(fieldRow);

  return obsInfo;
}
