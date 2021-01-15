#include "observationinfo.h"

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

void ObservationInfo::Serialize(aocommon::SerialOStream& stream) const {
  stream.Double(phaseCentreRA)
      .Double(phaseCentreDec)
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
      .Bool(hasDenormalPhaseCentre)
      .Double(phaseCentreDL)
      .Double(phaseCentreDM)
      .String(telescopeName)
      .String(observer)
      .String(fieldName);
}

struct ObservationInfo ReadObservationInfo(casacore::MeasurementSet& ms,
                                           size_t fieldId) {
  struct ObservationInfo obsInfo;

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
  if (fTable.keywordSet().isDefined("WSCLEAN_DL")) {
    obsInfo.phaseCentreDL =
        fTable.keywordSet().asDouble(casacore::RecordFieldId("WSCLEAN_DL"));
  } else {
    obsInfo.phaseCentreDL = 0.0;
  }
  if (fTable.keywordSet().isDefined("WSCLEAN_DM")) {
    obsInfo.phaseCentreDM =
        fTable.keywordSet().asDouble(casacore::RecordFieldId("WSCLEAN_DM"));
  } else {
    obsInfo.phaseCentreDM = 0.0;
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
