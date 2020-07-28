#include "lofarbeamterm.h"

#include <aocommon/imagecoordinates.h>

#include <EveryBeam/options.h>
#include <EveryBeam/element_response.h>

#include <EveryBeam/coords/ITRFDirection.h>
#include <EveryBeam/coords/ITRFConverter.h>
#include <EveryBeam/griddedresponse/griddedresponse.h>

#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/tables/Tables/TableRecord.h>

#include <algorithm>

using namespace everybeam;
using namespace aocommon;

LofarBeamTerm::LofarBeamTerm(casacore::MeasurementSet& ms,
                             const CoordinateSystem& coordinateSystem,
                             const std::string& dataColumnName)
    : _useDifferentialBeam(false), _useChannelFrequency(true) {
  // Telescope options
  Options options;
  options.use_differential_beam = _useDifferentialBeam;
  options.use_channel_frequency = _useChannelFrequency;
  options.data_column_name = dataColumnName;

  // Response model
  ElementResponseModel response_model = ElementResponseModel::kHamaker;

  // Load telescope
  _telescope = Load(ms, options, response_model);

  casacore::MSAntenna aTable(ms.antenna());
  casacore::MSField fieldTable(ms.field());
  if (fieldTable.nrow() != 1)
    throw std::runtime_error("Set has multiple fields");

  casacore::MPosition::ScalarColumn antPosColumn(
      aTable, aTable.columnName(casacore::MSAntennaEnums::POSITION));
  casacore::MPosition arrayPos = antPosColumn(0);

  // Compute Ra/Dec pointing direction
  casacore::MEpoch::ScalarColumn timeColumn(
      ms, ms.columnName(casacore::MSMainEnums::TIME));
  casacore::MDirection::ScalarColumn phaseDirColumn(
      fieldTable, fieldTable.columnName(casacore::MSFieldEnums::PHASE_DIR));
  casacore::MDirection phaseDir = phaseDirColumn(0);
  casacore::MEpoch curtime = timeColumn(0);
  casacore::MeasFrame frame(arrayPos, curtime);
  casacore::MDirection::Ref j2000Ref(casacore::MDirection::J2000, frame);
  casacore::MDirection j2000 =
      casacore::MDirection::Convert(phaseDir, j2000Ref)();
  casacore::Vector<casacore::Double> j2000Val = j2000.getValue().get();
  double phaseCentreRA = j2000Val[0];
  double phaseCentreDec = j2000Val[1];

  // Convert to everybeam::CoordinateSystem
  _coordinate_system = {coordinateSystem.width,
                        coordinateSystem.height,
                        phaseCentreRA,
                        phaseCentreDec,
                        coordinateSystem.dl,
                        coordinateSystem.dm,
                        coordinateSystem.phaseCentreDL,
                        coordinateSystem.phaseCentreDM};
}

bool LofarBeamTerm::calculateBeam(std::complex<float>* buffer, double time,
                                  double frequency, size_t) {
  // Get the gridded response
  std::unique_ptr<griddedresponse::GriddedResponse> grid_response =
      _telescope->GetGriddedResponse(_coordinate_system);
  grid_response->CalculateAllStations(buffer, time, frequency, 0);

  saveATermsIfNecessary(buffer, _telescope->GetNrStations(),
                        _coordinate_system.width, _coordinate_system.height);
  return true;
}