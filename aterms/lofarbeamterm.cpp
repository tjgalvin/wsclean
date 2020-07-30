#include "lofarbeamterm.h"

#include <EveryBeam/options.h>
#include <EveryBeam/elementresponse.h>
#include <EveryBeam/griddedresponse/griddedresponse.h>

using everybeam::ElementResponseModel;
using everybeam::Load;
using everybeam::Options;
using everybeam::coords::CoordinateSystem;
using everybeam::griddedresponse::GriddedResponse;

LofarBeamTerm::LofarBeamTerm(casacore::MeasurementSet& ms,
                             const CoordinateSystem& coordinateSystem,
                             const std::string& dataColumnName)
    : _coordinate_system(
          {coordinateSystem.width, coordinateSystem.height, coordinateSystem.ra,
           coordinateSystem.dec, coordinateSystem.dl, coordinateSystem.dm,
           coordinateSystem.phaseCentreDL, coordinateSystem.phaseCentreDM}),
      _useDifferentialBeam(false),
      _useChannelFrequency(true) {
  // Telescope options
  Options options;
  options.use_differential_beam = _useDifferentialBeam;
  options.use_channel_frequency = _useChannelFrequency;
  options.data_column_name = dataColumnName;

  // Response model
  ElementResponseModel response_model = ElementResponseModel::kHamaker;

  // Load telescope
  _telescope = Load(ms, options, response_model);
}

bool LofarBeamTerm::calculateBeam(std::complex<float>* buffer, double time,
                                  double frequency, size_t) {
  // Get the gridded response
  std::unique_ptr<GriddedResponse> grid_response =
      _telescope->GetGriddedResponse(_coordinate_system);
  grid_response->CalculateAllStations(buffer, time, frequency, 0);

  saveATermsIfNecessary(buffer, _telescope->GetNrStations(),
                        _coordinate_system.width, _coordinate_system.height);
  return true;
}