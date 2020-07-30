#include "mwabeamterm.h"

#include <EveryBeam/griddedresponse/griddedresponse.h>
#include <EveryBeam/options.h>

#include "../wsclean/logger.h"
#include "../mwa/findcoefffile.h"

using everybeam::Load;
using everybeam::Options;
using everybeam::coords::CoordinateSystem;
using everybeam::griddedresponse::GriddedResponse;
using wsclean::mwa::FindCoeffFile;

MWABeamTerm::MWABeamTerm(casacore::MeasurementSet& ms,
                         const CoordinateSystem& coordinateSystem,
                         const std::string& search_path)
    : _coordinate_system(
          {coordinateSystem.width, coordinateSystem.height, coordinateSystem.ra,
           coordinateSystem.dec, coordinateSystem.dl, coordinateSystem.dm,
           coordinateSystem.phaseCentreDL, coordinateSystem.phaseCentreDM}),
      _frequencyInterpolation(true),
      _coeff_path(FindCoeffFile(search_path)) {
  // Telescope options
  Options options;
  options.coeff_path = _coeff_path;
  options.frequency_interpolation = _frequencyInterpolation;
  _telescope = Load(ms, options);
}

bool MWABeamTerm::calculateBeam(std::complex<float>* buffer, double time,
                                double frequency, size_t) {
  // Get the gridded response
  std::unique_ptr<GriddedResponse> grid_response =
      _telescope->GetGriddedResponse(_coordinate_system);
  grid_response->CalculateAllStations(buffer, time, frequency, 0);

  saveATermsIfNecessary(buffer, _telescope->GetNrStations(),
                        _coordinate_system.width, _coordinate_system.height);

  return true;
}
