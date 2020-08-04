#include "everybeamaterm.h"
#include "../mwa/findcoefffile.h"

#include <EveryBeam/options.h>
#include <EveryBeam/elementresponse.h>
#include <EveryBeam/griddedresponse/griddedresponse.h>

#include <memory>

using everybeam::Load;
using everybeam::griddedresponse::GriddedResponse;

EveryBeamATerm::EveryBeamATerm(casacore::MeasurementSet& ms,
                               const CoordinateSystem& coordinateSystem,
                               const everybeam::Options& options)
    : _telescope(Load(ms, options)),
      _coordinate_system(
          {coordinateSystem.width, coordinateSystem.height, coordinateSystem.ra,
           coordinateSystem.dec, coordinateSystem.dl, coordinateSystem.dm,
           coordinateSystem.phaseCentreDL, coordinateSystem.phaseCentreDM}) {}

bool EveryBeamATerm::calculateBeam(std::complex<float>* buffer, double time,
                                   double frequency, size_t fieldId) {
  if (!_telescope->GetIsTimeRelevant()) {
    if (fieldId == _cachedFieldId && _cachedFrequency == frequency) {
      // Exit calculation
      return false;
    } else {
      // Update cached values
      _cachedFieldId = fieldId;
      _cachedFrequency = frequency;
    }
  }

  // Get the gridded response
  std::unique_ptr<GriddedResponse> grid_response =
      _telescope->GetGriddedResponse(_coordinate_system);
  grid_response->CalculateAllStations(buffer, time, frequency, fieldId);

  saveATermsIfNecessary(buffer, _telescope->GetNrStations(),
                        _coordinate_system.width, _coordinate_system.height);
  return true;
}
