#include "dishaterm.h"
#include "../wsclean/logger.h"

#include <EveryBeam/griddedresponse/griddedresponse.h>

#include <limits>

#include <casacore/tables/Tables/ArrayColumn.h>

using everybeam::Load;
using everybeam::coords::CoordinateSystem;
using everybeam::griddedresponse::GriddedResponse;

DishATerm::DishATerm(casacore::MeasurementSet& ms,
                     const CoordinateSystem& coordinateSystem)
    : _coordinate_system(
          {coordinateSystem.width, coordinateSystem.height, coordinateSystem.ra,
           coordinateSystem.dec, coordinateSystem.dl, coordinateSystem.dm,
           coordinateSystem.phaseCentreDL, coordinateSystem.phaseCentreDM}) {
  // Load telescope
  _telescope = Load(ms);
}

bool DishATerm::calculateBeam(std::complex<float>* buffer, double,
                              double frequency, size_t fieldId) {
  if (fieldId == _cachedFieldId && _cachedFrequency == frequency) {
    return false;
  } else {
    Logger::Debug << "Calculating VLA beam for field " << fieldId
                  << ", frequency " << frequency * 1e-6 << " MHz.\n";
    // Get the gridded response
    std::unique_ptr<GriddedResponse> grid_response =
        _telescope->GetGriddedResponse(_coordinate_system);
    grid_response->CalculateAllStations(buffer, 0., frequency, fieldId);

    _cachedFieldId = fieldId;
    _cachedFrequency = frequency;

    saveATermsIfNecessary(buffer, _telescope->GetNrStations(),
                          _coordinate_system.width, _coordinate_system.height);
    return true;
  }
}
