#ifndef DISH_ATERM_H
#define DISH_ATERM_H

#include <complex>

#include "atermbeam.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#ifdef HAVE_EVERYBEAM
#include <EveryBeam/load.h>
#include <EveryBeam/coords/coordutils.h>

/**
 * This class calculates the a-terms for dishes with a circularly symmetric
 * response.
 */
class DishATerm : public ATermBeam {
 public:
  DishATerm(casacore::MeasurementSet& ms,
            const CoordinateSystem& coordinateSystem);

  virtual double AverageUpdateTime() const final override { return 60.0 * 5; }

 private:
  bool calculateBeam(std::complex<float>* buffer, double time, double frequency,
                     size_t fieldId) final override;

  std::unique_ptr<everybeam::telescope::Telescope> _telescope;
  size_t _cachedFieldId;
  double _cachedFrequency;
  everybeam::coords::CoordinateSystem _coordinate_system;
};
#else
using DishATerm = ATermStub;
#endif  // HAVE_EVERYBEAM

#endif
