#ifndef DISH_ATERM_H
#define DISH_ATERM_H

#include <complex>

#include "atermbeam.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

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

  size_t _width, _height, _nAntenna;
  double _phaseCentreRA, _phaseCentreDec, _dl, _dm, _phaseCentreDL,
      _phaseCentreDM;
  size_t _cachedFieldId;
  double _cachedFrequency;
  std::vector<std::pair<double, double>> _fieldPointing;
};

#endif
