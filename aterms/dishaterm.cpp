#include "dishaterm.h"

#include "../primarybeam/vlabeam.h"
#include "../primarybeam/voltagepattern.h"

#include "../wsclean/logger.h"

#include <limits>

#include <casacore/tables/Tables/ArrayColumn.h>

DishATerm::DishATerm(casacore::MeasurementSet& ms,
                     const CoordinateSystem& coordinateSystem)
    : _width(coordinateSystem.width),
      _height(coordinateSystem.height),
      _phaseCentreRA(coordinateSystem.ra),
      _phaseCentreDec(coordinateSystem.dec),
      _dl(coordinateSystem.dl),
      _dm(coordinateSystem.dm),
      _phaseCentreDL(coordinateSystem.phaseCentreDL),
      _phaseCentreDM(coordinateSystem.phaseCentreDM),
      _cachedFieldId(std::numeric_limits<size_t>::max()),
      _cachedFrequency(0.0) {
  _nAntenna = ms.antenna().nrow();
  casacore::MSField fieldTable = ms.field();
  casacore::ArrayColumn<double> pointingDirCol(
      fieldTable, casacore::MSField::columnName(casacore::MSField::DELAY_DIR));
  for (size_t fieldId = 0; fieldId != fieldTable.nrow(); ++fieldId) {
    casacore::Array<double> pDir = pointingDirCol(fieldId);
    double pDirRA = *pDir.cbegin();
    double pDirDec = *(pDir.cbegin() + 1);
    _fieldPointing.emplace_back(pDirRA, pDirDec);
  }
}

bool DishATerm::calculateBeam(std::complex<float>* buffer, double,
                              double frequency, size_t fieldId) {
  if (fieldId == _cachedFieldId && _cachedFrequency == frequency) {
    return false;
  } else {
    Logger::Debug << "Calculating VLA beam for field " << fieldId
                  << ", frequency " << frequency * 1e-6 << " MHz.\n";
    double pDirRA = _fieldPointing[fieldId].first,
           pDirDec = _fieldPointing[fieldId].second;
    std::array<double, 5> coefs = VLABeam::GetCoefficients("", frequency);
    VoltagePattern vp;
    vp.frequencies.assign(1, frequency);
    vp.maximumRadiusArcMin = 53.0;
    aocommon::UVector<double> coefsVec(coefs.begin(), coefs.end());
    vp.EvaluatePolynomial(coefsVec, false);
    vp.Render(buffer, _width, _height, _dl, _dm, _phaseCentreRA,
              _phaseCentreDec, pDirRA, pDirDec, _phaseCentreDL, _phaseCentreDM,
              frequency);

    _cachedFieldId = fieldId;
    _cachedFrequency = frequency;

    for (size_t i = 1; i != _nAntenna; ++i) {
      std::copy_n(buffer, _width * _height * 4,
                  buffer + i * _width * _height * 4);
    }

    saveATermsIfNecessary(buffer, _nAntenna, _width, _height);

    return true;
  }
}
