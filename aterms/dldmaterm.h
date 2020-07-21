#ifndef DLDM_ATERM_H
#define DLDM_ATERM_H

#include "fitsatermbase.h"

#include "../fitsreader.h"

#include <complex>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

class DLDMATerm final : public FitsATermBase {
 public:
  DLDMATerm(size_t nAntenna, const CoordinateSystem& coordinateSystem);

  void Open(const std::vector<std::string>& filenames);

  virtual bool Calculate(std::complex<float>* buffer, double time,
                         double frequency, size_t fieldId,
                         const double* uvwInM) override;

  void SetUpdateInterval(double updateInterval) {
    _updateInterval = updateInterval;
  }

  virtual double AverageUpdateTime() const override {
    return std::min(FitsATermBase::AverageUpdateTime(), _updateInterval);
  }

 private:
  std::vector<FitsReader> _readers;
  aocommon::UVector<double> _scratch, _dlImage, _dmImage;
  std::vector<std::array<double, 2>> _uvws;
  double _updateInterval, _previousTime;

  void readImages(std::complex<float>* buffer, size_t timeIndex,
                  double frequency, const double* uvwInM);
  void evaluateDLDM(std::complex<float>* dest, const double* dl,
                    const double* dm, const double* uvwInM);
};

#endif
