#ifndef FITS_ATERM_H
#define FITS_ATERM_H

#include "fitsatermbase.h"

#include "../fitsreader.h"

#include <aocommon/uvector.h>

#include <complex>
#include <map>
#include <memory>
#include <vector>

/**
 * Class that reads in FITS images and resamples them onto aterm grids.
 * The fits file is supposed to have a TIME, FREQ and ANTENNA axis.
 */
class FitsATerm final : public FitsATermBase {
 public:
  FitsATerm(size_t nAntenna, const CoordinateSystem& coordinateSystem);
  ~FitsATerm();

  void OpenTECFiles(const std::vector<std::string>& filenames);
  void OpenDiagGainFiles(const std::vector<std::string>& filenames);

  virtual bool Calculate(std::complex<float>* buffer, double time,
                         double frequency, size_t fieldId,
                         const double* uvwInM) override;

 private:
  enum Mode { TECMode, DiagonalMode } _mode;

  void readImages(std::complex<float>* buffer, size_t timeIndex,
                  double frequency);

  void resample(const FitsReader& reader, double* dest, const double* source);

  void evaluateTEC(std::complex<float>* dest, const double* source,
                   double frequency);

  void copyToRealPolarization(std::complex<float>* dest, const double* source,
                              size_t polIndex);
  void copyToImaginaryPolarization(std::complex<float>* dest,
                                   const double* source, size_t polIndex);
  void setPolarization(std::complex<float>* dest, size_t polIndex,
                       std::complex<float> value);

  aocommon::UVector<double> _scratchA, _scratchB;
  std::vector<FitsReader> _readers;
};

#endif
