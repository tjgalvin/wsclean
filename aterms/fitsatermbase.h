#ifndef FITS_ATERM_BASE_H
#define FITS_ATERM_BASE_H

#include "atermbase.h"
#include "atermresampler.h"
#include "cache.h"

#include "../fitsreader.h"
#include "../windowfunction.h"

class FitsATermBase : public ATermBase {
 public:
  FitsATermBase(size_t nAntenna, const CoordinateSystem& coordinateSystem);
  ~FitsATermBase();

  virtual double AverageUpdateTime() const override;

  void SetTukeyWindow(double padding) { _resampler.SetTukeyWindow(padding); }

  void SetWindow(WindowFunction::Type window) { _resampler.SetWindow(window); }

  void SetDownSample(bool downsample) { _resampler.SetDownSample(downsample); }

  ATermResampler& Resampler() { return _resampler; }

 protected:
  void initializeFromFiles(std::vector<FitsReader>& readers);

  bool findFilePosition(std::complex<float>* buffer, double time,
                        double frequency, size_t& timeindex,
                        bool& requiresRecalculation);

  void storeInCache(double frequency, const std::complex<float>* buffer);

  struct Timestep {
    double time;
    size_t readerIndex;
    size_t imgIndex;
  };
  std::vector<Timestep> _timesteps;

  size_t Width() const { return _width; }
  size_t Height() const { return _height; }
  size_t NAntenna() const { return _nAntenna; }
  size_t NFrequencies() const { return _nFrequencies; }
  double DL() const { return _dl; }
  double DM() const { return _dm; }
  double PhaseCentreDL() const { return _phaseCentreDL; }
  double PhaseCentreDM() const { return _phaseCentreDM; }

 private:
  Cache _cache;
  size_t _curTimeindex;
  double _curFrequency;
  size_t _nFrequencies, _nAntenna, _width, _height;
  double _ra, _dec, _dl, _dm, _phaseCentreDL, _phaseCentreDM;
  ATermResampler _resampler;
};

#endif
