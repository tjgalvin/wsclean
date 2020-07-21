#ifndef ATERM_RESAMPLER_H
#define ATERM_RESAMPLER_H

#include "atermbase.h"

#include "../fitsreader.h"
#include "../windowfunction.h"

#include <aocommon/uvector.h>

class ATermResampler {
 public:
  ATermResampler(const ATermBase::CoordinateSystem& coordinateSystem);
  ~ATermResampler();

  /**
   * @param scratch vector of size at least ScratchASize()
   * @param output vector of size at least ScratchBSize()
   */
  void ReadAndResample(FitsReader& reader, size_t fileIndex,
                       aocommon::UVector<double>& scratch,
                       aocommon::UVector<double>& output, double stretchFactor);

  void SetTukeyWindow(double padding) {
    _window = WindowFunction::Tukey;
    _padding = padding;
  }

  void SetWindow(WindowFunction::Type window) { _window = window; }

  void SetDownSample(bool downsample) { _downsample = downsample; }

  size_t AllocatedWidth() const { return _allocatedWidth; }
  size_t AllocatedHeight() const { return _allocatedHeight; }

  size_t ScratchASize() const { return _allocatedWidth * _allocatedHeight; }

  size_t ScratchBSize(const FitsReader& reader) const {
    return std::max(_coordinateSystem.width * _coordinateSystem.height,
                    reader.ImageWidth() * reader.ImageHeight());
  }

  void OverrideFitsPhaseCentre(double ra, double dec) {
    _overrideFitsPhaseCentre = true;
    _overrideRa = ra;
    _overrideDec = dec;
  }

 private:
  void regrid(const FitsReader& reader, double* dest, const double* source,
              double stretchFactor);

  const ATermBase::CoordinateSystem _coordinateSystem;
  size_t _allocatedWidth, _allocatedHeight;
  std::unique_ptr<class FFTResampler> _resampler;
  bool _downsample;
  WindowFunction::Type _window;
  double _padding;
  bool _overrideFitsPhaseCentre;
  double _overrideRa, _overrideDec;
};

#endif
