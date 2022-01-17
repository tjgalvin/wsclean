#include "deconvolutionalgorithm.h"

#include <aocommon/system.h>

DeconvolutionAlgorithm::DeconvolutionAlgorithm()
    : _threshold(0.0),
      _majorIterThreshold(0.0),
      _gain(0.1),
      _mGain(1.0),
      _cleanBorderRatio(0.05),
      _maxIter(500),
      _iterationNumber(0),
      _threadCount(aocommon::system::ProcessorCount()),
      _allowNegativeComponents(true),
      _stopOnNegativeComponent(false),
      _cleanMask(nullptr),
      _logReceiver(nullptr),
      _spectralFitter(SpectralFittingMode::NoFitting, 0) {}

void DeconvolutionAlgorithm::ResizeImage(float* dest, size_t newWidth,
                                         size_t newHeight, const float* source,
                                         size_t width, size_t height) {
  size_t srcStartX = (width - newWidth) / 2,
         srcStartY = (height - newHeight) / 2;
  for (size_t y = 0; y != newHeight; ++y) {
    float* destPtr = dest + y * newWidth;
    const float* srcPtr = source + (y + srcStartY) * width + srcStartX;
    std::copy_n(srcPtr, newWidth, destPtr);
  }
}

void DeconvolutionAlgorithm::RemoveNaNsInPSF(float* psf, size_t width,
                                             size_t height) {
  float* endPtr = psf + width * height;
  while (psf != endPtr) {
    if (!std::isfinite(*psf)) *psf = 0.0;
    ++psf;
  }
}

void DeconvolutionAlgorithm::PerformSpectralFit(float* values, size_t x,
                                                size_t y) const {
  _spectralFitter.FitAndEvaluate(values, x, y, _fittingScratch);
}
