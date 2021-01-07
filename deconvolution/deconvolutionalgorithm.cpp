#include "deconvolutionalgorithm.h"

#include "../system/system.h"

DeconvolutionAlgorithm::DeconvolutionAlgorithm()
    : _threshold(0.0),
      _majorIterThreshold(0.0),
      _gain(0.1),
      _mGain(1.0),
      _cleanBorderRatio(0.05),
      _maxIter(500),
      _iterationNumber(0),
      _threadCount(System::ProcessorCount()),
      _allowNegativeComponents(true),
      _stopOnNegativeComponent(false),
      _cleanMask(0),
      _logReceiver(nullptr),
      _spectralFitter(NoSpectralFitting, 0) {}

void DeconvolutionAlgorithm::ResizeImage(float* dest, size_t newWidth,
                                         size_t newHeight, const float* source,
                                         size_t width, size_t height) {
  size_t srcStartX = (width - newWidth) / 2,
         srcStartY = (height - newHeight) / 2;
  for (size_t y = 0; y != newHeight; ++y) {
    float* destPtr = dest + y * newWidth;
    const float* srcPtr = source + (y + srcStartY) * width + srcStartX;
    memcpy(destPtr, srcPtr, newWidth * sizeof(double));
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

void DeconvolutionAlgorithm::PerformSpectralFit(float* values) {
  _spectralFitter.FitAndEvaluate(values);
}
