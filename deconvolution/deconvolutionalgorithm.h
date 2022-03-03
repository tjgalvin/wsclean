#ifndef CLEAN_ALGORITHM_H
#define CLEAN_ALGORITHM_H

#include <string>
#include <cmath>

#include "spectralfitter.h"

#include <aocommon/image.h>
#include <aocommon/logger.h>
#include <aocommon/polarization.h>
#include <aocommon/uvector.h>

namespace ao {
template <typename T>
class lane;
}

class DeconvolutionAlgorithm {
 public:
  virtual ~DeconvolutionAlgorithm() {}

  virtual float ExecuteMajorIteration(
      class ImageSet& dataImage, class ImageSet& modelImage,
      const std::vector<aocommon::Image>& psfImages,
      bool& reachedMajorThreshold) = 0;

  virtual std::unique_ptr<DeconvolutionAlgorithm> Clone() const = 0;

  void SetMaxNIter(size_t nIter) { _maxIter = nIter; }

  void SetThreshold(float threshold) { _threshold = threshold; }

  void SetMajorIterThreshold(float mThreshold) {
    _majorIterThreshold = mThreshold;
  }

  void SetGain(float gain) { _gain = gain; }

  void SetMGain(float mGain) { _mGain = mGain; }

  void SetAllowNegativeComponents(bool allowNegativeComponents) {
    _allowNegativeComponents = allowNegativeComponents;
  }

  void SetStopOnNegativeComponents(bool stopOnNegative) {
    _stopOnNegativeComponent = stopOnNegative;
  }

  void SetCleanBorderRatio(float borderRatio) {
    _cleanBorderRatio = borderRatio;
  }

  void SetThreadCount(size_t threadCount) { _threadCount = threadCount; }

  void SetLogReceiver(aocommon::LogReceiver& receiver) {
    _logReceiver = &receiver;
  }

  size_t MaxNIter() const { return _maxIter; }
  float Threshold() const { return _threshold; }
  float MajorIterThreshold() const { return _majorIterThreshold; }
  float Gain() const { return _gain; }
  float MGain() const { return _mGain; }
  float CleanBorderRatio() const { return _cleanBorderRatio; }
  bool AllowNegativeComponents() const { return _allowNegativeComponents; }
  bool StopOnNegativeComponents() const { return _stopOnNegativeComponent; }

  void SetCleanMask(const bool* cleanMask) { _cleanMask = cleanMask; }

  size_t IterationNumber() const { return _iterationNumber; }

  void SetIterationNumber(size_t iterationNumber) {
    _iterationNumber = iterationNumber;
  }

  static void ResizeImage(float* dest, size_t newWidth, size_t newHeight,
                          const float* source, size_t width, size_t height);

  static void RemoveNaNsInPSF(float* psf, size_t width, size_t height);

  void CopyConfigFrom(const DeconvolutionAlgorithm& source) {
    _threshold = source._threshold;
    _gain = source._gain;
    _mGain = source._mGain;
    _cleanBorderRatio = source._cleanBorderRatio;
    _maxIter = source._maxIter;
    // skip _iterationNumber
    _allowNegativeComponents = source._allowNegativeComponents;
    _stopOnNegativeComponent = source._stopOnNegativeComponent;
    _cleanMask = source._cleanMask;
    _spectralFitter = source._spectralFitter;
  }

  void SetSpectralFittingMode(SpectralFittingMode mode, size_t nTerms) {
    _spectralFitter.SetMode(mode, nTerms);
  }

  void SetSpectrallyForcedImages(std::vector<aocommon::Image>&& images) {
    _spectralFitter.SetForcedImages(std::move(images));
  }

  void InitializeFrequencies(const aocommon::UVector<double>& frequencies,
                             const aocommon::UVector<float>& weights) {
    _spectralFitter.SetFrequencies(frequencies.data(), weights.data(),
                                   frequencies.size());
  }

  const SpectralFitter& Fitter() const { return _spectralFitter; }

  void SetRMSFactorImage(aocommon::Image&& image) {
    _rmsFactorImage = std::move(image);
  }
  const aocommon::Image& RMSFactorImage() const { return _rmsFactorImage; }

 protected:
  DeconvolutionAlgorithm();

  DeconvolutionAlgorithm(const DeconvolutionAlgorithm&) = default;
  DeconvolutionAlgorithm& operator=(const DeconvolutionAlgorithm&) = default;

  void PerformSpectralFit(float* values, size_t x, size_t y) const;

  float _threshold, _majorIterThreshold, _gain, _mGain, _cleanBorderRatio;
  size_t _maxIter, _iterationNumber, _threadCount;
  bool _allowNegativeComponents, _stopOnNegativeComponent;
  const bool* _cleanMask;
  aocommon::Image _rmsFactorImage;
  mutable aocommon::UVector<float> _fittingScratch;

  aocommon::LogReceiver* _logReceiver;

  SpectralFitter _spectralFitter;
};

#endif
