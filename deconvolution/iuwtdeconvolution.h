#ifndef IUWT_DECONVOLUTION_H
#define IUWT_DECONVOLUTION_H

#include "deconvolutionalgorithm.h"
#include "imageset.h"

#include "../iuwt/iuwtdeconvolutionalgorithm.h"

#include "../wsclean/imagingtable.h"

#include <aocommon/uvector.h>

#include <memory>
#include <string>

class IUWTDeconvolution : public DeconvolutionAlgorithm {
 public:
  IUWTDeconvolution(class FFTWManager& fftwManager)
      : _fftwManager(fftwManager), _useSNRTest(false) {}

  virtual double ExecuteMajorIteration(
      ImageSet& dataImage, ImageSet& modelImage,
      const aocommon::UVector<const double*>& psfImages, size_t width,
      size_t height, bool& reachedMajorThreshold) final override {
    IUWTDeconvolutionAlgorithm algorithm(
        _fftwManager, width, height, _gain, _mGain, _cleanBorderRatio,
        _allowNegativeComponents, _cleanMask, _threshold, _useSNRTest);
    double val = algorithm.PerformMajorIteration(
        _iterationNumber, MaxNIter(), modelImage, dataImage, psfImages,
        reachedMajorThreshold);
    if (_iterationNumber >= MaxNIter()) reachedMajorThreshold = false;
    return val;
  }

  std::unique_ptr<DeconvolutionAlgorithm> Clone() const final override {
    return std::unique_ptr<IUWTDeconvolution>(new IUWTDeconvolution(*this));
  }

  void SetUseSNRTest(bool useSNRTest) { _useSNRTest = useSNRTest; }

 private:
  class FFTWManager& _fftwManager;
  bool _useSNRTest;
};

#endif
