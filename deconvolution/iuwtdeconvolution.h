#ifndef IUWT_DECONVOLUTION_H
#define IUWT_DECONVOLUTION_H

#include "deconvolutionalgorithm.h"
#include "imageset.h"

#include "../iuwt/iuwtdeconvolutionalgorithm.h"

#include "../structures/imagingtable.h"

#include <aocommon/uvector.h>

#include <memory>
#include <string>

class IUWTDeconvolution : public DeconvolutionAlgorithm {
 public:
  IUWTDeconvolution() : _useSNRTest(false) {}

  float ExecuteMajorIteration(ImageSet& dataImage, ImageSet& modelImage,
                              const std::vector<aocommon::Image>& psfImages,
                              bool& reachedMajorThreshold) final override {
    IUWTDeconvolutionAlgorithm algorithm(
        dataImage.Width(), dataImage.Height(), _gain, _mGain, _cleanBorderRatio,
        _allowNegativeComponents, _cleanMask, _threshold, _useSNRTest);
    float val = algorithm.PerformMajorIteration(
        _iterationNumber, MaxNIter(), modelImage, dataImage, psfImages,
        reachedMajorThreshold);
    if (_iterationNumber >= MaxNIter()) reachedMajorThreshold = false;
    return val;
  }

  std::unique_ptr<DeconvolutionAlgorithm> Clone() const final override {
    return std::make_unique<IUWTDeconvolution>(*this);
  }

  void SetUseSNRTest(bool useSNRTest) { _useSNRTest = useSNRTest; }

 private:
  bool _useSNRTest;
};

#endif
