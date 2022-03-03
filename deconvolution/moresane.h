#ifndef MORESANE_H
#define MORESANE_H

#include <string>

#include "deconvolutionalgorithm.h"
#include "imageset.h"

class MoreSane : public DeconvolutionAlgorithm {
 public:
  MoreSane(const std::string& moreSaneLocation,
           const std::string& moresaneArguments,
           const std::vector<double>& moresaneSigmaLevels,
           const std::string& prefixName, class FFTWManager& fftwManager)
      : _moresaneLocation(moreSaneLocation),
        _moresaneArguments(moresaneArguments),
        _moresaneSigmaLevels(moresaneSigmaLevels),
        _prefixName(prefixName),
        _fftwManager(fftwManager) {}

  float ExecuteMajorIteration(ImageSet& dataImage, ImageSet& modelImage,
                              const std::vector<aocommon::Image>& psfImages,
                              bool& reachedMajorThreshold) final override;

  virtual std::unique_ptr<DeconvolutionAlgorithm> Clone() const final override {
    return std::unique_ptr<DeconvolutionAlgorithm>(new MoreSane(*this));
  }

  void ExecuteMajorIteration(float* residualData, float* modelData,
                             const aocommon::Image& psfImage);

 private:
  const std::string _moresaneLocation, _moresaneArguments;

  const std::vector<double> _moresaneSigmaLevels;
  const std::string _prefixName;

  class FFTWManager& _fftwManager;
};

#endif
