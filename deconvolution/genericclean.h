#ifndef GENERIC_CLEAN_H
#define GENERIC_CLEAN_H

#include "deconvolutionalgorithm.h"
#include "imageset.h"
#include "simpleclean.h"

#include <aocommon/uvector.h>

#include <optional>

/**
 * This class implements a generalized version of Högbom clean. It performs a
 * single-channel or joined cleaning, depending on the number of images
 * provided. It can use a Clark-like optimization to speed up the cleaning. When
 * multiple frequencies are provided, it can perform spectral fitting.
 */
class GenericClean : public DeconvolutionAlgorithm {
 public:
  explicit GenericClean(class FFTWManager& fftwManager,
                        bool useSubMinorOptimization);

  virtual float ExecuteMajorIteration(
      ImageSet& dirtySet, ImageSet& modelSet,
      const aocommon::UVector<const float*>& psfs, size_t width, size_t height,
      bool& reachedMajorThreshold) final override;

  virtual std::unique_ptr<DeconvolutionAlgorithm> Clone() const final override {
    return std::unique_ptr<DeconvolutionAlgorithm>(new GenericClean(*this));
  }

 private:
  size_t _width, _height, _convolutionWidth, _convolutionHeight;
  float _convolutionPadding;
  bool _useSubMinorOptimization;

  std::optional<float> findPeak(const float* image, float* scratch, size_t& x,
                                size_t& y);

  std::string peakDescription(const float* image, size_t& x, size_t& y);

  void subtractImage(float* image, const float* psf, size_t x, size_t y,
                     float factor, size_t startY, size_t endY) const {
    SimpleClean::PartialSubtractImage(image, _width, _height, psf, _width,
                                      _height, x, y, factor, startY, endY);
  }

  class FFTWManager& _fftwManager;
};

#endif
