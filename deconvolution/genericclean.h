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

  float ExecuteMajorIteration(ImageSet& dirtySet, ImageSet& modelSet,
                              const std::vector<aocommon::Image>& psfs,
                              bool& reachedMajorThreshold) final override;

  virtual std::unique_ptr<DeconvolutionAlgorithm> Clone() const final override {
    return std::unique_ptr<DeconvolutionAlgorithm>(new GenericClean(*this));
  }

 private:
  size_t _convolutionWidth;
  size_t _convolutionHeight;
  const float _convolutionPadding;
  bool _useSubMinorOptimization;

  // Scratch buffer should at least accomodate space for image.Size() floats
  // and is only used to avoid unnecessary memory allocations.
  std::optional<float> findPeak(const aocommon::Image& image,
                                float* scratch_buffer, size_t& x, size_t& y);

  class FFTWManager& _fftwManager;
};

#endif
