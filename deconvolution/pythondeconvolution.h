#ifndef PYTHON_DECONVOLUTION_H
#define PYTHON_DECONVOLUTION_H

#include "deconvolutionalgorithm.h"
#include "imageset.h"

#include <aocommon/uvector.h>

#include <pybind11/embed.h>

class PythonDeconvolution : public DeconvolutionAlgorithm {
 public:
  PythonDeconvolution(const std::string& filename);

  virtual float ExecuteMajorIteration(
      ImageSet& dirtySet, ImageSet& modelSet,
      const aocommon::UVector<const float*>& psfs, size_t width, size_t height,
      bool& reachedMajorThreshold) final override;

  virtual std::unique_ptr<DeconvolutionAlgorithm> Clone() const final override {
    return std::unique_ptr<DeconvolutionAlgorithm>(
        new PythonDeconvolution(*this));
  }

 private:
  std::string _filename;
  // A Python interpreter can not be restarted, so the interpreter
  // needs to live for the entire run
  std::shared_ptr<pybind11::scoped_interpreter> _guard;
  pybind11::function _deconvolveFunction;

  void setBuffer(const ImageSet& imageSet, double* pyPtr, size_t width,
                 size_t height);
  void setPsf(const aocommon::UVector<const float*>& psfs, double* pyPtr,
              size_t width, size_t height);
  void getBuffer(ImageSet& imageSet, const double* pyPtr, size_t width,
                 size_t height);
};

#endif  // PYTHON_DECONVOLUTION_H
