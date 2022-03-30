#ifndef PYTHON_DECONVOLUTION_H
#define PYTHON_DECONVOLUTION_H

#include "deconvolutionalgorithm.h"
#include "imageset.h"

#include <aocommon/uvector.h>

#include <schaapcommon/fitters/spectralfitter.h>

#include <pybind11/embed.h>

class PythonDeconvolution : public DeconvolutionAlgorithm {
 public:
  PythonDeconvolution(const std::string& filename);

  float ExecuteMajorIteration(ImageSet& dirtySet, ImageSet& modelSet,
                              const std::vector<aocommon::Image>& psfs,
                              bool& reachedMajorThreshold) final override;

  virtual std::unique_ptr<DeconvolutionAlgorithm> Clone() const final override {
    return std::make_unique<PythonDeconvolution>(*this);
  }

 private:
  std::string _filename;
  // A Python interpreter can not be restarted, so the interpreter
  // needs to live for the entire run
  std::shared_ptr<pybind11::scoped_interpreter> _guard;
  pybind11::function _deconvolveFunction;

  void setBuffer(const ImageSet& imageSet, double* pyPtr);
  void setPsf(const std::vector<aocommon::Image>& psfs, double* pyPtr,
              size_t width, size_t height);
  void getBuffer(ImageSet& imageSet, const double* pyPtr);
};

#endif  // PYTHON_DECONVOLUTION_H
