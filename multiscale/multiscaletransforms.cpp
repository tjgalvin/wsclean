#include "multiscaletransforms.h"

#include "../math/fftconvolver.h"

void MultiScaleTransforms::Transform(std::vector<ImageF>& images,
                                     ImageF& scratch, float scale) {
  size_t kernelSize;
  ImageF shape = MakeShapeFunction(scale, kernelSize);

  scratch = 0.0;

  FFTConvolver::PrepareSmallKernel(scratch.data(), _width, _height,
                                   shape.data(), kernelSize);
  for (ImageF& image : images)
    FFTConvolver::ConvolveSameSize(_fftwManager, image.data(), scratch.data(),
                                   _width, _height);
}

void MultiScaleTransforms::PrepareTransform(float* kernel, float scale) {
  size_t kernelSize;
  ImageF shape = MakeShapeFunction(scale, kernelSize);

  std::fill_n(kernel, _width * _height, 0.0);

  FFTConvolver::PrepareSmallKernel(kernel, _width, _height, shape.data(),
                                   kernelSize);
}

void MultiScaleTransforms::FinishTransform(float* image, const float* kernel) {
  FFTConvolver::ConvolveSameSize(_fftwManager, image, kernel, _width, _height);
}
