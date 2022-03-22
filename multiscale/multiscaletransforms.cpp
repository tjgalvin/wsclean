#include "multiscaletransforms.h"

#include <schaapcommon/fft/convolution.h>

using aocommon::Image;

void MultiScaleTransforms::Transform(std::vector<Image>& images, Image& scratch,
                                     float scale) {
  size_t kernelSize;
  Image shape = MakeShapeFunction(scale, kernelSize);

  scratch = 0.0;

  schaapcommon::fft::PrepareSmallConvolutionKernel(
      scratch.Data(), _width, _height, shape.Data(), kernelSize, _threadCount);
  for (Image& image : images)
    schaapcommon::fft::Convolve(image.Data(), scratch.Data(), _width, _height,
                                _threadCount);
}

void MultiScaleTransforms::PrepareTransform(float* kernel, float scale) {
  size_t kernelSize;
  Image shape = MakeShapeFunction(scale, kernelSize);

  std::fill_n(kernel, _width * _height, 0.0);

  schaapcommon::fft::PrepareSmallConvolutionKernel(
      kernel, _width, _height, shape.Data(), kernelSize, _threadCount);
}

void MultiScaleTransforms::FinishTransform(float* image, const float* kernel) {
  schaapcommon::fft::Convolve(image, kernel, _width, _height, _threadCount);
}
