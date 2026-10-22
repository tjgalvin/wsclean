#include "tophatconvolution.h"

#include <schaapcommon/math/convolution.h>

#include <mutex>

#include <cmath>

using aocommon::Image;

namespace wsclean::tophat_convolution {

namespace {

void EnsureFftwfPlannerThreadSafe() {
  static std::once_flag once;
  std::call_once(once, []() { schaapcommon::math::MakeFftwfPlannerThreadSafe(); });
}

/// The size of the tophat kernel that Convolve() uses for a given radius.
size_t KernelSize(double radius) {
  size_t kernel_size = static_cast<size_t>(std::floor(radius)) * 2 + 1;
  if (kernel_size % 2 == 0) ++kernel_size;
  return kernel_size;
}

}  // namespace

Image MakeTopHatImage(size_t width, size_t height, double tophat_radius) {
  const ssize_t radius_squared = static_cast<size_t>(tophat_radius * tophat_radius);
  const ssize_t x_mid = (width - 1) / 2;
  const ssize_t y_mid = (height - 1) / 2;
  Image result(width, height, 0.0);
  size_t sum = 0;
  for (ssize_t y = 0; y != static_cast<ssize_t>(height); ++y) {
    float* result_row = &result[y * width];
    const ssize_t dy = y - y_mid;
    for (ssize_t x = 0; x != static_cast<ssize_t>(width); ++x) {
      const ssize_t dx = x - x_mid;
      if (dx * dx + dy * dy <= radius_squared) {
        result_row[x] = 1.0;
        ++sum;
      }
    }
  }
  result *= 1.0f / sum;
  return result;
}

void Convolve(Image& input, double radius) {
  EnsureFftwfPlannerThreadSafe();
  const size_t kernel_size = KernelSize(radius);
  const Image kernel = MakeTopHatImage(kernel_size, kernel_size, radius);
  schaapcommon::math::ResizeAndConvolve(input.Data(), input.Width(), input.Height(), kernel.Data(),
                                        kernel_size);
}

float MaskNoiseThreshold(double radius) {
  const float kernel_size = static_cast<float>(KernelSize(radius));
  return 0.5f / (kernel_size * kernel_size);
}

void RemoveMaskNoise(Image& mask, double radius) {
  const float threshold = MaskNoiseThreshold(radius);
  for (float& value : mask) {
    if (value < threshold) value = 0.0f;
  }
}

}  // namespace wsclean::tophat_convolution
