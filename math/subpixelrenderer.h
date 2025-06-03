#ifndef WSCLEAN_INTERPOLATING_RENDERER_H_
#define WSCLEAN_INTERPOLATING_RENDERER_H_

#include <cstring>

#include <aocommon/coordinatesystem.h>
#include <aocommon/image.h>
#include <aocommon/uvector.h>

namespace wsclean::math {

class SubPixelRenderer {
 public:
  SubPixelRenderer(size_t kernel_size)
      : x_kernel_(kernel_size + (kernel_size + 1) % 2),
        y_kernel_(kernel_size + (kernel_size + 1) % 2) {}

  /**
   * Render a source and convolve it with a sinc, which means it can be on
   * non-integer positions.
   */
  static void RenderSource(float* image, size_t width, size_t height,
                           float flux, double x, double y);

  /**
   * Render a source and convolve it with a sinc, which means it can be on
   * non-integer positions. The sinc is windowed to increase the performance. A
   * Hann-window is used. No corrections are made to correct for the window,
   * which implies it needs to be sufficient big to not cause errors.
   */
  void RenderWindowedSource(float* image, size_t width, size_t height,
                            float brightness, float x, float y);

 private:
  size_t kernel_size_;
  aocommon::UVector<float> x_kernel_;
  aocommon::UVector<float> y_kernel_;
};

std::vector<aocommon::Image> RenderSubPixelModel(
    const std::string& model_filename,
    const aocommon::CoordinateSystem& coordinate_system, double frequency,
    double bandwidth, size_t window_size, size_t n_terms, double mem_fraction,
    double mem_limit);

}  // namespace wsclean::math

#endif
