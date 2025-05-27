#ifndef WSCLEAN_FOURIER_DOMAIN_RENDERER_H_
#define WSCLEAN_FOURIER_DOMAIN_RENDERER_H_

#include <aocommon/uvector.h>

#include <fftw3.h>

#include "../model/modelcomponent.h"

namespace wsclean::math {

/**
 * Class that can render a list of sources by sampling sky model in
 * the Fourier domain. A sub-pixel shift is aplied in Fourier space, hence,
 * sources can be rendered on non-integer pixel positions.
 */
class FourierDomainRenderer {
 public:
  FourierDomainRenderer(size_t subgrid_width, double pixel_scale);

  FourierDomainRenderer() = delete;

  ~FourierDomainRenderer() { fftwf_destroy_plan(fourier_transform_plan_); }

  /**
   * @brief Samples a source in the Fourier domain, applies a sub-pixel shift
   * and renders the source onto an image.
   * @details Beware, this function is not thread safe.
   */
  void RenderModelComponent(float* image, const ModelComponent& component,
                            size_t width, size_t height, float flux, float x,
                            float y);

 private:
  /**
   * @brief Determines whether a source is an outlier. This is the case when the
   * sampled sub-grid is outside of the boundary of the image.
   */
  bool IsOutlier(int rounded_x, int rounded_y, int image_width,
                 int image_height) const;

  size_t subgrid_width_;
  size_t subgrid_size_;
  double pixel_scale_;
  aocommon::UVector<float> window_;
  std::vector<std::complex<float>> subgrid_;
  fftwf_plan fourier_transform_plan_;
};

}  // namespace wsclean::math

#endif
