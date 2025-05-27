#include "fourierdomainrenderer.h"

#include <iostream>
#include <numbers>
#include <complex>

#include <aocommon/image.h>

namespace {

aocommon::UVector<float> MakeTukeyWindow(size_t window_size, float inset_size) {
  aocommon::UVector<float> window(window_size);

  for (size_t i = 0; i < window_size; ++i) {
    const float xi = (0.5f + i) * 2.0f;
    if (xi < window_size - inset_size) {
      const float x = xi / (window_size - inset_size);
      window[i] =
          0.5f * (std::cos((x + 1.0f) * std::numbers::pi_v<float>) + 1.0f);
    } else if (xi < window_size + inset_size) {
      window[i] = 1.0f;
    } else {
      const float x =
          (xi - (window_size + inset_size)) / (window_size - inset_size);
      window[i] = 0.5f * (std::cos(x * std::numbers::pi_v<float>) + 1.0f);
    }
  }

  return window;
}

double EvaluateGaussianCoherence(double u, double v, double major_axis,
                                 double minor_axis, double position_angle) {
  constexpr double pi = std::numbers::pi_v<double>;
  const double phi = 3.0 * pi / 2.0 - position_angle;
  const double fwhm_constant =
      2.0 * std::sqrt(2.0 * std::numbers::ln2_v<double>);

  const double u_prime =
      major_axis / fwhm_constant * (u * std::cos(phi) - v * std::sin(phi));
  const double v_prime =
      minor_axis / fwhm_constant * (u * std::sin(phi) + v * std::cos(phi));

  return std::exp(-2.0 * pi * pi * (u_prime * u_prime + v_prime * v_prime));
}

}  // namespace

namespace wsclean::math {

FourierDomainRenderer::FourierDomainRenderer(size_t subgrid_width,
                                             double pixel_scale)
    : subgrid_width_(subgrid_width + (subgrid_width + 1) % 2),
      pixel_scale_(pixel_scale) {
  subgrid_size_ = subgrid_width_ * subgrid_width_;
  subgrid_.resize(subgrid_size_);

  const float inset = 0.8f * subgrid_width_;
  window_ = MakeTukeyWindow(subgrid_width_, inset);

  fftwf_complex* dummy_fftw_data = reinterpret_cast<fftwf_complex*>(
      fftwf_malloc(subgrid_size_ * sizeof(fftwf_complex)));
  fourier_transform_plan_ =
      fftwf_plan_dft_2d(subgrid_width_, subgrid_width_, dummy_fftw_data,
                        dummy_fftw_data, -1, FFTW_ESTIMATE);
  fftwf_free(dummy_fftw_data);
}

bool FourierDomainRenderer::IsOutlier(int rounded_x, int rounded_y,
                                      int image_width, int image_height) const {
  const int subgrid_half_width = subgrid_width_ / 2;

  const bool left = rounded_x + subgrid_half_width < 0;
  const bool right = rounded_x - subgrid_half_width >= image_width;
  const bool bottom = rounded_y + subgrid_half_width < 0;
  const bool top = rounded_y - subgrid_half_width >= image_height;

  return left || right || bottom || top;
}

void FourierDomainRenderer::RenderModelComponent(
    float* image, const ModelComponent& component, size_t width, size_t height,
    float flux, float x, float y) {
  using namespace std::complex_literals;

  const int rounded_x = std::round(x);
  const int rounded_y = std::round(y);

  const float subpixel_shift_x = x - rounded_x;
  const float subpixel_shift_y = y - rounded_y;

  if (IsOutlier(rounded_x, rounded_y, width, height)) {
    return;
  }

  for (size_t yi = 0; yi < subgrid_width_; ++yi) {
    float y_shift =
        (static_cast<int>(yi) - static_cast<int>(subgrid_width_) / 2) /
        static_cast<float>(subgrid_width_);
    double v = y_shift / pixel_scale_;
    for (size_t xi = 0; xi < subgrid_width_; ++xi) {
      const float x_shift =
          (static_cast<int>(xi) - static_cast<int>(subgrid_width_) / 2) /
          static_cast<float>(subgrid_width_);
      const double u = x_shift / pixel_scale_;

      // Shift input to appropriate quadrant of the image in preparation of the
      // Fourier transform
      const size_t shifted_xi =
          (xi - subgrid_width_ / 2 + subgrid_width_) % subgrid_width_;
      const size_t shifted_yi =
          (yi - subgrid_width_ / 2 + subgrid_width_) % subgrid_width_;
      const size_t shifted_index = shifted_yi * subgrid_width_ + shifted_xi;

      const float source_coherence =
          component.Type() == ModelComponent::GaussianSource
              ? static_cast<float>(EvaluateGaussianCoherence(
                    u, v, component.MajorAxis(), component.MinorAxis(),
                    component.PositionAngle()))
              : 1.0f;
      const float pixel_shift =
          2.0f * std::numbers::pi_v<float> *
          (x_shift * subpixel_shift_x + y_shift * subpixel_shift_y);
      const std::complex<float> subgrid_value =
          window_[yi] * window_[xi] * source_coherence *
          std::complex<float>(std::cos(pixel_shift), std::sin(pixel_shift));
      subgrid_[shifted_index] = subgrid_value;
    }
  }

  fftwf_execute_dft(fourier_transform_plan_,
                    reinterpret_cast<fftwf_complex*>(subgrid_.data()),
                    reinterpret_cast<fftwf_complex*>(subgrid_.data()));

  const float normalization_factor = 1.0f / subgrid_size_;

  const int subgrid_half_width = subgrid_width_ / 2;
  const int x_offset = rounded_x - subgrid_half_width;
  const int y_offset = rounded_y - subgrid_half_width;
  const size_t start_x = std::max<int>(0, x_offset);
  const size_t start_y = std::max<int>(0, y_offset);
  const size_t end_x = std::min<int>(x_offset + subgrid_width_, width);
  const size_t end_y = std::min<int>(y_offset + subgrid_width_, height);
  for (size_t yi = start_y; yi != end_y; ++yi) {
    float* image_ptr = &image[yi * width];
    size_t sub_yi = yi - y_offset;
    for (size_t xi = start_x; xi != end_x; ++xi) {
      size_t sub_xi = xi - x_offset;

      // Shift pixels back to the appropriate quadrant
      size_t shifted_sub_yi =
          (sub_yi - subgrid_half_width + subgrid_width_) % subgrid_width_;
      size_t shifted_sub_xi =
          (sub_xi - subgrid_half_width + subgrid_width_) % subgrid_width_;
      size_t shifted_subgrid_index =
          shifted_sub_yi * subgrid_width_ + shifted_sub_xi;

      float pixel_value =
          flux * normalization_factor * subgrid_[shifted_subgrid_index].real();
      image_ptr[xi] += std::isfinite(pixel_value) ? pixel_value : 0.0f;
    }
  }
}

}  // namespace wsclean::math