
#include "renderer.h"

#include <aocommon/imagecoordinates.h>

#include <schaapcommon/fft/convolution.h>
#include <schaapcommon/fft/restoreimage.h>

#include <cmath>
#include <array>

#include <boost/algorithm/clamp.hpp>

using aocommon::ImageCoordinates;
using boost::algorithm::clamp;

namespace renderer {
namespace {
long double Gaussian(long double x, long double sigma) {
  // Evaluate unnormalized Gaussian (unit valued peak, but integral over <-inf,
  // inf> doesn't evaluate to unity)
  const long double xi = x / sigma;
  return std::exp(static_cast<long double>(-0.5) * xi * xi);
}

void RenderPointComponent(float* image_data, size_t image_width,
                          size_t image_height, long double ra, long double dec,
                          long double pixel_scale_l, long double pixel_scale_m,
                          long double l_shift, long double m_shift,
                          long double position_ra, long double position_dec,
                          long double flux) {
  long double source_l;
  long double source_m;
  ImageCoordinates::RaDecToLM(position_ra, position_dec, ra, dec, source_l,
                              source_m);
  source_l -= l_shift;
  source_m -= m_shift;

  int source_x;
  int source_y;
  ImageCoordinates::LMToXY<long double>(source_l, source_m, pixel_scale_l,
                                        pixel_scale_m, image_width,
                                        image_height, source_x, source_y);

  if (source_x >= 0 && source_x < (int)image_width && source_y >= 0 &&
      source_y < (int)image_height) {
    float* image_data_ptr = image_data + source_y * image_width + source_x;
    (*image_data_ptr) += static_cast<double>(flux);
  }
}

void RenderGaussianComponent(
    float* image_data, size_t image_width, size_t image_height, long double ra,
    long double dec, long double pixel_scale_l, long double pixel_scale_m,
    long double l_shift, long double m_shift, long double position_ra,
    long double position_dec, long double gaus_major_axis,
    long double gaus_minor_axis, long double gaus_position_angle,
    long double flux) {
  // Using the FWHM formula for a Gaussian:
  const long double fwhm_constant = (2.0L * std::sqrt(2.0L * std::log(2.0L)));
  const long double sigma_major_axis = gaus_major_axis / fwhm_constant;
  const long double sigma_minor_axis = gaus_minor_axis / fwhm_constant;
  // TODO this won't work for non-equally spaced dimensions
  const long double min_pixel_scale = std::min(pixel_scale_l, pixel_scale_m);

  // Position angle is angle from North:
  const long double angle = gaus_position_angle + M_PI_2;
  const long double cos_angle = std::cos(angle);
  const long double sin_angle = std::sin(angle);

  // Make rotation matrix
  std::array<long double, 4> transf;
  transf[0] = cos_angle;
  transf[1] = -sin_angle;
  transf[2] = sin_angle;
  transf[3] = cos_angle;

  const double sigma_max = std::max(std::fabs(sigma_major_axis * transf[0]),
                                    std::fabs(sigma_major_axis * transf[1]));
  // Multiply with scaling matrix to make variance 1.
  transf[0] /= sigma_major_axis;
  transf[1] /= sigma_major_axis;
  transf[2] /= sigma_minor_axis;
  transf[3] /= sigma_minor_axis;
  const int bounding_box_size = std::ceil(sigma_max * 20.0 / min_pixel_scale);
  long double source_l;
  long double source_m;
  ImageCoordinates::RaDecToLM(position_ra, position_dec, ra, dec, source_l,
                              source_m);

  // Calculate the bounding box
  int source_x;
  int source_y;
  ImageCoordinates::LMToXY<long double>(
      source_l - l_shift, source_m - m_shift, pixel_scale_l, pixel_scale_m,
      image_width, image_height, source_x, source_y);
  const int x_left = clamp(source_x - bounding_box_size, 0, int(image_width));
  const int x_right =
      clamp(source_x + bounding_box_size, x_left, int(image_width));
  const int y_top = clamp(source_y - bounding_box_size, 0, int(image_height));
  const int y_bottom =
      clamp(source_y + bounding_box_size, y_top, int(image_height));

  std::vector<double> values;
  double flux_sum = 0.0;
  for (int y = y_top; y != y_bottom; ++y) {
    for (int x = x_left; x != x_right; ++x) {
      long double l, m;
      ImageCoordinates::XYToLM<long double>(x, y, pixel_scale_l, pixel_scale_m,
                                            image_width, image_height, l, m);
      l += l_shift;
      m += m_shift;
      const long double l_transf =
          (l - source_l) * transf[0] + (m - source_m) * transf[1];
      const long double m_transf =
          (l - source_l) * transf[2] + (m - source_m) * transf[3];
      const long double dist =
          std::sqrt(l_transf * l_transf + m_transf * m_transf);
      long double v = Gaussian(dist, 1.0L) * flux;
      flux_sum += static_cast<double>(v);
      values.emplace_back(v);
    }
  }
  const double* iter = values.data();
  const double factor = flux / flux_sum;
  for (int y = y_top; y != y_bottom; ++y) {
    float* image_data_ptr = image_data + y * image_width + x_left;
    for (int x = x_left; x != x_right; ++x) {
      (*image_data_ptr) += *iter * factor;
      ++image_data_ptr;
      ++iter;
    }
  }
}

void RenderModel(aocommon::Image& image,
                 const ImageCoordinateSettings& image_settings,
                 const Model& model, long double start_frequency,
                 long double end_frequency,
                 aocommon::PolarizationEnum polarization) {
  for (const ModelSource& source : model) {
    for (const ModelComponent& component : source) {
      const long double position_ra = component.PosRA(),
                        position_dec = component.PosDec();
      const long double integrated_flux = component.SED().IntegratedFlux(
          start_frequency, end_frequency, polarization);

      if (component.Type() == ModelComponent::GaussianSource) {
        const long double gaus_major_axis = component.MajorAxis(),
                          gaus_minor_axis = component.MinorAxis();
        const long double gaus_position_angle = component.PositionAngle();
        RenderGaussianComponent(
            image.Data(), image.Width(), image.Height(), image_settings.ra,
            image_settings.dec, image_settings.pixel_scale_l,
            image_settings.pixel_scale_m, image_settings.l_shift,
            image_settings.m_shift, position_ra, position_dec, gaus_major_axis,
            gaus_minor_axis, gaus_position_angle, integrated_flux);
      } else
        RenderPointComponent(
            image.Data(), image.Width(), image.Height(), image_settings.ra,
            image_settings.dec, image_settings.pixel_scale_l,
            image_settings.pixel_scale_m, image_settings.l_shift,
            image_settings.m_shift, position_ra, position_dec, integrated_flux);
    }
  }
}

}  // namespace

void RestoreWithEllipticalBeam(
    aocommon::Image& image, const ImageCoordinateSettings& image_settings,
    const Model& model, long double beam_major_axis,
    long double beam_minor_axis, long double beam_position_angle,
    long double start_frequency, long double end_frequency,
    aocommon::PolarizationEnum polarization, size_t thread_count) {
  aocommon::Image rendered_without_beam(image.Width(), image.Height(), 0.0);
  RenderModel(rendered_without_beam, image_settings, model, start_frequency,
              end_frequency, polarization);
  schaapcommon::fft::RestoreImage(
      image.Data(), rendered_without_beam.Data(), image.Width(), image.Height(),
      beam_major_axis, beam_minor_axis, beam_position_angle,
      image_settings.pixel_scale_l, image_settings.pixel_scale_m, thread_count);
}
}  // namespace renderer