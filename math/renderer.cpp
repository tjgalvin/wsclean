
#include "renderer.h"

#include <algorithm>
#include <array>
#include <cmath>

#include <aocommon/imagecoordinates.h>

#include <schaapcommon/math/convolution.h>
#include <schaapcommon/math/restoreimage.h>

using aocommon::ImageCoordinates;

namespace wsclean::renderer {
namespace {
long double Gaussian(long double x, long double sigma) {
  // Evaluate unnormalized Gaussian (unit valued peak, but integral over <-inf,
  // inf> doesn't evaluate to unity)
  const long double xi = x / sigma;
  return std::exp(static_cast<long double>(-0.5) * xi * xi);
}

void RenderPointComponent(aocommon::Image& image,
                          const ImageCoordinateSettings& image_settings,
                          long double position_ra, long double position_dec,
                          long double flux) {
  long double source_l;
  long double source_m;
  ImageCoordinates::RaDecToLM(position_ra, position_dec, image_settings.ra,
                              image_settings.dec, source_l, source_m);
  source_l -= image_settings.l_shift;
  source_m -= image_settings.m_shift;

  int source_x;
  int source_y;
  ImageCoordinates::LMToXY<long double>(
      source_l, source_m, image_settings.pixel_scale_l,
      image_settings.pixel_scale_m, image.Width(), image.Height(), source_x,
      source_y);

  if (source_x >= 0 && source_x < int(image.Width()) && source_y >= 0 &&
      source_y < int(image.Height())) {
    image.Value(source_x, source_y) += static_cast<double>(flux);
  }
}

void RenderGaussianComponent(aocommon::Image& image,
                             const ImageCoordinateSettings& image_settings,
                             long double position_ra, long double position_dec,
                             long double gaus_major_axis,
                             long double gaus_minor_axis,
                             long double gaus_position_angle,
                             long double flux) {
  // Using the FWHM formula for a Gaussian:
  const long double fwhm_constant = (2.0L * std::sqrt(2.0L * std::log(2.0L)));
  const long double sigma_major_axis = gaus_major_axis / fwhm_constant;
  const long double sigma_minor_axis = gaus_minor_axis / fwhm_constant;
  // TODO this won't work for non-equally spaced dimensions
  const long double min_pixel_scale =
      std::min(image_settings.pixel_scale_l, image_settings.pixel_scale_m);

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
  ImageCoordinates::RaDecToLM(position_ra, position_dec, image_settings.ra,
                              image_settings.dec, source_l, source_m);

  // Calculate the bounding box
  int source_x;
  int source_y;
  ImageCoordinates::LMToXY<long double>(
      source_l - image_settings.l_shift, source_m - image_settings.m_shift,
      image_settings.pixel_scale_l, image_settings.pixel_scale_m, image.Width(),
      image.Height(), source_x, source_y);
  const int x_left =
      std::clamp(source_x - bounding_box_size, 0, int(image.Width()));
  const int x_right =
      std::clamp(source_x + bounding_box_size, x_left, int(image.Width()));
  const int y_top =
      std::clamp(source_y - bounding_box_size, 0, int(image.Height()));
  const int y_bottom =
      std::clamp(source_y + bounding_box_size, y_top, int(image.Height()));

  std::vector<double> values;
  double flux_sum = 0.0;
  for (int y = y_top; y != y_bottom; ++y) {
    for (int x = x_left; x != x_right; ++x) {
      long double l, m;
      ImageCoordinates::XYToLM<long double>(
          x, y, image_settings.pixel_scale_l, image_settings.pixel_scale_m,
          image.Width(), image.Height(), l, m);
      l += image_settings.l_shift;
      m += image_settings.m_shift;
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
  // flux_sum can be zero for small or faint sources.
  // Render those as a point source instead.
  if (flux_sum > 0.0) {
    const double* iter = values.data();
    const double factor = flux / flux_sum;
    for (int y = y_top; y != y_bottom; ++y) {
      for (int x = x_left; x != x_right; ++x) {
        image.Value(x, y) += *iter * factor;
        ++iter;
      }
    }
  } else {
    RenderPointComponent(image, image_settings, position_ra, position_dec,
                         flux);
  }
}

void RenderModel(aocommon::Image& image,
                 const ImageCoordinateSettings& image_settings,
                 const Model& model, long double start_frequency,
                 long double end_frequency,
                 aocommon::PolarizationEnum polarization) {
  for (const ModelSource& source : model) {
    for (const ModelComponent& component : source) {
      const long double position_ra = component.PosRA();
      const long double position_dec = component.PosDec();
      const long double integrated_flux = component.SED().IntegratedFlux(
          start_frequency, end_frequency, polarization);

      if (component.Type() == ModelComponent::GaussianSource) {
        const long double gaus_major_axis = component.MajorAxis();
        const long double gaus_minor_axis = component.MinorAxis();
        const long double gaus_position_angle = component.PositionAngle();
        RenderGaussianComponent(image, image_settings, position_ra,
                                position_dec, gaus_major_axis, gaus_minor_axis,
                                gaus_position_angle, integrated_flux);
      } else
        RenderPointComponent(image, image_settings, position_ra, position_dec,
                             integrated_flux);
    }
  }
}

}  // namespace

void RestoreWithEllipticalBeam(aocommon::Image& image,
                               const ImageCoordinateSettings& image_settings,
                               const Model& model, long double beam_major_axis,
                               long double beam_minor_axis,
                               long double beam_position_angle,
                               long double start_frequency,
                               long double end_frequency,
                               aocommon::PolarizationEnum polarization) {
  aocommon::Image rendered_without_beam(image.Width(), image.Height(), 0.0);
  RenderModel(rendered_without_beam, image_settings, model, start_frequency,
              end_frequency, polarization);
  schaapcommon::math::RestoreImage(
      image.Data(), rendered_without_beam.Data(), image.Width(), image.Height(),
      beam_major_axis, beam_minor_axis, beam_position_angle,
      image_settings.pixel_scale_l, image_settings.pixel_scale_m);
}
}  // namespace wsclean::renderer
