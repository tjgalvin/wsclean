#ifndef RENDERER_H_
#define RENDERER_H_

#include <cstring>

#include <aocommon/fits/fitsreader.h>
#include <aocommon/image.h>
#include <aocommon/polarization.h>

#include "../model/model.h"

namespace renderer {

/**
 * @brief Struct collecting the relevant image coordinate settings for
 * rendereing the source components.
 *
 */
struct ImageCoordinateSettings {
  ImageCoordinateSettings() = default;

  /**
   * @brief Extract coordinate settings from an aocommon::FitsReader object.
   *
   * @param fits_reader aocommon::FitsReader object.
   */
  ImageCoordinateSettings(const aocommon::FitsReader& fits_reader)
      : ra(fits_reader.PhaseCentreRA()),
        dec(fits_reader.PhaseCentreDec()),
        pixel_scale_l(fits_reader.PixelSizeX()),
        pixel_scale_m(fits_reader.PixelSizeY()),
        l_shift(fits_reader.PhaseCentreDL()),
        m_shift(fits_reader.PhaseCentreDM()) {}

  long double ra;
  long double dec;
  long double pixel_scale_l;
  long double pixel_scale_m;
  long double l_shift;
  long double m_shift;
};

/**
 * @brief Restore a model image by convolving it with an elliptical Gaussian.
 *
 * @param image Image to which restored sources are written.
 * @param image_settings Image coordinate settings.
 * @param model Modeled source components.
 * @param beam_major_axis Major axis of elliptical beam to be applied [rad].
 * @param beam_minor_axis Minor axis of elliptical beam to be applied [rad].
 * @param beam_position_angle Position angle of beam [rad].
 * @param start_frequency Start frequency [Hz].
 * @param end_frequency End frequency [Hz].
 * @param polarization Polarization enum.
 * @param thread_count Numbers of threads to use.
 */
void RestoreWithEllipticalBeam(
    aocommon::Image& image, const ImageCoordinateSettings& image_settings,
    const Model& model, long double beam_major_axis,
    long double beam_minor_axis, long double beam_position_angle,
    long double start_frequency, long double end_frequency,
    aocommon::PolarizationEnum polarization, size_t thread_count);

}  // namespace renderer

#endif
