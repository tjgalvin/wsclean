#ifndef WSCLEAN_WTOWERS_GRIDDER_IMPL_H_
#define WSCLEAN_WTOWERS_GRIDDER_IMPL_H_

#include "wtowers_gridder.h"

#include <complex>
#include <cstddef>
#include <vector>

#include <aocommon/image.h>
#include <aocommon/logger.h>
#include <ska-sdp-func/utility/sdp_errors.h>
#include <ska-sdp-func/grid_data/sdp_grid_wstack_wtower.h>
#include <ska-sdp-func/grid_data/sdp_gridder_utils.h>
#include <ska-sdp-func/grid_data/sdp_gridder_wtower_height.h>
#include <ska-sdp-func/fourier_transforms/sdp_fft_padded_size.h>
#include <ska-sdp-func/utility/sdp_mem.h>

#include "../gridding/msgridder.h"

using aocommon::Logger;

/*
 * This file contains the implementation of various template methods from @ref
 * WTowersGridder
 */

namespace wsclean {

template <typename NumT>
WTowersGridder<NumT>::WTowersGridder(size_t width, size_t height,
                                     size_t trimmed_width,
                                     size_t trimmed_height, double pixel_size_x,
                                     double pixel_size_y, double l_shift,
                                     double m_shift, size_t n_threads,
                                     double accuracy, size_t verbosity)
    : width_(width),
      height_(height),
      trimmed_width_(trimmed_width),
      trimmed_height_(trimmed_height),
      n_threads_(n_threads),
      pixel_size_x_(pixel_size_x),
      pixel_size_y_(pixel_size_y),
      l_shift_(l_shift),
      m_shift_(m_shift),
      verbosity_(verbosity) {
  assert(verbosity <= 2);

  if (l_shift_ != 0.0f) {
    throw std::runtime_error("w-towers does not yet support l shift");
  }
  if (m_shift_ != 0.0f) {
    throw std::runtime_error("w-towers does not yet support m shift");
  }
  /// Various configuration paramaters for W-towers, we calculate a few of them
  /// and hardcode the rest Future implementations would make these user
  /// configurable or do a better job of au to calculation.

  wtowers_parameters_.cell_size_rad = pixel_size_x_;
  wtowers_parameters_.image_size = trimmed_width_;
  wtowers_parameters_.verbosity = verbosity > 0 ? 1 : 0;

  // W-towers code uses "cell_size_rad"; in WSClean code this is the equivalent
  // of:
  //   - pixel_size_x_; The angular width of a pixel in radians.
  //   - pixel_size_y_; The angular height of a pixel in radians.
  // W-towers does not support rectangular images so these values must match.
  if (pixel_size_x_ != pixel_size_y_) {
    throw std::runtime_error("w-towers does not yet support non-square images");
  }

  // W-towers uses "image_size"; in WSClean code that is the equivalent of:
  //   - width; The width of the untrimmed image in pixels
  //   - height; The height of the untrimmed image in pixels.
  // Or potentially:
  //   - width_t; The width of the trimmed image in pixels
  //   - trimmed_height_; The height of the trimmed image in pixels.
  // W-towers does not support rectangular images so these values must match.
  if (width_ != height_) {
    throw std::runtime_error("w-towers does not yet support non-square images");
  }

  wtowers_parameters_.field_of_view =
      std::sin(wtowers_parameters_.cell_size_rad) *
      wtowers_parameters_.image_size;
  wtowers_parameters_.grid_size =
      2 * (sdp_fft_padded_size(int(wtowers_parameters_.image_size * 0.5),
                               wtowers_parameters_.padding_factor));
  wtowers_parameters_.grid_resolution =
      std::sin(wtowers_parameters_.cell_size_rad) *
      wtowers_parameters_.grid_size;
  wtowers_parameters_.w_step = sdp_gridder_determine_w_step(
      wtowers_parameters_.grid_resolution, wtowers_parameters_.field_of_view,
      0.0, 0.0, 0.0);

  sdp_Error status = SDP_SUCCESS;
  // NB! We intentionally use 2 * subgrid_size instead of image_size
  // as despite the parameter name this is what the documentation suggests.
  wtowers_parameters_.w_towers_height =
      sdp_gridder_determine_max_w_tower_height(
          2 * wtowers_parameters_.subgrid_size,
          wtowers_parameters_.subgrid_size, wtowers_parameters_.grid_resolution,
          wtowers_parameters_.w_step, wtowers_parameters_.shear_u,
          wtowers_parameters_.shear_v, wtowers_parameters_.support,
          wtowers_parameters_.oversampling, wtowers_parameters_.w_support,
          wtowers_parameters_.w_oversampling, wtowers_parameters_.field_of_view,
          wtowers_parameters_.subgrid_frac, 3, accuracy, &status);

  if (wtowers_parameters_.w_towers_height == 0.0f) {
    Logger::Info << "WARNING: Unable to determine a suitable w_towers_height. "
                    "Using default of 150.\n";
    wtowers_parameters_.w_towers_height = 150;
  }
  if (status != SDP_SUCCESS) {
    throw std::runtime_error("Error computing w_towers_height");
  }
}

template <typename NumT>
size_t WTowersGridder<NumT>::ConstantMemoryUsage() const {
  // Storage for "grid": pessimistically assume an oversampling factor of 2
  size_t constant =
      4.0 * trimmed_width_ * trimmed_height_ * sizeof(std::complex<float>);
  // For prediction, we also need a copy of the dirty image
  constant +=
      trimmed_width_ * trimmed_height_ * sizeof(NumT);  // trimmed dirty image
  return constant;
}

template <typename NumT>
size_t WTowersGridder<NumT>::PerVisibilityMemoryUsage() const {
  // For now we assume this is the same as wgridder.
  // See comments in wgridder/wgridder_implementation.h for how this size
  // is picked. This and ConstantMemoryUsage() should be reworked to more
  // precise w-towers specific estimates.
  return 8;
}

template <typename NumT>
void WTowersGridder<NumT>::InitializeInversion() {
  image_.assign(trimmed_width_ * trimmed_height_, 0);
}

template <typename NumT>
void WTowersGridder<NumT>::LogParameters() const {
  Logger::Debug << "field_of_view: " << wtowers_parameters_.field_of_view
                << "\n";
  Logger::Debug << "grid_resolution / theta: "
                << wtowers_parameters_.grid_resolution << "\n";
  Logger::Debug << "image_size: " << wtowers_parameters_.image_size << "\n";
  Logger::Debug << "grid_size: " << wtowers_parameters_.grid_size << "\n";
  Logger::Debug << "subgrid_size: " << wtowers_parameters_.subgrid_size << "\n";
  Logger::Debug << "cell_size_rad: " << wtowers_parameters_.cell_size_rad
                << "\n";
  Logger::Debug << "w_step: " << wtowers_parameters_.w_step << "\n";
  Logger::Debug << "shear_u: " << wtowers_parameters_.shear_u << "\n";
  Logger::Debug << "shear_v: " << wtowers_parameters_.shear_v << "\n";
  Logger::Debug << "support: " << wtowers_parameters_.support << "\n";
  Logger::Debug << "oversampling: " << wtowers_parameters_.oversampling << "\n";
  Logger::Debug << "w_support: " << wtowers_parameters_.w_support << "\n";
  Logger::Debug << "w_oversampling: " << wtowers_parameters_.w_oversampling
                << "\n";
  Logger::Debug << "subgrid_frac: " << wtowers_parameters_.subgrid_frac << "\n";
  Logger::Debug << "w_towers_height: " << wtowers_parameters_.w_towers_height
                << "\n";
}

template <typename NumT>
void WTowersGridder<NumT>::AddInversionData(
    size_t n_rows, size_t n_channels, const double *uvws,
    const double *frequencies, const std::complex<float> *visibilities) {
  const bool decreasing_freq =
      (n_channels > 1) && (frequencies[1] < frequencies[0]);
  if (decreasing_freq) {
    throw std::runtime_error(
        "W-towers does not currently support frequencies that aren't in "
        "ascending order\n");
  }

  const double frequency_step = frequencies[1] - frequencies[0];

  sdp_MemType image_data_type = SDP_MEM_FLOAT;
  if constexpr (std::is_same_v<NumT, double>) {
    image_data_type = SDP_MEM_DOUBLE;
  }

  sdp_Error status = SDP_SUCCESS;

  std::vector<NumT> dirty_image;
  dirty_image.assign(
      wtowers_parameters_.grid_size * wtowers_parameters_.grid_size, 0);

  const std::vector<int64_t> visibilities_shape{
      static_cast<int64_t>(n_rows), static_cast<int64_t>(n_channels)};
  sdp_Mem *wrapped_visibilities = sdp_mem_create_wrapper(
      (void *)visibilities, SDP_MEM_COMPLEX_FLOAT, SDP_MEM_CPU, 2,
      visibilities_shape.data(), nullptr, &status);
  if (status != SDP_SUCCESS) {
    throw std::runtime_error(
        "w-towers: Fatal error while wrapping memory for visibilities");
  }
  const std::vector<int64_t> uvws_shape{static_cast<int64_t>(n_rows), 3};
  sdp_Mem *wrapped_uvws =
      sdp_mem_create_wrapper((void *)uvws, SDP_MEM_DOUBLE, SDP_MEM_CPU, 2,
                             uvws_shape.data(), nullptr, &status);
  if (status != SDP_SUCCESS) {
    throw std::runtime_error(
        "w-towers: Fatal error while wrapping memory for uvws");
  }
  const std::vector<int64_t> image_shape{
      static_cast<int64_t>(wtowers_parameters_.grid_size),
      static_cast<int64_t>(wtowers_parameters_.grid_size)};
  sdp_Mem *wrapped_dirty = sdp_mem_create_wrapper(
      (void *)dirty_image.data(), image_data_type, SDP_MEM_CPU, 2,
      image_shape.data(), nullptr, &status);
  if (status != SDP_SUCCESS) {
    throw std::runtime_error(
        "w-towers: Fatal error while wrapping memory for dirty image");
  }

  LogParameters();

  sdp_grid_wstack_wtower_grid_all(
      wrapped_visibilities, frequencies[0], frequency_step, wrapped_uvws,
      wtowers_parameters_.subgrid_size, wtowers_parameters_.grid_resolution,
      wtowers_parameters_.w_step, wtowers_parameters_.shear_u,
      wtowers_parameters_.shear_v, wtowers_parameters_.support,
      wtowers_parameters_.oversampling, wtowers_parameters_.w_support,
      wtowers_parameters_.w_oversampling, wtowers_parameters_.subgrid_frac,
      wtowers_parameters_.w_towers_height, wtowers_parameters_.verbosity,
      wrapped_dirty, n_threads_, &status);

  if (status != SDP_SUCCESS) {
    throw std::runtime_error(
        "w-towers: Failure inside sdp_grid_wstack_wtower_grid_all");
  }

  sdp_mem_free(wrapped_dirty);
  sdp_mem_free(wrapped_uvws);
  sdp_mem_free(wrapped_visibilities);

  aocommon::ImageBase<NumT>::Trim(
      dirty_image.data(), wtowers_parameters_.image_size,
      wtowers_parameters_.image_size, dirty_image.data(),
      wtowers_parameters_.grid_size, wtowers_parameters_.grid_size);
  for (size_t i = 0;
       i < wtowers_parameters_.image_size * wtowers_parameters_.image_size; ++i)
    image_[i] += dirty_image[i];
}

template <typename NumT>
void WTowersGridder<NumT>::FinalizeImage(double multiplication_factor) {
  for (auto &pix : image_) pix *= multiplication_factor;
}

template <typename NumT>
std::vector<float> WTowersGridder<NumT>::RealImage() {
  const size_t dx = (width_ - trimmed_width_) / 2;
  const size_t dy = (height_ - trimmed_height_) / 2;
  std::vector<float> image(width_ * height_,
                           std::numeric_limits<float>::quiet_NaN());
  for (size_t i = 0; i < trimmed_width_; ++i)
    for (size_t j = 0; j < trimmed_height_; ++j)
      image[(i + dx) + (j + dy) * width_] = image_[i * trimmed_height_ + j];
  return image;
}

template <typename NumT>
void WTowersGridder<NumT>::InitializePrediction(const float *image_data) {
  const size_t dx = (width_ - wtowers_parameters_.image_size) / 2;
  const size_t dy = (height_ - trimmed_height_) / 2;
  image_.resize(wtowers_parameters_.image_size * trimmed_height_);
  for (size_t i = 0; i < wtowers_parameters_.image_size; ++i)
    for (size_t j = 0; j < trimmed_height_; ++j)
      image_[i * trimmed_height_ + j] =
          image_data[(i + dx) + (j + dy) * width_];
}

template <typename NumT>
void WTowersGridder<NumT>::PredictVisibilities(
    size_t n_rows, size_t n_channels, const double *uvws,
    const double *frequencies, std::complex<float> *visibilities) const {
  const bool decreasing_freq =
      (n_channels > 1) && (frequencies[1] < frequencies[0]);
  if (decreasing_freq) {
    throw std::runtime_error(
        "W-towers does not currently support frequencies that aren't in "
        "ascending order\n");
  }

  const double frequency_step = frequencies[1] - frequencies[0];

  sdp_MemType image_data_type = SDP_MEM_FLOAT;
  if constexpr (std::is_same_v<NumT, double>) {
    image_data_type = SDP_MEM_DOUBLE;
  }

  sdp_Error status = SDP_SUCCESS;

  aocommon::ImageBase<NumT> untrimmed_image(wtowers_parameters_.grid_size,
                                            wtowers_parameters_.grid_size);
  aocommon::ImageBase<NumT>::Untrim(
      untrimmed_image.Data(), wtowers_parameters_.grid_size,
      wtowers_parameters_.grid_size, image_.data(),
      wtowers_parameters_.image_size, wtowers_parameters_.image_size);

  const std::vector<int64_t> visibilities_shape{
      static_cast<int64_t>(n_rows), static_cast<int64_t>(n_channels)};
  sdp_Mem *wrapped_visibilities = sdp_mem_create_wrapper(
      (void *)visibilities, SDP_MEM_COMPLEX_FLOAT, SDP_MEM_CPU, 2,
      visibilities_shape.data(), nullptr, &status);
  if (status != SDP_SUCCESS) {
    throw std::runtime_error(
        "w-towers: Fatal error while wrapping memory for visibilities");
  }
  const std::vector<int64_t> image_shape{
      static_cast<int64_t>(wtowers_parameters_.grid_size),
      static_cast<int64_t>(wtowers_parameters_.grid_size)};
  const std::vector<int64_t> uvws_shape{static_cast<int64_t>(n_rows), 3};
  sdp_Mem *wrapped_uvws =
      sdp_mem_create_wrapper((void *)uvws, SDP_MEM_DOUBLE, SDP_MEM_CPU, 2,
                             uvws_shape.data(), nullptr, &status);
  if (status != SDP_SUCCESS) {
    throw std::runtime_error(
        "w-towers: Fatal error while wrapping memory for uvws");
  }
  sdp_Mem *wrapped_dirty = sdp_mem_create_wrapper(
      (void *)untrimmed_image.Data(), image_data_type, SDP_MEM_CPU, 2,
      image_shape.data(), nullptr, &status);
  if (status != SDP_SUCCESS) {
    throw std::runtime_error(
        "w-towers: Fatal error while wrapping memory for dirty image");
  }

  LogParameters();

  sdp_grid_wstack_wtower_degrid_all(
      wrapped_dirty, frequencies[0], frequency_step, wrapped_uvws,
      wtowers_parameters_.subgrid_size, wtowers_parameters_.grid_resolution,
      wtowers_parameters_.w_step, wtowers_parameters_.shear_u,
      wtowers_parameters_.shear_v, wtowers_parameters_.support,
      wtowers_parameters_.oversampling, wtowers_parameters_.w_support,
      wtowers_parameters_.w_oversampling, wtowers_parameters_.subgrid_frac,
      wtowers_parameters_.w_towers_height, wtowers_parameters_.verbosity,
      wrapped_visibilities, n_threads_, &status);

  if (status != SDP_SUCCESS) {
    throw std::runtime_error(
        "w-towers: Failure inside sdp_grid_wstack_wtower_degrid_all");
  }

  sdp_mem_free(wrapped_dirty);
  sdp_mem_free(wrapped_uvws);
  sdp_mem_free(wrapped_visibilities);
}

}  // namespace wsclean

#endif  // #ifndef WSCLEAN_WTOWERS_GRIDDER_SIMPLE_IMPL_H_
