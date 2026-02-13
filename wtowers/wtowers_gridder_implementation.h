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
WTowersGridder<NumT>::WTowersGridder(
    size_t width, size_t height, size_t trimmed_width, size_t trimmed_height,
    double pixel_size_x, double pixel_size_y, double l_shift, double m_shift,
    int subgrid_size, int support, int w_support, double padding,
    double w_padding, size_t n_threads, double accuracy, double max_abs_w,
    size_t verbosity)
    : width_(width),
      height_(height),
      trimmed_width_(trimmed_width),
      trimmed_height_(trimmed_height),
      n_threads_(n_threads),
      pixel_size_x_(pixel_size_x),
      pixel_size_y_(pixel_size_y),
      l_shift_(l_shift),
      m_shift_(m_shift),
      subgrid_size_(subgrid_size),
      support_(support),
      w_support_(w_support),
      padding_factor_(padding),
      padding_factor_w_(w_padding),
      verbosity_(verbosity > 0 ? 1 : 0) {
  assert(verbosity <= 2);

  /// Various configuration paramaters for W-towers, we calculate a few of them
  /// and hardcode the rest. Future implementations should try to choose optimal
  /// parameters automatically based on desired precision (likely need to do
  /// systematic tests and compile a table)

  // Determine projection to use for the image we are generating
  double cell_size_l = std::sin(pixel_size_x);
  double cell_size_m = std::sin(pixel_size_y);
  int image_size_l =
      2 * sdp_fft_padded_size(int(trimmed_width / 2), padding_factor_);
  int image_size_m =
      2 * sdp_fft_padded_size(int(trimmed_height / 2), padding_factor_);
  wtowers_parameters_.projection =
      sdp_GridProjection(image_size_l, image_size_m, cell_size_l * image_size_l,
                         cell_size_m * image_size_m, -l_shift, m_shift,
                         0,    // Overwritten below
                         0, 0  // Shears not supported by WSClean. One day!
      );

  // Derive how we're dealing with the w-axis.
  const double fov_l = cell_size_l * trimmed_width;
  const double fov_m = cell_size_m * trimmed_height;
  wtowers_parameters_.field_of_view_l = fov_l;
  wtowers_parameters_.field_of_view_m = fov_m;
  wtowers_parameters_.projection.shift_n = sdp_gridder_determine_n_shift(
      wtowers_parameters_.projection, fov_l, fov_m);
  wtowers_parameters_.w_step = sdp_gridder_determine_w_step(
      wtowers_parameters_.projection, fov_l, fov_m, 1. / padding_factor_w_);

  // Determine maximum w-tower height (i.e. at what point it needs to fall
  // back to w-stacking). This does a bunch of test-imaging internally. To
  // make this go fast we reduce the test image size to about twice the
  // subgrid size. The result is typically that the w-tower-height gets
  // somewhat underestimated, but not to a degree that really matters.
  // On the flip side, this check can currently be a bit slow...
  sdp_Error status = SDP_SUCCESS;
  sdp_GridProjection test_projection = wtowers_parameters_.projection;
  test_projection.size_l = test_projection.size_m = 2 * subgrid_size_;

  wtowers_parameters_.w_towers_height =
      sdp_gridder_determine_max_w_tower_height(
          test_projection, subgrid_size_, wtowers_parameters_.w_step, max_abs_w,
          support_, wtowers_parameters_.oversampling, w_support_,
          wtowers_parameters_.w_oversampling, fov_l, fov_m,
          wtowers_parameters_.subgrid_frac, 3, accuracy, &status);
  if (status != SDP_SUCCESS) {
    throw std::runtime_error("Error computing w_towers_height");
  }

  // This used to return zero due to a semi-rare oversampling bug that was
  // tickled by the height calculation... Better safe than sorry. Also if it
  // happens legitimately likely want to give a hint what needs changing.
  // An alternative strategy here would actually be to reduce w_step. But
  // I have yet to see a realistic imaging configuration that would trip this.
  if (wtowers_parameters_.w_towers_height == 0) {
    throw std::runtime_error(
        "Calculated w_towers_height not enough to cover "
        "w support. Consider bigger subgrids?");
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
  Logger::Debug << "field_of_view:  " << wtowers_parameters_.field_of_view_l
                << ", " << wtowers_parameters_.field_of_view_m << "\n";
  Logger::Debug << "image_size:     " << wtowers_parameters_.projection.size_l
                << " x " << wtowers_parameters_.projection.size_m << "\n";
  Logger::Debug << "subgrid_size:   " << subgrid_size_ << "\n";
  Logger::Debug << "grid_resolution / theta: "
                << wtowers_parameters_.projection.theta_l << ", "
                << wtowers_parameters_.projection.theta_m << "\n";
  double cell_size_l = wtowers_parameters_.projection.theta_l /
                       wtowers_parameters_.projection.size_l;
  double cell_size_m = wtowers_parameters_.projection.theta_m /
                       wtowers_parameters_.projection.size_m;
  Logger::Debug << "cell_size_rad:  " << cell_size_l << ", " << cell_size_m
                << "\n";
  Logger::Debug << "w_step: " << wtowers_parameters_.w_step << "\n";
  Logger::Debug << "w_towers_height: " << wtowers_parameters_.w_towers_height
                << "\n";
  Logger::Debug << "support:        " << support_ << "\n";
  Logger::Debug << "oversampling:   " << wtowers_parameters_.oversampling
                << "\n";
  Logger::Debug << "w_support:      " << w_support_ << "\n";
  Logger::Debug << "w_oversampling: " << wtowers_parameters_.w_oversampling
                << "\n";
  Logger::Debug << "subgrid_frac:   " << wtowers_parameters_.subgrid_frac
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
  dirty_image.assign(wtowers_parameters_.projection.size_l *
                         wtowers_parameters_.projection.size_m,
                     0);

  const int64_t visibilities_shape[2] = {static_cast<int64_t>(n_rows),
                                         static_cast<int64_t>(n_channels)};
  sdp_Mem *wrapped_visibilities =
      sdp_mem_create_wrapper(visibilities, SDP_MEM_COMPLEX_FLOAT, SDP_MEM_CPU,
                             2, visibilities_shape, nullptr, &status);
  if (status != SDP_SUCCESS) {
    throw std::runtime_error(
        "w-towers: Fatal error while wrapping memory for visibilities");
  }
  const int64_t uvws_shape[2] = {static_cast<int64_t>(n_rows), 3};
  sdp_Mem *wrapped_uvws = sdp_mem_create_wrapper(
      uvws, SDP_MEM_DOUBLE, SDP_MEM_CPU, 2, uvws_shape, nullptr, &status);
  if (status != SDP_SUCCESS) {
    throw std::runtime_error(
        "w-towers: Fatal error while wrapping memory for uvws");
  }
  const int64_t image_shape[2] = {
      static_cast<int64_t>(wtowers_parameters_.projection.size_l),
      static_cast<int64_t>(wtowers_parameters_.projection.size_m)};
  const int64_t image_strides[2] = {sdp_mem_type_size(image_data_type),
                                    sdp_mem_type_size(image_data_type) *
                                        wtowers_parameters_.projection.size_l};
  sdp_Mem *wrapped_dirty =
      sdp_mem_create_wrapper(dirty_image.data(), image_data_type, SDP_MEM_CPU,
                             2, image_shape, image_strides, &status);
  if (status != SDP_SUCCESS) {
    throw std::runtime_error(
        "w-towers: Fatal error while wrapping memory for dirty image");
  }

  LogParameters();

  sdp_grid_wstack_wtower_grid_all(
      wrapped_visibilities, frequencies[0], frequency_step, wrapped_uvws,
      wtowers_parameters_.projection, subgrid_size_, wtowers_parameters_.w_step,
      support_, wtowers_parameters_.oversampling, w_support_,
      wtowers_parameters_.w_oversampling, wtowers_parameters_.subgrid_frac,
      wtowers_parameters_.w_towers_height, verbosity_, wrapped_dirty,
      n_threads_, &status);

  if (status != SDP_SUCCESS) {
    throw std::runtime_error(
        "w-towers: Failure inside sdp_grid_wstack_wtower_grid_all");
  }

  sdp_mem_free(wrapped_dirty);
  sdp_mem_free(wrapped_uvws);
  sdp_mem_free(wrapped_visibilities);

  aocommon::ImageBase<NumT>::Trim(dirty_image.data(), trimmed_width_,
                                  trimmed_height_, dirty_image.data(),
                                  wtowers_parameters_.projection.size_l,
                                  wtowers_parameters_.projection.size_m);
  for (size_t i = 0; i < trimmed_width_ * trimmed_height_; ++i)
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
  for (size_t j = 0; j < trimmed_height_; ++j)
    for (size_t i = 0; i < trimmed_width_; ++i)
      image[(i + dx) + (j + dy) * width_] = image_[i + j * trimmed_width_];
  return image;
}

template <typename NumT>
void WTowersGridder<NumT>::InitializePrediction(const float *image_data) {
  const size_t dx = (width_ - trimmed_width_) / 2;
  const size_t dy = (height_ - trimmed_height_) / 2;
  image_.resize(trimmed_width_ * trimmed_height_);
  for (size_t j = 0; j < trimmed_height_; ++j)
    for (size_t i = 0; i < trimmed_width_; ++i)
      image_[i + j * trimmed_width_] = image_data[(i + dx) + (j + dy) * width_];
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

  aocommon::ImageBase<NumT> untrimmed_image(
      wtowers_parameters_.projection.size_l,
      wtowers_parameters_.projection.size_m);
  aocommon::ImageBase<NumT>::Untrim(
      untrimmed_image.Data(), wtowers_parameters_.projection.size_l,
      wtowers_parameters_.projection.size_m, image_.data(), trimmed_width_,
      trimmed_height_);

  const std::vector<int64_t> visibilities_shape{
      static_cast<int64_t>(n_rows), static_cast<int64_t>(n_channels)};
  sdp_Mem *wrapped_visibilities =
      sdp_mem_create_wrapper(visibilities, SDP_MEM_COMPLEX_FLOAT, SDP_MEM_CPU,
                             2, visibilities_shape.data(), nullptr, &status);
  if (status != SDP_SUCCESS) {
    throw std::runtime_error(
        "w-towers: Fatal error while wrapping memory for visibilities");
  }
  const int64_t uvws_shape[2] = {static_cast<int64_t>(n_rows), 3};
  sdp_Mem *wrapped_uvws = sdp_mem_create_wrapper(
      uvws, SDP_MEM_DOUBLE, SDP_MEM_CPU, 2, uvws_shape, nullptr, &status);
  if (status != SDP_SUCCESS) {
    throw std::runtime_error(
        "w-towers: Fatal error while wrapping memory for uvws");
  }
  const int64_t image_shape[2] = {
      static_cast<int64_t>(wtowers_parameters_.projection.size_l),
      static_cast<int64_t>(wtowers_parameters_.projection.size_m)};
  const int64_t image_strides[2] = {sdp_mem_type_size(image_data_type),
                                    sdp_mem_type_size(image_data_type) *
                                        wtowers_parameters_.projection.size_l};
  sdp_Mem *wrapped_dirty = sdp_mem_create_wrapper(
      untrimmed_image.Data(), image_data_type, SDP_MEM_CPU, 2, image_shape,
      image_strides, &status);
  if (status != SDP_SUCCESS) {
    throw std::runtime_error(
        "w-towers: Fatal error while wrapping memory for dirty image");
  }

  LogParameters();

  sdp_grid_wstack_wtower_degrid_all(
      wrapped_dirty, wtowers_parameters_.projection, frequencies[0],
      frequency_step, wrapped_uvws, subgrid_size_, wtowers_parameters_.w_step,
      support_, wtowers_parameters_.oversampling, w_support_,
      wtowers_parameters_.w_oversampling, wtowers_parameters_.subgrid_frac,
      wtowers_parameters_.w_towers_height, verbosity_, wrapped_visibilities,
      n_threads_, &status);

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
