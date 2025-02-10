#ifndef WSCLEAN_WGRIDDER_SIMPLE_IMPL_H_
#define WSCLEAN_WGRIDDER_SIMPLE_IMPL_H_

#include "gridder_simple.h"

#include <complex>
#include <cstddef>
#include <vector>

#include "ducc0/wgridder/wgridder.h"
#include "ducc0/fft/fftnd_impl.h"

#include <LRUCache11.hpp>

#include "../gridding/msgridder.h"

using namespace ducc0;

/*
 * This file contains the implementation of various template methods from @ref
 * WGriddingGridder_Simple They are implemented here instead of directly in the
 * header or a single source file in order to be able to instantiate them in two
 * different source files. This is done because each instantiation must compile
 * the entire class as well as ducc0, which is quite resource intensive even for
 * a single compile. By breaking it into two source files, one for float and one
 * for double we keep time and memory requirements more manageable and allow for
 * better build parrelilisation.
 */

namespace wsclean {

template <typename NumT>
WGriddingGridder_Simple<NumT>::WGriddingGridder_Simple(
    size_t width, size_t height, size_t trimmed_width, size_t trimmed_height,
    double pixel_size_x, double pixel_size_y, double l_shift, double m_shift,
    size_t n_threads, double epsilon, size_t verbosity, bool tuning)
    : width_(width),
      height_(height),
      trimmed_width_(trimmed_width),
      trimmed_height_(trimmed_height),
      n_threads_(n_threads),
      pixel_size_x_(pixel_size_x),
      pixel_size_y_(pixel_size_y),
      l_shift_(l_shift),
      m_shift_(m_shift),
      epsilon_(epsilon),
      verbosity_(verbosity),
      tuning_(tuning) {
  MR_assert(verbosity <= 2, "verbosity must be 0, 1, or 2");
}

template <typename NumT>
size_t WGriddingGridder_Simple<NumT>::ConstantMemoryUsage() const {
  // Storage for "grid": pessimistically assume an oversampling factor of 2
  size_t constant = sigma_max * sigma_max * trimmed_width_ * trimmed_height_ *
                    sizeof(std::complex<float>);
  // For prediction, we also need a copy of the dirty image
  constant +=
      trimmed_width_ * trimmed_height_ * sizeof(NumT);  // trimmed dirty image
  return constant;
}

template <typename NumT>
size_t WGriddingGridder_Simple<NumT>::PerVisibilityMemoryUsage() const {
  // Storage for the indexing information is really hard to estimate ...
  // it can go up to 8 bytes per visibility, but this is a really pathological
  // scenario; should typically be below 1 byte/visibility
  return 8;  // Overestimation, but the best we can do here
}

template <typename NumT>
void WGriddingGridder_Simple<NumT>::InitializeInversion() {
  image_.assign(trimmed_width_ * trimmed_height_, 0);
}

template <typename NumT>
void WGriddingGridder_Simple<NumT>::AddInversionData(
    size_t n_rows, size_t n_channels, const double *uvws,
    const double *frequencies, const std::complex<float> *visibilities) {
  const bool decreasing_freq =
      (n_channels > 1) && (frequencies[1] < frequencies[0]);
  auto wrapped_frequencies(
      decreasing_freq
          ? cmav<double, 1>(frequencies + n_channels - 1, {n_channels}, {-1})
          : cmav<double, 1>(frequencies, {n_channels}));
  auto ms(
      decreasing_freq
          ? cmav<std::complex<float>, 2>(visibilities + n_channels - 1,
                                         {n_rows, n_channels},
                                         {ptrdiff_t(n_channels), -1})
          : cmav<std::complex<float>, 2>(visibilities, {n_rows, n_channels}));

  AddInversionMs(n_rows, uvws, wrapped_frequencies, ms);
}

template <typename NumT>
void WGriddingGridder_Simple<NumT>::AddInversionDataWithCorrectionCallback(
    GainMode mode, size_t n_polarizations, size_t n_rows, const double *uvws,
    const double *frequencies, size_t n_channels,
    const aocommon::BandData &selected_band,
    const std::pair<size_t, size_t> *antennas,
    const std::complex<float> *visibilities, const size_t *time_offsets,
    MsGridder *gridder, size_t n_antenna) {
  assert((selected_band.ChannelCount() <= 1) ||
         (frequencies[1] >= frequencies[0]));

  const cmav<double, 1> frequencies2(frequencies,
                                     {selected_band.ChannelCount()});
  // Construct a templated ms:
  //     VisibilityCallbackBuffer<mode, n_polarizations>
  // populated with visibilities and other data and call AddInversionMs(n_rows,
  // uvws, frequencies2, ms) on it.
  AddInversionMs(mode, n_polarizations, n_rows, uvws, std::ref(frequencies2),
                 n_channels, std::ref(selected_band), antennas, visibilities,
                 time_offsets, gridder, n_antenna);
}

template <typename NumT>
template <typename... Params>
void WGriddingGridder_Simple<NumT>::AddInversionMs(GainMode mode,
                                                   Params... params) {
  switch (mode) {
    case GainMode::kXX: {
      AddInversionMs<GainMode::kXX>(params...);
      break;
    }
    case GainMode::kYY: {
      AddInversionMs<GainMode::kYY>(params...);
      break;
    }
    case GainMode::k2VisDiagonal: {
      AddInversionMs<GainMode::k2VisDiagonal>(params...);
      break;
    }
    case GainMode::kTrace: {
      AddInversionMs<GainMode::kTrace>(params...);
      break;
    }
    case GainMode::kFull: {
      AddInversionMs<GainMode::kFull>(params...);
      break;
    }
  }
}

template <typename NumT>
void WGriddingGridder_Simple<NumT>::FinalizeImage(
    double multiplication_factor) {
  for (auto &pix : image_) pix *= multiplication_factor;
}

template <typename NumT>
std::vector<float> WGriddingGridder_Simple<NumT>::RealImage() {
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
void WGriddingGridder_Simple<NumT>::InitializePrediction(
    const float *image_data) {
  const size_t dx = (width_ - trimmed_width_) / 2;
  const size_t dy = (height_ - trimmed_height_) / 2;
  image_.resize(trimmed_width_ * trimmed_height_);
  for (size_t i = 0; i < trimmed_width_; ++i)
    for (size_t j = 0; j < trimmed_height_; ++j)
      image_[i * trimmed_height_ + j] =
          image_data[(i + dx) + (j + dy) * width_];
}

template <typename NumT>
void WGriddingGridder_Simple<NumT>::PredictVisibilities(
    size_t n_rows, size_t n_channels, const double *uvws,
    const double *frequencies, std::complex<float> *visibilities) const {
  cmav<double, 2> wrapped_uvws(uvws, {n_rows, 3});
  bool decreasing_freq = (n_channels > 1) && (frequencies[1] < frequencies[0]);
  auto wrapped_frequencies(
      decreasing_freq
          ? cmav<double, 1>(frequencies + n_channels - 1, {n_channels}, {-1})
          : cmav<double, 1>(frequencies, {n_channels}));
  auto ms(
      decreasing_freq
          ? vmav<std::complex<float>, 2>(visibilities + n_channels - 1,
                                         {n_rows, n_channels},
                                         {ptrdiff_t(n_channels), -1})
          : vmav<std::complex<float>, 2>(visibilities, {n_rows, n_channels}));
  cmav<NumT, 2> tdirty(image_.data(), {trimmed_width_, trimmed_height_});
  cmav<float, 2> twgt(nullptr, {0, 0});
  cmav<std::uint8_t, 2> tmask(nullptr, {0, 0});
  if (!tuning_)
    dirty2ms<NumT, NumT>(wrapped_uvws, wrapped_frequencies, tdirty, twgt, tmask,
                         pixel_size_x_, pixel_size_y_, epsilon_, true,
                         n_threads_, ms, verbosity_, true, false, sigma_min,
                         sigma_max, -l_shift_, -m_shift_);
  else
    dirty2ms_tuning<NumT, NumT>(wrapped_uvws, wrapped_frequencies, tdirty, twgt,
                                tmask, pixel_size_x_, pixel_size_y_, epsilon_,
                                true, n_threads_, ms, verbosity_, true, false,
                                sigma_min, sigma_max, -l_shift_, -m_shift_);
}

}  // namespace wsclean

#endif  // #ifndef WSCLEAN_WGRIDDER_SIMPLE_IMPL_H_
