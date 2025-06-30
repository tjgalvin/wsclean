#ifndef WSCLEAN_WGRIDDER_IMPL_H_
#define WSCLEAN_WGRIDDER_IMPL_H_

#include "wgridder.h"

#include <complex>
#include <cstddef>
#include <vector>

#include "ducc0/wgridder/wgridder.h"
#include "ducc0/fft/fftnd_impl.h"

#include "../gridding/msgridder.h"

using namespace ducc0;

/*
 * This file contains the implementation of various template methods from @ref
 * WGridder. They are implemented here instead of directly in wgridder.h or
 * wgridder.cpp in order to be able to instantiate them in two different source
 * files.
 * This is done when wgridder_double.cpp and wgridder_float.cpp include this
 * file.
 * This is done because each instantiation must compile the entire class as well
 * as ducc0, which is quite resource intensive.
 * By breaking compilation in two, time and memory requirements are kept more
 * manageable and the build can make better use of parallel processing.
 */

namespace wsclean {

template <typename NumT>
WGridder<NumT>::WGridder(size_t width, size_t height, size_t trimmed_width,
                         size_t trimmed_height, double pixel_size_x,
                         double pixel_size_y, double l_shift, double m_shift,
                         size_t n_threads, double epsilon, size_t verbosity,
                         bool tuning)
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
size_t WGridder<NumT>::ConstantMemoryUsage() const {
  // Storage for "grid": pessimistically assume an oversampling factor of 2
  size_t constant = sigma_max * sigma_max * trimmed_width_ * trimmed_height_ *
                    sizeof(std::complex<float>);
  // For prediction, we also need a copy of the dirty image
  constant +=
      trimmed_width_ * trimmed_height_ * sizeof(NumT);  // trimmed dirty image
  return constant;
}

template <typename NumT>
size_t WGridder<NumT>::PerVisibilityMemoryUsage() const {
  // Storage for the indexing information is really hard to estimate ...
  // it can go up to 8 bytes per visibility, but this is a really pathological
  // scenario; should typically be below 1 byte/visibility
  return 8;  // Overestimation, but the best we can do here
}

template <typename NumT>
void WGridder<NumT>::InitializeInversion() {
  image_.assign(trimmed_width_ * trimmed_height_, 0);
}

template <typename NumT>
void WGridder<NumT>::AddInversionData(size_t n_rows, size_t n_channels,
                                      const double *uvws,
                                      const double *frequencies,
                                      const std::complex<float> *visibilities) {
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
void WGridder<NumT>::AddInversionDataWithCorrectionCallback(
    GainMode mode, size_t n_polarizations, size_t n_rows, const double *uvws,
    const double *frequencies, VisibilityCallbackData &data) {
  assert((data.selected_band.ChannelCount() <= 1) ||
         (frequencies[1] >= frequencies[0]));

  const cmav<double, 1> frequencies2(frequencies,
                                     {data.selected_band.ChannelCount()});

  MsGridder *gridder = data.gridder;
  const size_t n_parms = gridder->NumValuesPerSolution();

  // Construct a templated ms:
  //     VisibilityCallbackBuffer<mode, n_polarizations>
  // populated with visibilities and other data and call
  // CreateAndAddInversionMs(n_rows, uvws, frequencies2, ms) on it.
  const bool apply_beam = gridder->WillApplyBeam();
  const bool apply_forward =
      gridder->GetPsfMode() == PsfMode::kDirectionDependent;
  const bool has_h5_parm = gridder->GetVisibilityModifier().HasH5Parm();
  CreateAndAddInversionMs(mode, n_polarizations, n_parms, apply_beam,
                          apply_forward, has_h5_parm, n_rows, uvws,
                          frequencies2, data);
}

template <typename NumT>
void WGridder<NumT>::CreateAndAddInversionMs(
    GainMode mode, size_t n_polarizations, size_t n_parms, bool apply_beam,
    bool apply_forward, bool has_h5_parm, size_t n_rows, const double *uvws,
    const ducc0::cmav<double, 1> &frequencies, VisibilityCallbackData &data) {
  switch (mode) {
    case GainMode::kXX: {
      CreateAndAddInversionMs2<GainMode::kXX>(
          n_polarizations, n_parms, apply_beam, apply_forward, has_h5_parm,
          n_rows, uvws, frequencies, data);
      break;
    }
    case GainMode::kYY: {
      CreateAndAddInversionMs2<GainMode::kYY>(
          n_polarizations, n_parms, apply_beam, apply_forward, has_h5_parm,
          n_rows, uvws, frequencies, data);
      break;
    }
    case GainMode::k2VisDiagonal: {
      CreateAndAddInversionMs2<GainMode::k2VisDiagonal>(
          n_polarizations, n_parms, apply_beam, apply_forward, has_h5_parm,
          n_rows, uvws, frequencies, data);
      break;
    }
    case GainMode::kTrace: {
      CreateAndAddInversionMs2<GainMode::kTrace>(
          n_polarizations, n_parms, apply_beam, apply_forward, has_h5_parm,
          n_rows, uvws, frequencies, data);
      break;
    }
    case GainMode::kFull: {
      CreateAndAddInversionMs2<GainMode::kFull>(
          n_polarizations, n_parms, apply_beam, apply_forward, has_h5_parm,
          n_rows, uvws, frequencies, data);
      break;
    }
  }
}

template <typename NumT>
template <GainMode Mode>
void WGridder<NumT>::CreateAndAddInversionMs2(
    size_t n_polarizations, size_t n_parms, bool apply_beam, bool apply_forward,
    bool has_h5_parm, size_t n_rows, const double *uvws,
    const ducc0::cmav<double, 1> &frequencies, VisibilityCallbackData &data) {
  switch (n_polarizations) {
    case 1: {
      if constexpr (GetNVisibilities(Mode) == 1) {
        CreateAndAddInversionMs3<Mode, 1>(n_parms, apply_beam, apply_forward,
                                          has_h5_parm, n_rows, uvws,
                                          frequencies, data);
      }
      break;
    }
    case 2: {
      if constexpr (GetNVisibilities(Mode) == 2) {
        CreateAndAddInversionMs3<Mode, 2>(n_parms, apply_beam, apply_forward,
                                          has_h5_parm, n_rows, uvws,
                                          frequencies, data);
      }
      break;
    }
    case 4: {
      if constexpr (GetNVisibilities(Mode) == 4) {
        CreateAndAddInversionMs3<Mode, 4>(n_parms, apply_beam, apply_forward,
                                          has_h5_parm, n_rows, uvws,
                                          frequencies, data);
      }
      break;
    }
    default:
      assert(false);
  }
}

template <typename NumT>
template <GainMode Mode, size_t NPolarizations>
void WGridder<NumT>::CreateAndAddInversionMs3(
    size_t n_parms, bool apply_beam, bool apply_forward, bool has_h5_parm,
    size_t n_rows, const double *uvws,
    const ducc0::cmav<double, 1> &frequencies, VisibilityCallbackData &data) {
  if (n_parms == 2) {
    CreateAndAddInversionMs4<Mode, NPolarizations, 2>(apply_beam, apply_forward,
                                                      has_h5_parm, n_rows, uvws,
                                                      frequencies, data);
  } else {
    CreateAndAddInversionMs4<Mode, NPolarizations, 4>(apply_beam, apply_forward,
                                                      has_h5_parm, n_rows, uvws,
                                                      frequencies, data);
  }
}

template <typename NumT>
template <GainMode Mode, size_t NPolarizations, size_t NParms>
void WGridder<NumT>::CreateAndAddInversionMs4(
    bool apply_beam, bool apply_forward, bool has_h5_parm, size_t n_rows,
    const double *uvws, const ducc0::cmav<double, 1> &frequencies,
    VisibilityCallbackData &data) {
  if (apply_beam) {
    CreateAndAddInversionMs5<Mode, NPolarizations, NParms, true>(
        apply_forward, has_h5_parm, n_rows, uvws, frequencies, data);
  } else {
    CreateAndAddInversionMs5<Mode, NPolarizations, NParms, false>(
        apply_forward, has_h5_parm, n_rows, uvws, frequencies, data);
  }
}

template <typename NumT>
template <GainMode Mode, size_t NPolarizations, size_t NParms, bool ApplyBeam>
void WGridder<NumT>::CreateAndAddInversionMs5(
    bool apply_forward, bool has_h5_parm, size_t n_rows, const double *uvws,
    const ducc0::cmav<double, 1> &frequencies, VisibilityCallbackData &data) {
  if (apply_forward) {
    CreateAndAddInversionMs6<Mode, NPolarizations, NParms, ApplyBeam, true>(
        has_h5_parm, n_rows, uvws, frequencies, data);
  } else {
    CreateAndAddInversionMs6<Mode, NPolarizations, NParms, ApplyBeam, false>(
        has_h5_parm, n_rows, uvws, frequencies, data);
  }
}

template <typename NumT>
template <GainMode Mode, size_t NPolarizations, size_t NParms, bool ApplyBeam,
          bool ApplyForward>
void WGridder<NumT>::CreateAndAddInversionMs6(
    bool has_h5_parm, size_t n_rows, const double *uvws,
    const ducc0::cmav<double, 1> &frequencies, VisibilityCallbackData &data) {
  if (has_h5_parm) {
    CreateAndAddInversionMs7<Mode, NPolarizations, NParms, ApplyBeam,
                             ApplyForward, true>(n_rows, uvws, frequencies,
                                                 data);
  } else {
    CreateAndAddInversionMs7<Mode, NPolarizations, NParms, ApplyBeam,
                             ApplyForward, false>(n_rows, uvws, frequencies,
                                                  data);
  }
}

template <typename NumT>
template <GainMode Mode, size_t NPolarizations, size_t NParms, bool ApplyBeam,
          bool ApplyForward, bool HasH5Parm>
void WGridder<NumT>::CreateAndAddInversionMs7(
    size_t n_rows, const double *uvws,
    const ducc0::cmav<double, 1> &frequencies, VisibilityCallbackData &data) {
  if (data.gridder->GetPsfMode() == PsfMode::kDirectionDependent) {
    CreateAndAddInversionMs8<Mode, NPolarizations, NParms, ApplyBeam,
                             ApplyForward, HasH5Parm, true>(n_rows, uvws,
                                                            frequencies, data);
  } else {
    CreateAndAddInversionMs8<Mode, NPolarizations, NParms, ApplyBeam,
                             ApplyForward, HasH5Parm, false>(n_rows, uvws,
                                                             frequencies, data);
  }
}

template <typename NumT>
template <GainMode Mode, size_t NPolarizations, size_t NParms, bool ApplyBeam,
          bool ApplyForward, bool HasH5Parm, bool ApplyRotation>
void WGridder<NumT>::CreateAndAddInversionMs8(
    size_t n_rows, const double *uvws,
    const ducc0::cmav<double, 1> &frequencies, VisibilityCallbackData &data) {
  const std::function visibility_callback =
      internal::VisibilityCallback<Mode, NPolarizations, NParms, ApplyBeam,
                                   ApplyForward, HasH5Parm, ApplyRotation>;
  const VisibilityCallbackBuffer<std::complex<float>> ms(n_rows, data,
                                                         visibility_callback);
  AddInversionMs(n_rows, uvws, frequencies, ms);
}

template <typename NumT>
void WGridder<NumT>::FinalizeImage(double multiplication_factor) {
  for (auto &pix : image_) pix *= multiplication_factor;
}

template <typename NumT>
std::vector<float> WGridder<NumT>::RealImage() {
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
void WGridder<NumT>::InitializePrediction(const float *image_data) {
  const size_t dx = (width_ - trimmed_width_) / 2;
  const size_t dy = (height_ - trimmed_height_) / 2;
  image_.resize(trimmed_width_ * trimmed_height_);
  for (size_t i = 0; i < trimmed_width_; ++i)
    for (size_t j = 0; j < trimmed_height_; ++j)
      image_[i * trimmed_height_ + j] =
          image_data[(i + dx) + (j + dy) * width_];
}

template <typename NumT>
void WGridder<NumT>::PredictVisibilities(
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

#endif  // #ifndef WSCLEAN_WGRIDDER_IMPL_H_
