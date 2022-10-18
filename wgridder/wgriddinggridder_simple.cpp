#include "wgriddinggridder_simple.h"
#include "ducc0/wgridder/wgridder.h"

using namespace ducc0;

template <typename NumT>
WGriddingGridder_Simple<NumT>::WGriddingGridder_Simple(
    size_t width, size_t height, size_t width_t, size_t height_t,
    double pixelSizeX, double pixelSizeY, double l_shift, double m_shift,
    size_t nthreads, double epsilon, size_t verbosity, bool tuning)
    : width_(width),
      height_(height),
      width_t_(width_t),
      height_t_(height_t),
      nthreads_(nthreads),
      pixelSizeX_(pixelSizeX),
      pixelSizeY_(pixelSizeY),
      l_shift_(l_shift),
      m_shift_(m_shift),
      epsilon_(epsilon),
      verbosity_(verbosity),
      tuning_(tuning) {
  MR_assert(verbosity <= 2, "verbosity must be 0, 1, or 2");
}

template <typename NumT>
void WGriddingGridder_Simple<NumT>::MemoryUsage(size_t &constant,
                                                size_t &per_vis) const {
  // storage for "grid": pessimistically assume an oversampling factor of 2
  constant = sigma_max * sigma_max * width_t_ * height_t_ *
             sizeof(std::complex<float>);
  // for prediction, we also need a copy of the dirty image
  constant += width_t_ * height_t_ * sizeof(NumT);  // trimmed dirty image
  // Storage for the indexing information is really hard to estimate ...
  // it can go up to 8 bytes per visibility, but this is a really pathological
  // scenario; should typically be below 1 byte/visibility
  per_vis = 8;  // overestimation, but the best we can do here
}

template <typename NumT>
void WGriddingGridder_Simple<NumT>::InitializeInversion() {
  img.assign(width_t_ * height_t_, 0);
}

template <typename NumT>
void WGriddingGridder_Simple<NumT>::AddInversionData(
    size_t nrows, size_t nchan, const double *uvw, const double *freq,
    const std::complex<float> *vis) {
  cmav<double, 2> uvw2(uvw, {nrows, 3});
  cmav<double, 1> freq2(freq, {nchan});
  cmav<std::complex<float>, 2> ms(vis, {nrows, nchan});
  vmav<NumT, 2> tdirty({width_t_, height_t_});
  cmav<float, 2> twgt(nullptr, {0, 0});
  cmav<std::uint8_t, 2> tmask(nullptr, {0, 0});
  if (!tuning_)
    ms2dirty<NumT, NumT>(uvw2, freq2, ms, twgt, tmask, pixelSizeX_, pixelSizeY_,
                         epsilon_, true, nthreads_, tdirty, verbosity_, true,
                         false, sigma_min, sigma_max, -l_shift_, -m_shift_);
  else
    ms2dirty_tuning<NumT, NumT>(uvw2, freq2, ms, twgt, tmask, pixelSizeX_,
                                pixelSizeY_, epsilon_, true, nthreads_, tdirty,
                                verbosity_, true, false, sigma_min, sigma_max,
                                -l_shift_, -m_shift_);
  for (size_t i = 0; i < width_t_ * height_t_; ++i) img[i] += tdirty.raw(i);
}

template <typename NumT>
void WGriddingGridder_Simple<NumT>::FinalizeImage(double multiplicationFactor) {
  for (auto &pix : img) pix *= multiplicationFactor;
}

template <typename NumT>
std::vector<float> WGriddingGridder_Simple<NumT>::RealImage() {
  size_t dx = (width_ - width_t_) / 2;
  size_t dy = (height_ - height_t_) / 2;
  std::vector<float> image(width_ * height_,
                           std::numeric_limits<float>::quiet_NaN());
  for (size_t i = 0; i < width_t_; ++i)
    for (size_t j = 0; j < height_t_; ++j)
      image[(i + dx) + (j + dy) * width_] = img[i * height_t_ + j];
  return image;
}

template <typename NumT>
void WGriddingGridder_Simple<NumT>::InitializePrediction(
    const float *image_data) {
  size_t dx = (width_ - width_t_) / 2;
  size_t dy = (height_ - height_t_) / 2;
  img.resize(width_t_ * height_t_);
  for (size_t i = 0; i < width_t_; ++i)
    for (size_t j = 0; j < height_t_; ++j)
      img[i * height_t_ + j] = image_data[(i + dx) + (j + dy) * width_];
}

template <typename NumT>
void WGriddingGridder_Simple<NumT>::PredictVisibilities(
    size_t nrows, size_t nchan, const double *uvw, const double *freq,
    std::complex<float> *vis) const {
  cmav<double, 2> uvw2(uvw, {nrows, 3});
  cmav<double, 1> freq2(freq, {nchan});
  vmav<std::complex<float>, 2> ms(vis, {nrows, nchan});
  cmav<NumT, 2> tdirty(img.data(), {width_t_, height_t_});
  cmav<float, 2> twgt(nullptr, {0, 0});
  cmav<std::uint8_t, 2> tmask(nullptr, {0, 0});
  if (!tuning_)
    dirty2ms<NumT, NumT>(uvw2, freq2, tdirty, twgt, tmask, pixelSizeX_,
                         pixelSizeY_, epsilon_, true, nthreads_, ms, verbosity_,
                         true, false, sigma_min, sigma_max, -l_shift_,
                         -m_shift_);
  else
    dirty2ms_tuning<NumT, NumT>(uvw2, freq2, tdirty, twgt, tmask, pixelSizeX_,
                                pixelSizeY_, epsilon_, true, nthreads_, ms,
                                verbosity_, true, false, sigma_min, sigma_max,
                                -l_shift_, -m_shift_);
}

template class WGriddingGridder_Simple<float>;
template class WGriddingGridder_Simple<double>;
