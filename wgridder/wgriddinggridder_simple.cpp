#include "wgriddinggridder_simple.h"
#include "ducc0/wgridder/wgridder.h"

using namespace ducc0;

WGriddingGridder_Simple::WGriddingGridder_Simple(
    size_t width, size_t height, size_t width_t, size_t height_t,
    double pixelSizeX, double pixelSizeY, double shiftL, double shiftM,
    size_t nthreads, double epsilon, size_t verbosity)
    : width_(width),
      height_(height),
      width_t_(width_t),
      height_t_(height_t),
      nthreads_(nthreads),
      pixelSizeX_(pixelSizeX),
      pixelSizeY_(pixelSizeY),
      shiftL_(shiftL),
      shiftM_(shiftM),
      epsilon_(epsilon),
      verbosity_(verbosity) {
  MR_assert(verbosity <= 2, "verbosity must be 0, 1, or 2");
}

void WGriddingGridder_Simple::memUsage(size_t &constant,
                                       size_t &per_vis) const {
  // storage for "grid": pessimistically assume an oversampling factor of 2
  constant = sigma_max * sigma_max * width_t_ * height_t_ *
             sizeof(std::complex<float>);
  // for prediction, we also need a copy of the dirty image
  constant += width_t_ * height_t_ * sizeof(float);  // trimmed dirty image
  // Storage for the indexing information is really hard to estimate ...
  // it can go up to 8 bytes per visibility, but this is a really pathological
  // scenario; should typically be below 1 byte/visibility
  per_vis = 8;  // overestimation, but the best we can do here
}

void WGriddingGridder_Simple::InitializeInversion() {
  img.assign(width_t_ * height_t_, 0);
}

void WGriddingGridder_Simple::AddInversionData(size_t nrows, size_t nchan,
                                               const double *uvw,
                                               const double *freq,
                                               const std::complex<float> *vis) {
  mav<double, 2> uvw2(uvw, {nrows, 3});
  mav<double, 1> freq2(freq, {nchan});
  mav<std::complex<float>, 2> ms(vis, {nrows, nchan});
  mav<float, 2> tdirty({width_t_, height_t_});
  mav<float, 2> twgt(nullptr, {0, 0}, false);
  mav<std::uint8_t, 2> tmask(nullptr, {0, 0}, false);
  ms2dirty<float, float>(uvw2, freq2, ms, twgt, tmask, pixelSizeX_, pixelSizeY_,
                         epsilon_, true, nthreads_, tdirty, verbosity_, true,
                         false, sigma_min, sigma_max, -shiftL_, -shiftM_);
  for (size_t i = 0; i < width_t_ * height_t_; ++i) img[i] += tdirty.craw(i);
}

void WGriddingGridder_Simple::FinalizeImage(double multiplicationFactor) {
  for (auto &pix : img) pix *= multiplicationFactor;
}

std::vector<float> WGriddingGridder_Simple::RealImage() {
  size_t dx = (width_ - width_t_) / 2;
  size_t dy = (height_ - height_t_) / 2;
  std::vector<float> image(width_ * height_, 0);
  for (size_t i = 0; i < width_t_; ++i)
    for (size_t j = 0; j < height_t_; ++j)
      image[(i + dx) + (j + dy) * width_] = img[i * height_t_ + j];
  return image;
}

void WGriddingGridder_Simple::InitializePrediction(const float *image_data) {
  size_t dx = (width_ - width_t_) / 2;
  size_t dy = (height_ - height_t_) / 2;
  img.resize(width_t_ * height_t_);
  for (size_t i = 0; i < width_t_; ++i)
    for (size_t j = 0; j < height_t_; ++j)
      img[i * height_t_ + j] = image_data[(i + dx) + (j + dy) * width_];
}

void WGriddingGridder_Simple::PredictVisibilities(
    size_t nrows, size_t nchan, const double *uvw, const double *freq,
    std::complex<float> *vis) const {
  mav<double, 2> uvw2(uvw, {nrows, 3});
  mav<double, 1> freq2(freq, {nchan});
  mav<std::complex<float>, 2> ms(vis, {nrows, nchan}, true);
  mav<float, 2> tdirty(img.data(), {width_t_, height_t_});
  mav<float, 2> twgt(nullptr, {0, 0}, false);
  mav<std::uint8_t, 2> tmask(nullptr, {0, 0}, false);
  dirty2ms<float, float>(uvw2, freq2, tdirty, twgt, tmask, pixelSizeX_,
                         pixelSizeY_, epsilon_, true, nthreads_, ms, verbosity_,
                         true, false, sigma_min, sigma_max, -shiftL_, -shiftM_);
}
