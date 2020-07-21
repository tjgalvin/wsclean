#include "wgriddinggridder_simple.h"
#include "gridder_cxx.h"

using gridder::detail::CodeLocation;

using namespace std;

WGriddingGridder_Simple::WGriddingGridder_Simple(
    size_t width, size_t height, size_t width_t, size_t height_t,
    double pixelSizeX, double pixelSizeY, size_t nthreads, double epsilon,
    size_t verbosity)
    : width_(width),
      height_(height),
      width_t_(width_t),
      height_t_(height_t),
      nthreads_(nthreads),
      pixelSizeX_(pixelSizeX),
      pixelSizeY_(pixelSizeY),
      epsilon_(epsilon),
      verbosity_(verbosity) {
  myassert(verbosity <= 2, "verbosity must be 0, 1, or 2");
}

void WGriddingGridder_Simple::memUsage(size_t &constant,
                                       size_t &per_vis) const {
  constant = width_ * height_ * sizeof(complex<float>) *
                 2  // sometimes we have two arrays in memory
             + width_t_ * height_t_ * sizeof(float);  // trimmed dirty image
  per_vis = sizeof(gridder::detail::idx_t) *
            2;  // overestimation, but the best we can do here
}

void WGriddingGridder_Simple::InitializeInversion() {
  img.assign(width_t_ * height_t_, 0);
}

void WGriddingGridder_Simple::AddInversionData(size_t nrows, size_t nchan,
                                               const double *uvw,
                                               const double *freq,
                                               const complex<float> *vis) {
  gridder::const_mav<double, 2> uvw2(uvw, {nrows, 3});
  gridder::const_mav<double, 1> freq2(freq, {nchan});
  gridder::const_mav<complex<float>, 2> ms(vis, {nrows, nchan});
  vector<float> tmp(width_t_ * height_t_, 0);
  gridder::mav<float, 2> tdirty(tmp.data(), {width_t_, height_t_});
  gridder::const_mav<float, 2> twgt(nullptr, {0, 0});
  gridder::ms2dirty_general(uvw2, freq2, ms, twgt, pixelSizeX_, pixelSizeY_,
                            width_, height_, epsilon_, true, nthreads_, tdirty,
                            verbosity_, true);
  for (size_t i = 0; i < width_t_ * height_t_; ++i) img[i] += tmp[i];
}

void WGriddingGridder_Simple::FinalizeImage(double multiplicationFactor,
                                            bool correctFFTFactor) {
  if (correctFFTFactor) multiplicationFactor /= sqrt(width_ * height_);
  for (auto &pix : img) pix *= multiplicationFactor;
}

vector<float> WGriddingGridder_Simple::RealImage() {
  size_t dx = (width_ - width_t_) / 2;
  size_t dy = (height_ - height_t_) / 2;
  vector<float> image(width_ * height_, 0);
  for (size_t i = 0; i < width_t_; ++i)
    for (size_t j = 0; j < height_t_; ++j)
      image[(i + dx) + (j + dy) * width_] = img[i * height_t_ + j];
  return image;
}

void WGriddingGridder_Simple::InitializePrediction(vector<float> &&image) {
  size_t dx = (width_ - width_t_) / 2;
  size_t dy = (height_ - height_t_) / 2;
  img.resize(width_t_ * height_t_);
  myassert(image.size() == width_ * height_, "bad image dimensions");
  for (size_t i = 0; i < width_t_; ++i)
    for (size_t j = 0; j < height_t_; ++j)
      img[i * height_t_ + j] = image[(i + dx) + (j + dy) * width_];
}

void WGriddingGridder_Simple::PredictVisibilities(
    size_t nrows, size_t nchan, const double *uvw, const double *freq,
    std::complex<float> *vis) const {
  gridder::const_mav<double, 2> uvw2(uvw, {nrows, 3});
  gridder::const_mav<double, 1> freq2(freq, {nchan});
  gridder::mav<complex<float>, 2> ms(vis, {nrows, nchan});
  gridder::const_mav<float, 2> tdirty(img.data(), {width_t_, height_t_});
  gridder::const_mav<float, 2> twgt(nullptr, {0, 0});
  gridder::dirty2ms_general(uvw2, freq2, tdirty, twgt, pixelSizeX_, pixelSizeY_,
                            width_, height_, epsilon_, true, nthreads_, ms,
                            verbosity_, true);
}
