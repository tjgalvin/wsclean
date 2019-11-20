#include "wgriddinggridder_simple.h"
#include "gridder_cxx.h"

using gridder::detail::CodeLocation;

using namespace std;

WGriddingGridder_Simple::WGriddingGridder_Simple(size_t width, size_t height,
  double pixelSizeX, double pixelSizeY, size_t nthreads, double epsilon)
  : width_(width), height_(height), nthreads_(nthreads),
    pixelSizeX_(pixelSizeX), pixelSizeY_(pixelSizeY), epsilon_(epsilon) {}

void WGriddingGridder_Simple::memUsage(size_t &constant, size_t &per_vis) const
  {
  constant = pixelSizeX_*pixelSizeY_
             *4 // padding
             *sizeof(complex<float>)
             *2 // sometimes we have two arrays in memory
           + pixelSizeX_*pixelSizeY_
             *sizeof(float); // dirty image
  per_vis = sizeof(gridder::detail::idx_t)*2; // overestimation, but the best we can do here
  }

void WGriddingGridder_Simple::InitializeInversion()
  {
  img.resize(pixelSizeX_*pixelSizeY_,0);
  }

void WGriddingGridder_Simple::AddInversionData(size_t nrows, size_t nchan,
  const double *uvw, const double *freq, const complex<float> *vis)
  {
  gridder::const_mav<double,2> uvw2(uvw, {nrows, 3});
  gridder::const_mav<double,1> freq2(freq, {nchan});
  gridder::const_mav<complex<float>,2> ms(vis, {nrows,nchan});
  vector<float> tmp(width_*height_, 0);
  gridder::mav<float,2> tdirty(tmp.data(), {width_, height_});
  gridder::const_mav<float,2> twgt(nullptr, {0,0});
  gridder::ms2dirty(uvw2, freq2, ms, twgt,
    pixelSizeX_, pixelSizeY_, epsilon_, true, nthreads_, tdirty, 0);
  for (size_t i=0; i<width_*height_; ++i)
    img[i] += tmp[i];
  }

void WGriddingGridder_Simple::FinalizeImage(double multiplicationFactor,
  bool correctFFTFactor)
  {
  if (correctFFTFactor)
    multiplicationFactor /= sqrt(4*width_*height_);
  for (auto &pix: img)
    pix *= multiplicationFactor;
  }

vector<float> WGriddingGridder_Simple::RealImage()
  {
  return move(img);
  }

void WGriddingGridder_Simple::InitializePrediction(vector<float> &&image)
  {
  img=image;
  myassert(img.size()==width_*height_, "bad image dimensions");
  }

void WGriddingGridder_Simple::PredictVisibilities(size_t nrows, size_t nchan,
  const double *uvw, const double *freq,
  std::complex<float> *vis) const
  {
  gridder::const_mav<double,2> uvw2(uvw, {nrows, 3});
  gridder::const_mav<double,1> freq2(freq, {nchan});
  gridder::mav<complex<float>,2> ms(vis, {nrows,nchan});
  gridder::const_mav<float,2> tdirty(img.data(), {width_, height_});
  gridder::const_mav<float,2> twgt(nullptr, {0,0});
  gridder::dirty2ms(uvw2, freq2, tdirty, twgt,
    pixelSizeX_, pixelSizeY_, epsilon_, true, nthreads_, ms, 0);
  }
