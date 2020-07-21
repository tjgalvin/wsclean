#ifndef GRIDDER_CXX_H
#define GRIDDER_CXX_H

/*
 *  This file is part of nifty_gridder.
 *
 *  nifty_gridder is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  nifty_gridder is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with nifty_gridder; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/* Copyright (C) 2019 Max-Planck-Society
   Author: Martin Reinecke */

#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <array>
#include <fftw3.h>

#include "threading.h"

#if defined(__GNUC__)
#define NOINLINE __attribute__((noinline))
#define ALIGNED(align) __attribute__((aligned(align)))
#define RESTRICT __restrict__
#else
#define NOINLINE
#define ALIGNED(align)
#define RESTRICT
#endif

namespace gridder {

namespace detail {

using namespace std;

template <typename T>
struct VLEN {
  static constexpr size_t val = 1;
};

#if (defined(__AVX512F__))
template <>
struct VLEN<float> {
  static constexpr size_t val = 16;
};
template <>
struct VLEN<double> {
  static constexpr size_t val = 8;
};
#elif (defined(__AVX__))
template <>
struct VLEN<float> {
  static constexpr size_t val = 8;
};
template <>
struct VLEN<double> {
  static constexpr size_t val = 4;
};
#elif (defined(__SSE2__))
template <>
struct VLEN<float> {
  static constexpr size_t val = 4;
};
template <>
struct VLEN<double> {
  static constexpr size_t val = 2;
};
#elif (defined(__VSX__))
template <>
struct VLEN<float> {
  static constexpr size_t val = 4;
};
template <>
struct VLEN<double> {
  static constexpr size_t val = 2;
};
#endif

// "mav" stands for "multidimensional array view"
template <typename T, size_t ndim>
class mav {
  static_assert((ndim > 0) && (ndim < 3), "only supports 1D and 2D arrays");

 private:
  T *d;
  array<size_t, ndim> shp;
  array<ptrdiff_t, ndim> str;

 public:
  mav(T *d_, const array<size_t, ndim> &shp_,
      const array<ptrdiff_t, ndim> &str_)
      : d(d_), shp(shp_), str(str_) {}
  mav(T *d_, const array<size_t, ndim> &shp_) : d(d_), shp(shp_) {
    str[ndim - 1] = 1;
    for (size_t d = 2; d <= ndim; ++d)
      str[ndim - d] = str[ndim - d + 1] * shp[ndim - d + 1];
  }
  T &operator[](size_t i) const { return operator()(i); }
  T &operator()(size_t i) const {
    static_assert(ndim == 1, "ndim must be 1");
    return d[str[0] * i];
  }
  T &operator()(size_t i, size_t j) const {
    static_assert(ndim == 2, "ndim must be 2");
    return d[str[0] * i + str[1] * j];
  }
  size_t shape(size_t i) const { return shp[i]; }
  const array<size_t, ndim> &shape() const { return shp; }
  size_t size() const {
    size_t res = 1;
    for (auto v : shp) res *= v;
    return res;
  }
  ptrdiff_t stride(size_t i) const { return str[i]; }
  T *data() const { return d; }
  bool last_contiguous() const {
    return (str[ndim - 1] == 1) || (str[ndim - 1] == 0);
  }
  void check_storage(const char *name) const {
    if (!last_contiguous())
      cout << "Array '" << name
           << "': last dimension is not contiguous.\n"
              "This may slow down computation significantly!\n";
  }
  bool contiguous() const {
    ptrdiff_t stride = 1;
    for (size_t i = 0; i < ndim; ++i) {
      if (str[ndim - 1 - i] != stride) return false;
      stride *= shp[ndim - 1 - i];
    }
    return true;
  }
  void fill(const T &val) const {
    // FIXME: special cases for contiguous arrays and/or zeroing?
    if (ndim == 1)
      for (size_t i = 0; i < shp[0]; ++i) d[str[0] * i] = val;
    else if (ndim == 2)
      for (size_t i = 0; i < shp[0]; ++i)
        for (size_t j = 0; j < shp[1]; ++j) d[str[0] * i + str[1] * j] = val;
  }
};

template <typename T, size_t ndim>
using const_mav = mav<const T, ndim>;
template <typename T, size_t ndim>
const_mav<T, ndim> cmav(const mav<T, ndim> &mav) {
  return const_mav<T, ndim>(mav.data(), mav.shape());
}
template <typename T, size_t ndim>
const_mav<T, ndim> nullmav() {
  array<size_t, ndim> shp;
  shp.fill(0);
  return const_mav<T, ndim>(nullptr, shp);
}

//
// basic utilities
//
#if defined(__GNUC__)
#define LOC_ CodeLocation(__FILE__, __LINE__, __PRETTY_FUNCTION__)
#else
#define LOC_ CodeLocation(__FILE__, __LINE__)
#endif

#define myfail(...)                                    \
  do {                                                 \
    std::ostringstream os;                             \
    streamDump__(os, LOC_, "\n", ##__VA_ARGS__, "\n"); \
    throw std::runtime_error(os.str());                \
  } while (0)

#define myassert(cond, ...)                         \
  do {                                              \
    if (cond)                                       \
      ;                                             \
    else {                                          \
      myfail("Assertion failure\n", ##__VA_ARGS__); \
    }                                               \
  } while (0)

template <typename T>
inline void streamDump__(std::ostream &os, const T &value) {
  os << value;
}

template <typename T, typename... Args>
inline void streamDump__(std::ostream &os, const T &value,
                         const Args &... args) {
  os << value;
  streamDump__(os, args...);
}
// to be replaced with std::source_location once available
class CodeLocation {
 private:
  const char *file, *func;
  int line;

 public:
  CodeLocation(const char *file_, int line_, const char *func_ = nullptr)
      : file(file_), func(func_), line(line_) {}

  ostream &print(ostream &os) const {
    os << "file: " << file << ", line: " << line;
    if (func) os << ", function: " << func;
    return os;
  }
};

inline std::ostream &operator<<(std::ostream &os, const CodeLocation &loc) {
  return loc.print(os);
}

template <size_t ndim>
void checkShape(const array<size_t, ndim> &shp1,
                const array<size_t, ndim> &shp2) {
  for (size_t i = 0; i < ndim; ++i)
    myassert(shp1[i] == shp2[i], "shape mismatch");
}

template <typename T>
inline T fmod1(T v) {
  return v - floor(v);
}

template <typename T>
class aligned_data {
 private:
  T *d;
  size_t s;

 public:
  aligned_data(size_t s_)
      : d(reinterpret_cast<T *>(fftw_malloc(s_ * sizeof(T)))), s(s_) {
    myassert(d != nullptr, "allocation failed");
  }
  ~aligned_data() { fftw_free(d); }
  T *data() { return d; }
  const T *data() const { return d; }
  size_t size() const { return s; }
  T *begin() { return d; }
  T *end() { return d + s; }
};

template <typename T, size_t ndim>
class tmpStorage {
 private:
  aligned_data<T> d;
  mav<T, ndim> mav_;

  static size_t prod(const array<size_t, ndim> &shp) {
    size_t res = 1;
    for (auto v : shp) res *= v;
    return res;
  }

 public:
  tmpStorage(const array<size_t, ndim> &shp)
      : d(prod(shp)), mav_(d.data(), shp) {}
  mav<T, ndim> &getMav() { return mav_; }
  const_mav<T, ndim> getCmav() { return cmav(mav_); }
  void fill(const T &val) { std::fill(d.begin(), d.end(), val); }
};

template <typename T>
void prep_fft(size_t s0, size_t s1, size_t nthreads, size_t verbosity) {
  myfail("not implemented");
}
template <>
void prep_fft<float>(size_t s0, size_t s1, size_t nthreads, size_t verbosity) {
  if (verbosity > 0) cout << "preparing FFT plans ... " << flush;
  aligned_data<fftwf_complex> tmp(s0 * s1);
  fftwf_plan_with_nthreads(nthreads);
  auto fplan = fftwf_plan_dft_2d(s0, s1, tmp.data(), tmp.data(), FFTW_FORWARD,
                                 FFTW_MEASURE);
  fftwf_destroy_plan(fplan);
  auto bplan = fftwf_plan_dft_2d(s0, s1, tmp.data(), tmp.data(), FFTW_BACKWARD,
                                 FFTW_MEASURE);
  fftwf_destroy_plan(bplan);
  if (verbosity > 0) cout << "done" << endl;
}
template <typename T>
void exec_fft(const mav<complex<T>, 2> &arr, bool fwd, size_t nthreads) {
  myfail("not implemented");
}
template <>
void exec_fft(const mav<complex<float>, 2> &arr, bool fwd, size_t nthreads) {
  auto d = reinterpret_cast<fftwf_complex *>(arr.data());
  fftwf_plan_with_nthreads(nthreads);
  auto plan =
      fftwf_plan_dft_2d(arr.shape(0), arr.shape(1), d, d,
                        fwd ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
}

//
// Utilities for Gauss-Legendre quadrature
//

static inline double one_minus_x2(double x) {
  return (fabs(x) > 0.1) ? (1. + x) * (1. - x) : 1. - x * x;
}

void legendre_prep(int n, vector<double> &x, vector<double> &w,
                   size_t nthreads) {
  constexpr double pi = 3.141592653589793238462643383279502884197;
  constexpr double eps = 3e-14;
  int m = (n + 1) >> 1;
  x.resize(m);
  w.resize(m);

  const double t0 = 1 - (1 - 1. / n) / (8. * n * n);
  const double t1 = 1. / (4. * n + 2.);

  threading::Scheduler<size_t> sched;
  sched.execDynamic(1, m + 1, nthreads, 1, 0., [&]() {
    size_t lo, hi;
    while (sched.getNext(lo, hi))
      for (size_t i = lo; i < hi; ++i) {
        double x0 = cos(pi * ((i << 2) - 1) * t1) * t0;

        int dobreak = 0;
        int j = 0;
        double dpdx;
        while (1) {
          double P_1 = 1.0;
          double P0 = x0;
          double dx, x1;

          for (int k = 2; k <= n; k++) {
            double P_2 = P_1;
            P_1 = P0;
            //            P0 = ((2*k-1)*x0*P_1-(k-1)*P_2)/k;
            P0 = x0 * P_1 + (k - 1.) / k * (x0 * P_1 - P_2);
          }

          dpdx = (P_1 - x0 * P0) * n / one_minus_x2(x0);

          /* Newton step */
          x1 = x0 - P0 / dpdx;
          dx = x0 - x1;
          x0 = x1;
          if (dobreak) break;

          if (abs(dx) <= eps) dobreak = 1;
          myassert(++j < 100, "convergence problem");
        }

        x[m - i] = x0;
        w[m - i] = 2. / (one_minus_x2(x0) * dpdx * dpdx);
      }
  });
}

//
// Start of real gridder functionality
//

class ES_Kernel {
 private:
  static constexpr double pi = 3.141592653589793238462643383279502884197;
  double beta;
  int p;
  vector<double> x, wgt, psi;
  size_t supp;

 public:
  ES_Kernel(size_t supp_, double ofactor, size_t nthreads)
      : beta(get_beta(supp_, ofactor) * supp_),
        p(int(1.5 * supp_ + 2)),
        supp(supp_) {
    legendre_prep(2 * p, x, wgt, nthreads);
    psi = x;
    for (auto &v : psi) v = operator()(v);
  }
  ES_Kernel(size_t supp_, size_t nthreads) : ES_Kernel(supp_, 2., nthreads) {}

  double operator()(double v) const {
    return exp(beta * (sqrt(1. - v * v) - 1.));
  }
  /* Compute correction factors for the ES gridding kernel
     This implementation follows eqs. (3.8) to (3.10) of Barnett et al. 2018 */
  double corfac(double v) const {
    double tmp = 0;
    for (int i = 0; i < p; ++i)
      tmp += wgt[i] * psi[i] * cos(pi * supp * v * x[i]);
    return 1. / (supp * tmp);
  }
  static double get_beta(size_t supp, double ofactor = 2) {
    myassert((supp >= 2) && (supp <= 15), "unsupported support size");
    if (ofactor >= 2) {
      static const vector<double> opt_beta{
          -1,    0.14,   1.70,   2.08,  2.205, 2.26,  2.29,   2.307,
          2.316, 2.3265, 2.3324, 2.282, 2.294, 2.304, 2.3138, 2.317};
      myassert(supp < opt_beta.size(), "bad support size");
      return opt_beta[supp];
    }
    if (ofactor >= 1.15) {
      // empirical, but pretty accurate approximation
      static const array<double, 16> betacorr{
          0, 0, -0.51, -0.21, -0.1, -0.05, -0.025, -0.0125,
          0, 0, 0,     0,     0,    0,     0,      0};
      auto x0 = 1. / (2 * ofactor);
      auto bcstrength = 1. + (x0 - 0.25) * 2.5;
      return 2.32 + bcstrength * betacorr[supp] + (0.25 - x0) * 3.1;
    }
    myfail("oversampling factor ", ofactor, " is too small");
  }

  static size_t get_supp(double epsilon, double ofactor = 2) {
    double epssq = epsilon * epsilon;
    if (ofactor >= 2) {
      static const vector<double> maxmaperr{
          1e8,      0.19,     2.98e-3,  5.98e-5,  1.11e-6,  2.01e-8,
          3.55e-10, 5.31e-12, 8.81e-14, 1.34e-15, 2.17e-17, 2.12e-19,
          2.88e-21, 3.92e-23, 8.21e-25, 7.13e-27};

      for (size_t i = 1; i < maxmaperr.size(); ++i)
        if (epssq > maxmaperr[i]) return i;
      myfail("requested epsilon too small - minimum is 1e-13");
    }
    if (ofactor >= 1.15) {
      for (size_t w = 2; w < 16; ++w) {
        auto estimate =
            12 *
            exp(-2. * w * ofactor);  // empirical, not very good approximation
        if (epssq > estimate) return w;
      }
      myfail("requested epsilon too small");
    }
    myfail("oversampling factor ", ofactor, " is too small");
  }
};

/* Compute correction factors for the ES gridding kernel
   This implementation follows eqs. (3.8) to (3.10) of Barnett et al. 2018 */
vector<double> correction_factors(size_t n, double ofactor, size_t nval,
                                  size_t supp, size_t nthreads) {
  ES_Kernel kernel(supp, ofactor, nthreads);
  vector<double> res(nval);
  double xn = 1. / n;
  threading::Scheduler<size_t> sched;
  sched.execStatic(0, nval, nthreads, 0, [&]() {
    size_t lo, hi;
    while (sched.getNext(lo, hi))
      for (size_t k = lo; k != hi; ++k) res[k] = kernel.corfac(k * xn);
  });
  return res;
}

using idx_t = uint32_t;

struct RowChan {
  idx_t row, chan;
};

struct UVW {
  double u, v, w;
  UVW() {}
  UVW(double u_, double v_, double w_) : u(u_), v(v_), w(w_) {}
  UVW operator*(double fct) const { return UVW(u * fct, v * fct, w * fct); }
  void Flip() {
    u = -u;
    v = -v;
    w = -w;
  }
  bool FixW() {
    bool flip = w < 0;
    if (flip) Flip();
    return flip;
  }
};

class Baselines {
 protected:
  vector<UVW> coord;
  vector<double> f_over_c;
  idx_t nrows, nchan;
  idx_t shift, mask;

 public:
  template <typename T>
  Baselines(const const_mav<T, 2> &coord_, const const_mav<T, 1> &freq,
            bool negate_v = false) {
    constexpr double speedOfLight = 299792458.;
    myassert(coord_.shape(1) == 3, "dimension mismatch");
    auto hugeval = size_t(~(idx_t(0)));
    myassert(coord_.shape(0) < hugeval, "too many entries in MS");
    myassert(coord_.shape(1) < hugeval, "too many entries in MS");
    myassert(coord_.size() < hugeval, "too many entries in MS");
    nrows = coord_.shape(0);
    nchan = freq.shape(0);
    shift = 0;
    while ((idx_t(1) << shift) < nchan) ++shift;
    mask = (idx_t(1) << shift) - 1;
    myassert(nrows * (mask + 1) < hugeval, "too many entries in MS");
    f_over_c.resize(nchan);
    for (size_t i = 0; i < nchan; ++i) {
      myassert(freq[i] > 0, "negative channel frequency encountered");
      f_over_c[i] = freq(i) / speedOfLight;
    }
    coord.resize(nrows);
    if (negate_v)
      for (size_t i = 0; i < coord.size(); ++i)
        coord[i] = UVW(coord_(i, 0), -coord_(i, 1), coord_(i, 2));
    else
      for (size_t i = 0; i < coord.size(); ++i)
        coord[i] = UVW(coord_(i, 0), coord_(i, 1), coord_(i, 2));
  }

  RowChan getRowChan(idx_t index) const {
    return RowChan{index >> shift, index & mask};
  }

  UVW effectiveCoord(const RowChan &rc) const {
    return coord[rc.row] * f_over_c[rc.chan];
  }
  UVW effectiveCoord(idx_t index) const {
    return effectiveCoord(getRowChan(index));
  }
  size_t Nrows() const { return nrows; }
  size_t Nchannels() const { return nchan; }
  idx_t getIdx(idx_t irow, idx_t ichan) const {
    return ichan + (irow << shift);
  }
};

class GridderConfig {
 protected:
  size_t nx_dirty, ny_dirty, nu, nv;
  double ofactor, eps, psx, psy;
  size_t supp, nsafe;
  double beta;
  size_t nthreads;
  double ushift, vshift;
  int maxiu0, maxiv0;

  template <typename T>
  complex<T> wscreen(T x, T y, T w, bool adjoint) const {
    constexpr T pi = T(3.141592653589793238462643383279502884197);
    T tmp = 1 - x - y;
    if (tmp <= 0) return 1;              // no phase factor beyond the horizon
    T nm1 = (-x - y) / (sqrt(tmp) + 1);  // more accurate form of sqrt(1-x-y)-1
    T phase = 2 * pi * w * nm1;
    if (adjoint) phase *= -1;
    return complex<T>(cos(phase), sin(phase));
  }

 public:
  GridderConfig(size_t nxdirty, size_t nydirty, size_t nu_, size_t nv_,
                double epsilon, double pixsize_x, double pixsize_y,
                size_t nthreads_)
      : nx_dirty(nxdirty),
        ny_dirty(nydirty),
        nu(nu_),
        nv(nv_),
        ofactor(min(double(nu) / nxdirty, double(nv) / nydirty)),
        eps(epsilon),
        psx(pixsize_x),
        psy(pixsize_y),
        supp(ES_Kernel::get_supp(epsilon, ofactor)),
        nsafe((supp + 1) / 2),
        beta(ES_Kernel::get_beta(supp, ofactor) * supp),
        nthreads(nthreads_),
        ushift(supp * (-0.5) + 1 + nu),
        vshift(supp * (-0.5) + 1 + nv),
        maxiu0((nu + nsafe) - supp),
        maxiv0((nv + nsafe) - supp) {
    myassert(nu >= 2 * nsafe, "nu too small");
    myassert(nv >= 2 * nsafe, "nv too small");
    myassert((nx_dirty & 1) == 0, "nx_dirty must be even");
    myassert((ny_dirty & 1) == 0, "ny_dirty must be even");
    myassert((nu & 1) == 0, "nu must be even");
    myassert((nv & 1) == 0, "nv must be even");
    myassert(epsilon > 0, "epsilon must be positive");
    myassert(pixsize_x > 0, "pixsize_x must be positive");
    myassert(pixsize_y > 0, "pixsize_y must be positive");
    myassert(ofactor >= 1.15, "oversampling factor ", ofactor,
             " smaller than 1.15");
  }
  GridderConfig(size_t nxdirty, size_t nydirty, double epsilon,
                double pixsize_x, double pixsize_y, size_t nthreads_)
      : GridderConfig(nxdirty, nydirty, max<size_t>(30, 2 * nxdirty),
                      max<size_t>(30, 2 * nydirty), epsilon, pixsize_x,
                      pixsize_y, nthreads_) {}
  size_t Nxdirty() const { return nx_dirty; }
  size_t Nydirty() const { return ny_dirty; }
  double Epsilon() const { return eps; }
  double Pixsize_x() const { return psx; }
  double Pixsize_y() const { return psy; }
  size_t Nu() const { return nu; }
  size_t Nv() const { return nv; }
  size_t Supp() const { return supp; }
  size_t Nsafe() const { return nsafe; }
  double Beta() const { return beta; }
  size_t Nthreads() const { return nthreads; }
  double Ofactor() const { return ofactor; }

  template <typename T>
  void grid2dirty_post2(const mav<complex<T>, 2> &tmav, const mav<T, 2> &dirty,
                        T w) const {
    checkShape(dirty.shape(), {nx_dirty, ny_dirty});
    double x0 = -0.5 * nx_dirty * psx, y0 = -0.5 * ny_dirty * psy;
    threading::Scheduler<size_t> sched;
    sched.execStatic(0, nx_dirty / 2 + 1, nthreads, 0, [&]() {
      size_t lo, hi;
      while (sched.getNext(lo, hi))
        for (size_t i = lo; i != hi; ++i) {
          T fx = T(x0 + i * psx);
          fx *= fx;
          for (size_t j = 0; j <= ny_dirty / 2; ++j) {
            T fy = T(y0 + j * psy);
            auto ws = wscreen(fx, fy * fy, w, true);
            size_t ix = nu - nx_dirty / 2 + i;
            if (ix >= nu) ix -= nu;
            size_t jx = nv - ny_dirty / 2 + j;
            if (jx >= nv) jx -= nv;
            dirty(i, j) += (tmav(ix, jx) * ws).real();  // lower left
            size_t i2 = nx_dirty - i, j2 = ny_dirty - j;
            size_t ix2 = nu - nx_dirty / 2 + i2;
            if (ix2 >= nu) ix2 -= nu;
            size_t jx2 = nv - ny_dirty / 2 + j2;
            if (jx2 >= nv) jx2 -= nv;
            if ((i > 0) && (i < i2)) {
              dirty(i2, j) += (tmav(ix2, jx) * ws).real();  // lower right
              if ((j > 0) && (j < j2))
                dirty(i2, j2) += (tmav(ix2, jx2) * ws).real();  // upper right
            }
            if ((j > 0) && (j < j2))
              dirty(i, j2) += (tmav(ix, jx2) * ws).real();  // upper left
          }
        }
    });
  }

  template <typename T>
  void grid2dirty_c_overwrite_wscreen_add(const mav<complex<T>, 2> &grid,
                                          const mav<T, 2> &dirty, T w) const {
    checkShape(grid.shape(), {nu, nv});
    exec_fft(grid, false, nthreads);
    grid2dirty_post2(grid, dirty, w);
  }

  template <typename T>
  void dirty2grid_pre2(const const_mav<T, 2> &dirty,
                       const mav<complex<T>, 2> &grid, T w) const {
    checkShape(dirty.shape(), {nx_dirty, ny_dirty});
    checkShape(grid.shape(), {nu, nv});
    grid.fill(0);

    double x0 = -0.5 * nx_dirty * psx, y0 = -0.5 * ny_dirty * psy;
    threading::Scheduler<size_t> sched;
    sched.execStatic(0, nx_dirty / 2 + 1, nthreads, 0, [&]() {
      size_t lo, hi;
      while (sched.getNext(lo, hi))
        for (size_t i = lo; i != hi; ++i) {
          T fx = T(x0 + i * psx);
          fx *= fx;
          for (size_t j = 0; j <= ny_dirty / 2; ++j) {
            T fy = T(y0 + j * psy);
            auto ws = wscreen(fx, fy * fy, w, false);
            size_t ix = nu - nx_dirty / 2 + i;
            if (ix >= nu) ix -= nu;
            size_t jx = nv - ny_dirty / 2 + j;
            if (jx >= nv) jx -= nv;
            grid(ix, jx) = dirty(i, j) * ws;  // lower left
            size_t i2 = nx_dirty - i, j2 = ny_dirty - j;
            size_t ix2 = nu - nx_dirty / 2 + i2;
            if (ix2 >= nu) ix2 -= nu;
            size_t jx2 = nv - ny_dirty / 2 + j2;
            if (jx2 >= nv) jx2 -= nv;
            if ((i > 0) && (i < i2)) {
              grid(ix2, jx) = dirty(i2, j) * ws;  // lower right
              if ((j > 0) && (j < j2))
                grid(ix2, jx2) = dirty(i2, j2) * ws;  // upper right
            }
            if ((j > 0) && (j < j2))
              grid(ix, jx2) = dirty(i, j2) * ws;  // upper left
          }
        }
    });
  }

  template <typename T>
  void dirty2grid_c_wscreen(const const_mav<T, 2> &dirty,
                            const mav<complex<T>, 2> &grid, T w) const {
    dirty2grid_pre2(dirty, grid, w);
    exec_fft(grid, true, nthreads);
  }

  void getpix(double u_in, double v_in, double &u, double &v, int &iu0,
              int &iv0) const {
    u = fmod1(u_in * psx) * nu;
    iu0 = min(int(u + ushift) - int(nu), maxiu0);
    v = fmod1(v_in * psy) * nv;
    iv0 = min(int(v + vshift) - int(nv), maxiv0);
  }

  template <typename T>
  void apply_wscreen(const const_mav<complex<T>, 2> &dirty,
                     mav<complex<T>, 2> &dirty2, double w, bool adjoint) const {
    checkShape(dirty.shape(), {nx_dirty, ny_dirty});
    checkShape(dirty2.shape(), {nx_dirty, ny_dirty});

    double x0 = -0.5 * nx_dirty * psx, y0 = -0.5 * ny_dirty * psy;
    threading::Scheduler<size_t> sched;
    sched.execStatic(0, nx_dirty / 2 + 1, nthreads, 0, [&]() {
      size_t lo, hi;
      while (sched.getNext(lo, hi))
        for (size_t i = lo; i != hi; ++i) {
          T fx = T(x0 + i * psx);
          fx *= fx;
          for (size_t j = 0; j <= ny_dirty / 2; ++j) {
            T fy = T(y0 + j * psy);
            auto ws = wscreen(fx, fy * fy, T(w), adjoint);
            dirty2(i, j) = dirty(i, j) * ws;  // lower left
            size_t i2 = nx_dirty - i, j2 = ny_dirty - j;
            if ((i > 0) && (i < i2)) {
              dirty2(i2, j) = dirty(i2, j) * ws;  // lower right
              if ((j > 0) && (j < j2))
                dirty2(i2, j2) = dirty(i2, j2) * ws;  // upper right
            }
            if ((j > 0) && (j < j2))
              dirty2(i, j2) = dirty(i, j2) * ws;  // upper left
          }
        }
    });
  }
};

constexpr int logsquare = 4;

template <typename T, typename T2 = complex<T>>
class Helper {
 private:
  const GridderConfig &gconf;
  int nu, nv, nsafe, supp;
  T beta;
  const T2 *grid_r;
  T2 *grid_w;
  int su, sv;
  int iu0, iv0;  // start index of the current visibility
  int bu0, bv0;  // start index of the current buffer

  vector<T2> rbuf, wbuf;
  bool do_w_gridding;
  double w0, xdw;
  size_t nexp;
  size_t nvecs;
  vector<std::mutex> &locks;

  void dump() const {
    if (bu0 < -nsafe) return;  // nothing written into buffer yet

    int idxu = (bu0 + nu) % nu;
    int idxv0 = (bv0 + nv) % nv;
    for (int iu = 0; iu < su; ++iu) {
      int idxv = idxv0;
      {
        std::lock_guard<std::mutex> lock(locks[idxu]);
        for (int iv = 0; iv < sv; ++iv) {
          grid_w[idxu * nv + idxv] += wbuf[iu * sv + iv];
          if (++idxv >= nv) idxv = 0;
        }
      }
      if (++idxu >= nu) idxu = 0;
    }
  }

  void load() {
    int idxu = (bu0 + nu) % nu;
    int idxv0 = (bv0 + nv) % nv;
    for (int iu = 0; iu < su; ++iu) {
      int idxv = idxv0;
      for (int iv = 0; iv < sv; ++iv) {
        rbuf[iu * sv + iv] = grid_r[idxu * nv + idxv];
        if (++idxv >= nv) idxv = 0;
      }
      if (++idxu >= nu) idxu = 0;
    }
  }

 public:
  const T2 *p0r;
  T2 *p0w;
  T kernel[64] ALIGNED(64);

  Helper(const GridderConfig &gconf_, const T2 *grid_r_, T2 *grid_w_,
         vector<std::mutex> &locks_, double w0_ = -1, double dw_ = -1)
      : gconf(gconf_),
        nu(gconf.Nu()),
        nv(gconf.Nv()),
        nsafe(gconf.Nsafe()),
        supp(gconf.Supp()),
        beta(T(gconf.Beta())),
        grid_r(grid_r_),
        grid_w(grid_w_),
        su(2 * nsafe + (1 << logsquare)),
        sv(2 * nsafe + (1 << logsquare)),
        bu0(-1000000),
        bv0(-1000000),
        rbuf(su * sv * (grid_r != nullptr), T(0)),
        wbuf(su * sv * (grid_w != nullptr), T(0)),
        do_w_gridding(dw_ > 0),
        w0(w0_),
        xdw(T(1) / dw_),
        nexp(2 * supp + do_w_gridding),
        nvecs(VLEN<T>::val * ((nexp + VLEN<T>::val - 1) / VLEN<T>::val)),
        locks(locks_) {}
  ~Helper() {
    if (grid_w) dump();
  }

  int lineJump() const { return sv; }
  T Wfac() const { return kernel[2 * supp]; }
  void prep(const UVW &in) {
    double u, v;
    gconf.getpix(in.u, in.v, u, v, iu0, iv0);
    double xsupp = 2. / supp;
    double x0 = xsupp * (iu0 - u);
    double y0 = xsupp * (iv0 - v);
    for (int i = 0; i < supp; ++i) {
      kernel[i] = T(x0 + i * xsupp);
      kernel[i + supp] = T(y0 + i * xsupp);
    }
    if (do_w_gridding)
      kernel[2 * supp] = min(T(1), T(xdw * xsupp * abs(w0 - in.w)));
    for (size_t i = nexp; i < nvecs; ++i) kernel[i] = 0;
    for (size_t i = 0; i < nvecs; ++i)
      kernel[i] = exp(beta * (sqrt(T(1) - kernel[i] * kernel[i]) - T(1)));
    if ((iu0 < bu0) || (iv0 < bv0) || (iu0 + supp > bu0 + su) ||
        (iv0 + supp > bv0 + sv)) {
      if (grid_w) {
        dump();
        fill(wbuf.begin(), wbuf.end(), T(0));
      }
      bu0 = ((((iu0 + nsafe) >> logsquare) << logsquare)) - nsafe;
      bv0 = ((((iv0 + nsafe) >> logsquare) << logsquare)) - nsafe;
      if (grid_r) load();
    }
    p0r = grid_r ? rbuf.data() + sv * (iu0 - bu0) + iv0 - bv0 : nullptr;
    p0w = grid_w ? wbuf.data() + sv * (iu0 - bu0) + iv0 - bv0 : nullptr;
  }
};

template <class T, class Serv>
class SubServ {
 private:
  const Serv &srv;
  const_mav<idx_t, 1> subidx;

 public:
  SubServ(const Serv &orig, const const_mav<idx_t, 1> &subidx_)
      : srv(orig), subidx(subidx_) {}
  size_t Nvis() const { return subidx.size(); }
  const Baselines &getBaselines() const { return srv.getBaselines(); }
  UVW getCoord(size_t i) const { return srv.getCoord(subidx[i]); }
  complex<T> getVis(size_t i) const { return srv.getVis(subidx[i]); }
  idx_t getIdx(size_t i) const { return srv.getIdx(subidx[i]); }
  void setVis(size_t i, const complex<T> &v) const { srv.setVis(subidx[i], v); }
  void addVis(size_t i, const complex<T> &v) const { srv.addVis(subidx[i], v); }
};

template <class T, class T2>
class MsServ {
 private:
  const Baselines &baselines;
  const_mav<idx_t, 1> idx;
  T2 ms;
  const_mav<T, 2> wgt;
  size_t nvis;
  bool have_wgt;

 public:
  using Tsub = SubServ<T, MsServ>;

  MsServ(const Baselines &baselines_, const const_mav<idx_t, 1> &idx_, T2 ms_,
         const const_mav<T, 2> &wgt_)
      : baselines(baselines_),
        idx(idx_),
        ms(ms_),
        wgt(wgt_),
        nvis(idx.shape(0)),
        have_wgt(wgt.size() != 0) {
    wgt.check_storage("wgt");
    ms.check_storage("ms");
    checkShape(ms.shape(), {baselines.Nrows(), baselines.Nchannels()});
    if (have_wgt) checkShape(wgt.shape(), ms.shape());
  }
  Tsub getSubserv(const const_mav<idx_t, 1> &subidx) const {
    return Tsub(*this, subidx);
  }
  size_t Nvis() const { return nvis; }
  const Baselines &getBaselines() const { return baselines; }
  UVW getCoord(size_t i) const { return baselines.effectiveCoord(idx(i)); }
  complex<T> getVis(size_t i) const {
    auto rc = baselines.getRowChan(idx(i));
    return have_wgt ? ms(rc.row, rc.chan) * wgt(rc.row, rc.chan)
                    : ms(rc.row, rc.chan);
  }
  idx_t getIdx(size_t i) const { return idx[i]; }
  void setVis(size_t i, const complex<T> &v) const {
    auto rc = baselines.getRowChan(idx(i));
    ms(rc.row, rc.chan) = have_wgt ? v * wgt(rc.row, rc.chan) : v;
  }
  void addVis(size_t i, const complex<T> &v) const {
    auto rc = baselines.getRowChan(idx(i));
    ms(rc.row, rc.chan) += have_wgt ? v * wgt(rc.row, rc.chan) : v;
  }
};
template <class T, class T2>
MsServ<T, T2> makeMsServ(const Baselines &baselines,
                         const const_mav<idx_t, 1> &idx, const T2 &ms,
                         const const_mav<T, 2> &wgt) {
  return MsServ<T, T2>(baselines, idx, ms, wgt);
}

template <typename T, typename Serv>
void x2grid_c(const GridderConfig &gconf, const Serv &srv,
              mav<complex<T>, 2> &grid, double w0 = -1, double dw = -1) {
  checkShape(grid.shape(), {gconf.Nu(), gconf.Nv()});
  myassert(grid.contiguous(), "grid is not contiguous");
  size_t supp = gconf.Supp();
  size_t nthreads = gconf.Nthreads();
  bool do_w_gridding = dw > 0;
  vector<std::mutex> locks(gconf.Nu());

  size_t np = srv.Nvis();
  threading::Scheduler<size_t> sched;
  //  sched.execStatic(0, np, nthreads, 0, [&]()
  sched.execDynamic(0, np, nthreads, 100, 0.2, [&]() {
    Helper<T> hlp(gconf, nullptr, grid.data(), locks, w0, dw);
    int jump = hlp.lineJump();
    const T *RESTRICT ku = hlp.kernel;
    const T *RESTRICT kv = hlp.kernel + supp;

    size_t lo, hi;
    while (sched.getNext(lo, hi))
      for (size_t ipart = lo; ipart != hi; ++ipart) {
        UVW coord = srv.getCoord(ipart);
        auto flip = coord.FixW();
        hlp.prep(coord);
        auto *RESTRICT ptr = hlp.p0w;
        auto v(srv.getVis(ipart));
        if (do_w_gridding) v *= hlp.Wfac();
        if (flip) v = conj(v);
        for (size_t cu = 0; cu < supp; ++cu) {
          complex<T> tmp(v * ku[cu]);
          size_t cv = 0;
          for (; cv + 3 < supp; cv += 4) {
            ptr[cv] += tmp * kv[cv];
            ptr[cv + 1] += tmp * kv[cv + 1];
            ptr[cv + 2] += tmp * kv[cv + 2];
            ptr[cv + 3] += tmp * kv[cv + 3];
          }
          for (; cv < supp; ++cv) ptr[cv] += tmp * kv[cv];
          ptr += jump;
        }
      }
  });
}

template <typename T, typename Serv>
void grid2x_c(const GridderConfig &gconf, const const_mav<complex<T>, 2> &grid,
              const Serv &srv, double w0 = -1, double dw = -1) {
  checkShape(grid.shape(), {gconf.Nu(), gconf.Nv()});
  myassert(grid.contiguous(), "grid is not contiguous");
  size_t supp = gconf.Supp();
  size_t nthreads = gconf.Nthreads();
  bool do_w_gridding = dw > 0;
  vector<std::mutex> locks(gconf.Nu());

  // Loop over sampling points
  size_t np = srv.Nvis();
  threading::Scheduler<size_t> sched;
  sched.execDynamic(0, np, nthreads, 1000, 0.5, [&]() {
    Helper<T> hlp(gconf, grid.data(), nullptr, locks, w0, dw);
    int jump = hlp.lineJump();
    const T *RESTRICT ku = hlp.kernel;
    const T *RESTRICT kv = hlp.kernel + supp;

    size_t lo, hi;
    while (sched.getNext(lo, hi))
      for (size_t ipart = lo; ipart != hi; ++ipart) {
        UVW coord = srv.getCoord(ipart);
        auto flip = coord.FixW();
        hlp.prep(coord);
        complex<T> r = 0;
        const auto *RESTRICT ptr = hlp.p0r;
        for (size_t cu = 0; cu < supp; ++cu) {
          complex<T> tmp(0);
          size_t cv = 0;
          for (; cv + 3 < supp; cv += 4)
            tmp += ptr[cv] * kv[cv] + ptr[cv + 1] * kv[cv + 1] +
                   ptr[cv + 2] * kv[cv + 2] + ptr[cv + 3] * kv[cv + 3];
          for (; cv < supp; ++cv) tmp += ptr[cv] * kv[cv];
          r += tmp * ku[cu];
          ptr += jump;
        }
        if (flip) r = conj(r);
        if (do_w_gridding) r *= hlp.Wfac();
        srv.addVis(ipart, r);
      }
  });
}

template <typename T>
void apply_global_corrections(const GridderConfig &gconf,
                              const mav<T, 2> &dirty, const ES_Kernel &kernel,
                              double dw, bool divide_by_n) {
  auto nx_dirty = gconf.Nxdirty();
  auto ny_dirty = gconf.Nydirty();
  size_t nthreads = gconf.Nthreads();
  auto psx = gconf.Pixsize_x();
  auto psy = gconf.Pixsize_y();
  double x0 = -0.5 * nx_dirty * psx, y0 = -0.5 * ny_dirty * psy;
  auto cfu = correction_factors(gconf.Nu(), gconf.Ofactor(), nx_dirty / 2 + 1,
                                gconf.Supp(), nthreads);
  auto cfv = correction_factors(gconf.Nv(), gconf.Ofactor(), ny_dirty / 2 + 1,
                                gconf.Supp(), nthreads);
  threading::Scheduler<size_t> sched;
  sched.execStatic(0, nx_dirty / 2 + 1, nthreads, 0, [&]() {
    size_t lo, hi;
    while (sched.getNext(lo, hi))
      for (size_t i = lo; i != hi; ++i) {
        auto fx = T(x0 + i * psx);
        fx *= fx;
        for (size_t j = 0; j <= ny_dirty / 2; ++j) {
          auto fy = T(y0 + j * psy);
          fy *= fy;
          T fct = 0;
          auto tmp = 1 - fx - fy;
          if (tmp >= 0) {
            auto nm1 =
                (-fx - fy) / (sqrt(tmp) + 1);  // accurate form of sqrt(1-x-y)-1
            fct = T(kernel.corfac(nm1 * dw));
            if (divide_by_n) fct /= nm1 + 1;
          } else  // beyond the horizon, don't really know what to do here
          {
            if (divide_by_n)
              fct = 0;
            else {
              auto nm1 = sqrt(-tmp) - 1;
              fct = T(kernel.corfac(nm1 * dw));
            }
          }
          fct *= T(cfu[nx_dirty / 2 - i] * cfv[ny_dirty / 2 - j]);
          size_t i2 = nx_dirty - i, j2 = ny_dirty - j;
          dirty(i, j) *= fct;
          if ((i > 0) && (i < i2)) {
            dirty(i2, j) *= fct;
            if ((j > 0) && (j < j2)) dirty(i2, j2) *= fct;
          }
          if ((j > 0) && (j < j2)) dirty(i, j2) *= fct;
        }
      }
  });
}

template <typename Serv>
class WgridHelper {
 private:
  const Serv &srv;
  double wmin, dw;
  size_t nplanes, supp;
  vector<vector<idx_t>> minplane;
  size_t verbosity;

  int curplane;
  vector<idx_t> subidx;

  static void wminmax(const GridderConfig &, const Serv &srv, double &wmin,
                      double &wmax) {
    size_t nvis = srv.Nvis();
    if (nvis == 0) {
      wmin = 0.0;
      wmax = 0.0;
    } else {
      wmin = 1e38;
      wmax = -1e38;
      for (size_t ipart = 0; ipart < nvis; ++ipart) {
        auto wval = abs(srv.getCoord(ipart).w);
        wmin = min(wmin, wval);
        wmax = max(wmax, wval);
      }
    }
  }

  template <typename T>
  static void update_idx(vector<T> &v, const vector<T> &add,
                         const vector<T> &del) {
    myassert(v.size() >= del.size(), "must not happen");
    vector<T> res;
    res.reserve((v.size() + add.size()) - del.size());
    auto iin = v.begin(), ein = v.end();
    auto iadd = add.begin(), eadd = add.end();
    auto irem = del.begin(), erem = del.end();

    while (iin != ein) {
      if ((irem != erem) && (*iin == *irem)) {
        ++irem;
        ++iin;
      }  // skip removed entry
      else if ((iadd != eadd) && (*iadd < *iin))
        res.push_back(*(iadd++));  // add new entry
      else
        res.push_back(*(iin++));
    }
    myassert(irem == erem, "must not happen");
    while (iadd != eadd) res.push_back(*(iadd++));
    myassert(res.size() == (v.size() + add.size()) - del.size(),
             "must not happen");
    v.swap(res);
  }

 public:
  WgridHelper(const GridderConfig &gconf, const Serv &srv_, size_t verbosity_)
      : srv(srv_), verbosity(verbosity_), curplane(-1) {
    size_t nvis = srv.Nvis();
    size_t nthreads = gconf.Nthreads();
    double wmax;

    wminmax(gconf, srv, wmin, wmax);
    if (verbosity > 0) cout << "Using " << nthreads << " threads" << endl;
    if (verbosity > 0) cout << "W range: " << wmin << " to " << wmax << endl;

    double x0 = -0.5 * gconf.Nxdirty() * gconf.Pixsize_x(),
           y0 = -0.5 * gconf.Nydirty() * gconf.Pixsize_y();
    double nmin = sqrt(max(1. - x0 * x0 - y0 * y0, 0.)) - 1.;
    if (x0 * x0 + y0 * y0 > 1.) nmin = -sqrt(abs(1. - x0 * x0 - y0 * y0)) - 1.;
    dw = 0.25 / abs(nmin);
    nplanes = size_t((wmax - wmin) / dw + 2);
    if (nplanes <= 1)
      dw = 1;
    else
      dw = (1. + 1e-13) * (wmax - wmin) / (nplanes - 1);

    supp = gconf.Supp();
    wmin -= (0.5 * supp - 1) * dw;
    wmax += (0.5 * supp - 1) * dw;
    nplanes += supp - 2;
    if (verbosity > 0) cout << "Kernel support: " << supp << endl;
    if (verbosity > 0) cout << "nplanes: " << nplanes << endl;

    minplane.resize(nplanes);
#if 0
      // extra short, but potentially inefficient version:
      for (size_t ipart=0; ipart<nvis; ++ipart)
        {
        int plane0 = max(0,int(1+(abs(srv.getCoord(ipart).w)-(0.5*supp*dw)-wmin)/dw));
        minplane[plane0].push_back(idx_t(ipart));
        }
#else
    // more efficient: precalculate final vector sizes and avoid reallocations
    vector<size_t> cnt(nplanes, 0);
    for (size_t ipart = 0; ipart < nvis; ++ipart) {
      int plane0 = max(
          0, int(1 +
                 (abs(srv.getCoord(ipart).w) - (0.5 * supp * dw) - wmin) / dw));
      ++cnt[plane0];
    }

    // fill minplane
    for (size_t j = 0; j < nplanes; ++j) minplane[j].resize(cnt[j]);
    vector<size_t> ofs(nplanes, 0);
    for (size_t ipart = 0; ipart < nvis; ++ipart) {
      int plane0 = max(
          0, int(1 +
                 (abs(srv.getCoord(ipart).w) - (0.5 * supp * dw) - wmin) / dw));
      minplane[plane0][ofs[plane0]++] = idx_t(ipart);
    }
#endif
  }

  typename Serv::Tsub getSubserv() const {
    auto subidx2 = const_mav<idx_t, 1>(subidx.data(), {subidx.size()});
    return srv.getSubserv(subidx2);
  }
  double W() const { return wmin + curplane * dw; }
  size_t Nvis() const { return subidx.size(); }
  double DW() const { return dw; }
  bool advance() {
    if (++curplane >= int(nplanes)) return false;
    update_idx(
        subidx, minplane[curplane],
        curplane >= int(supp) ? minplane[curplane - supp] : vector<idx_t>());
    if (verbosity > 1)
      cout << "Working on plane " << curplane << " containing " << subidx.size()
           << " visibilities" << endl;
    return true;
  }
};

template <typename T, typename Serv>
void x2dirty(const GridderConfig &gconf, const Serv &srv,
             const mav<T, 2> &dirty, bool do_wstacking, size_t verbosity) {
  if (do_wstacking) {
    size_t nthreads = gconf.Nthreads();
    if (verbosity > 0) cout << "Gridding using improved w-stacking" << endl;
    WgridHelper<Serv> hlp(gconf, srv, verbosity);
    double dw = hlp.DW();
    dirty.fill(0);
    tmpStorage<complex<T>, 2> grid_({gconf.Nu(), gconf.Nv()});
    auto grid = grid_.getMav();
    while (hlp.advance())  // iterate over w planes
    {
      if (hlp.Nvis() == 0) continue;
      grid.fill(0);
      x2grid_c(gconf, hlp.getSubserv(), grid, hlp.W(), dw);
      gconf.grid2dirty_c_overwrite_wscreen_add(grid, dirty, T(hlp.W()));
    }
    // correct for w gridding etc.
    apply_global_corrections(gconf, dirty,
                             ES_Kernel(gconf.Supp(), gconf.Ofactor(), nthreads),
                             dw, false);
  } else
    myfail("not supported");
}

template <typename T, typename Serv>
void dirty2x(const GridderConfig &gconf, const const_mav<T, 2> &dirty,
             const Serv &srv, bool do_wstacking, size_t verbosity) {
  if (do_wstacking) {
    size_t nx_dirty = gconf.Nxdirty(), ny_dirty = gconf.Nydirty();
    size_t nthreads = gconf.Nthreads();
    if (verbosity > 0) cout << "Degridding using improved w-stacking" << endl;
    WgridHelper<Serv> hlp(gconf, srv, verbosity);
    double dw = hlp.DW();
    tmpStorage<T, 2> tdirty_({nx_dirty, ny_dirty});
    auto tdirty = tdirty_.getMav();
    for (size_t i = 0; i < nx_dirty; ++i)
      for (size_t j = 0; j < ny_dirty; ++j) tdirty(i, j) = dirty(i, j);
    // correct for w gridding etc.
    apply_global_corrections(gconf, tdirty,
                             ES_Kernel(gconf.Supp(), gconf.Ofactor(), nthreads),
                             dw, false);
    tmpStorage<complex<T>, 2> grid_({gconf.Nu(), gconf.Nv()});
    auto grid = grid_.getMav();
    while (hlp.advance())  // iterate over w planes
    {
      if (hlp.Nvis() == 0) continue;
      gconf.dirty2grid_c_wscreen(cmav(tdirty), grid, T(hlp.W()));
      grid2x_c(gconf, cmav(grid), hlp.getSubserv(), hlp.W(), dw);
    }
  } else
    myfail("not supported");
}

void calc_share(size_t nshares, size_t myshare, size_t nwork, size_t &lo,
                size_t &hi) {
  size_t nbase = nwork / nshares;
  size_t additional = nwork % nshares;
  lo = myshare * nbase + ((myshare < additional) ? myshare : additional);
  hi = lo + nbase + (myshare < additional);
}

template <typename T>
vector<idx_t> getWgtIndices(const Baselines &baselines,
                            const GridderConfig &gconf,
                            const const_mav<T, 2> &wgt,
                            const const_mav<complex<T>, 2> &ms) {
  size_t nrow = baselines.Nrows(), nchan = baselines.Nchannels(),
         nsafe = gconf.Nsafe();
  bool have_wgt = wgt.size() != 0;
  if (have_wgt) checkShape(wgt.shape(), {nrow, nchan});
  bool have_ms = ms.size() != 0;
  if (have_ms) checkShape(ms.shape(), {nrow, nchan});
  constexpr int side = 1 << logsquare;
  size_t nbu = (gconf.Nu() + 1 + side - 1) >> logsquare,
         nbv = (gconf.Nv() + 1 + side - 1) >> logsquare;
  vector<vector<idx_t>> acc;
  vector<idx_t> tmp(nrow * nchan);

  auto nthr = gconf.Nthreads();
  acc.resize(nthr);
  // aocommon::ParallelFor<size_t> loop(nthr);
  // loop.Run(0, nthr, [&](size_t, size_t thread) {
  for (size_t thread = 0; thread != nthr; ++thread) {
    size_t lo, hi;
    calc_share(nthr, thread, nrow, lo, hi);
    vector<idx_t> &lacc(acc[thread]);
    lacc.assign(nbu * nbv + 1, 0);
    for (idx_t irow = lo, idx = lo * nchan; irow < hi; ++irow)
      for (idx_t ichan = 0; ichan < nchan; ++ichan, ++idx)
        if (((!have_ms) || (norm(ms(irow, ichan)) != 0)) &&
            ((!have_wgt) || (wgt(irow, ichan) != 0))) {
          auto uvw = baselines.effectiveCoord(RowChan{irow, idx_t(ichan)});
          if (uvw.w < 0) uvw.Flip();
          double u, v;
          int iu0, iv0;
          gconf.getpix(uvw.u, uvw.v, u, v, iu0, iv0);
          iu0 = (iu0 + nsafe) >> logsquare;
          iv0 = (iv0 + nsafe) >> logsquare;
          ++lacc[nbv * iu0 + iv0 + 1];
          tmp[idx] = nbv * iu0 + iv0;  // not allowed in parallel
        } else
          tmp[idx] = ~idx_t(0);
  }
  for (size_t thread = 0; thread != nthr; ++thread) {
    size_t lo2, hi2;
    calc_share(nthr, thread, nbu * nbv, lo2, hi2);
    for (size_t i = lo2 + 1; i < hi2 + 1; ++i) {
      idx_t sum = 0;
      for (size_t j = 0; j != nthr; ++j) sum += acc[j][i];
      acc[0][i] = sum;  // Not allowed in parallel
    }
  }

  auto &acc0(acc[0]);
  for (size_t i = 1; i < acc0.size(); ++i) acc0[i] += acc0[i - 1];

  vector<idx_t> res(acc0.back());
  for (size_t irow = 0, idx = 0; irow < nrow; ++irow)
    for (size_t ichan = 0; ichan < nchan; ++ichan, ++idx)
      if (tmp[idx] != (~idx_t(0)))
        res[acc0[tmp[idx]]++] = baselines.getIdx(irow, ichan);
  return res;
}

template <typename T>
void ms2dirty_general(const const_mav<double, 2> &uvw,
                      const const_mav<double, 1> &freq,
                      const const_mav<complex<T>, 2> &ms,
                      const const_mav<T, 2> &wgt, double pixsize_x,
                      double pixsize_y, size_t nu, size_t nv, double epsilon,
                      bool do_wstacking, size_t nthreads,
                      const mav<T, 2> &dirty, size_t verbosity,
                      bool negate_v = false) {
  prep_fft<T>(nu, nv, nthreads, verbosity);
  Baselines baselines(uvw, freq, negate_v);
  GridderConfig gconf(dirty.shape(0), dirty.shape(1), nu, nv, epsilon,
                      pixsize_x, pixsize_y, nthreads);
  auto idx = getWgtIndices(baselines, gconf, wgt, ms);
  auto idx2 = const_mav<idx_t, 1>(idx.data(), {idx.size()});
  x2dirty(gconf, makeMsServ(baselines, idx2, ms, wgt), dirty, do_wstacking,
          verbosity);
}

template <typename T>
void dirty2ms_general(const const_mav<double, 2> &uvw,
                      const const_mav<double, 1> &freq,
                      const const_mav<T, 2> &dirty, const const_mav<T, 2> &wgt,
                      double pixsize_x, double pixsize_y, size_t nu, size_t nv,
                      double epsilon, bool do_wstacking, size_t nthreads,
                      const mav<complex<T>, 2> &ms, size_t verbosity,
                      bool negate_v = false) {
  prep_fft<T>(nu, nv, nthreads, verbosity);
  Baselines baselines(uvw, freq, negate_v);
  GridderConfig gconf(dirty.shape(0), dirty.shape(1), nu, nv, epsilon,
                      pixsize_x, pixsize_y, nthreads);
  const_mav<complex<T>, 2> null_ms(nullptr, {0, 0});
  auto idx = getWgtIndices(baselines, gconf, wgt, null_ms);
  auto idx2 = const_mav<idx_t, 1>(idx.data(), {idx.size()});
  ms.fill(0);
  dirty2x(gconf, dirty, makeMsServ(baselines, idx2, ms, wgt), do_wstacking,
          verbosity);
}

}  // namespace detail

// public names
using detail::CodeLocation;
using detail::const_mav;
using detail::dirty2ms_general;
using detail::mav;
using detail::ms2dirty_general;
using detail::streamDump__;
}  // namespace gridder

#endif
