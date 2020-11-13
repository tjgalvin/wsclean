#include <boost/test/unit_test.hpp>

#include "../deconvolution/simpleclean.h"

#include <aocommon/uvector.h>

#include <random>

BOOST_AUTO_TEST_SUITE(clean_algorithms)

template <typename NumT>
struct NoiseFixture {
  NoiseFixture() : n(2048), psf(n * n, 0.0), img(n * n), normal_dist(0.0, 1.0) {
    mt.seed(42);
    for (size_t i = 0; i != n * n; ++i) {
      img[i] = normal_dist(mt);
    }
  }

  size_t n;
  aocommon::UVector<NumT> psf, img;
  std::mt19937 mt;
  std::normal_distribution<NumT> normal_dist;
};

const size_t nRepeats =
    3; /* This should be set to 100 to assert the performance */

BOOST_AUTO_TEST_CASE(partialSubtractImagePerformance) {
  NoiseFixture<float> f;
  size_t x = f.n / 2, y = f.n / 2;
  for (size_t repeat = 0; repeat != nRepeats; ++repeat)
    SimpleClean::PartialSubtractImage(f.img.data(), f.n, f.n, f.psf.data(), f.n,
                                      f.n, x, y, 0.5, 0, f.n / 2);
  BOOST_CHECK(true);
}

#if defined __AVX__ && !defined FORCE_NON_AVX
BOOST_AUTO_TEST_CASE(partialSubtractImageAVXPerformance) {
  NoiseFixture<double> f;
  size_t x = f.n / 2, y = f.n / 2;
  for (size_t repeat = 0; repeat != nRepeats; ++repeat)
    SimpleClean::PartialSubtractImageAVX(f.img.data(), f.n, f.n, f.psf.data(),
                                         f.n, f.n, x, y, 0.5, 0, f.n / 2);
  BOOST_CHECK(true);
}
#endif

BOOST_AUTO_TEST_SUITE_END()
