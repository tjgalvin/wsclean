#include <boost/test/unit_test.hpp>

#include "../deconvolution/peakfinder.h"

#include <aocommon/uvector.h>

#include <random>

BOOST_AUTO_TEST_SUITE(peak_finder)

const size_t nRepeats =
    3; /* This should be set to 100 to assert the performance */

#if defined __AVX__ && !defined FORCE_NON_AVX
template <typename NumT>
struct CleanTestFixture {
  size_t x, y;
  aocommon::UVector<NumT> img;

  CleanTestFixture(size_t n = 16) : x(size_t(-1)), y(size_t(-1)), img(n, 0) {}
  void findPeak(size_t width = 4, size_t height = 2, size_t ystart = 0,
                size_t yend = 2) {
    PeakFinder::AVX(img.data(), width, height, x, y, true, ystart, yend, 0, 0);
  }
};
#endif

struct NoiseFixture {
  NoiseFixture() : n(2048), psf(n * n, 0.0), img(n * n), normal_dist(0.0, 1.0) {
    mt.seed(42);
    for (size_t i = 0; i != n * n; ++i) {
      img[i] = normal_dist(mt);
    }
  }

  size_t n;
  aocommon::UVector<double> psf, img;
  std::mt19937 mt;
  std::normal_distribution<double> normal_dist;
};

#if defined __AVX__ && !defined FORCE_NON_AVX
BOOST_AUTO_TEST_CASE(findPeakAVX1Double) {
  CleanTestFixture<double> f;
  f.img[0] = 1;
  f.findPeak();
  BOOST_CHECK_EQUAL(f.x, 0u);
  BOOST_CHECK_EQUAL(f.y, 0u);
}

BOOST_AUTO_TEST_CASE(findPeakAVX2Double) {
  CleanTestFixture<double> f;
  f.img[0] = 1;
  f.img[1] = 2;
  f.findPeak();
  BOOST_CHECK_EQUAL(f.x, 1u);
  BOOST_CHECK_EQUAL(f.y, 0u);
}

BOOST_AUTO_TEST_CASE(findPeakAVX3Double) {
  CleanTestFixture<double> f;
  f.img[0] = 1;
  f.img[1] = 2;
  f.img[4] = 3;
  f.findPeak();
  BOOST_CHECK_EQUAL(f.x, 0u);
  BOOST_CHECK_EQUAL(f.y, 1u);
}

BOOST_AUTO_TEST_CASE(findPeakAVX4Double) {
  CleanTestFixture<double> f;
  f.img[0] = 1;
  f.img[1] = 2;
  f.img[4] = 3;
  f.img[7] = 4;
  f.findPeak();
  BOOST_CHECK_EQUAL(f.x, 3u);
  BOOST_CHECK_EQUAL(f.y, 1u);
}

BOOST_AUTO_TEST_CASE(findPeakAVX5Double) {
  CleanTestFixture<double> f;
  f.img[0] = 1;
  f.img[1] = 2;
  f.img[4] = 3;
  f.img[7] = 4;
  f.img[15] = 6;
  f.findPeak(4, 4, 0, 4);
  BOOST_CHECK_EQUAL(f.x, 3u);
  BOOST_CHECK_EQUAL(f.y, 3u);
}

BOOST_AUTO_TEST_CASE(findPeakAVX6Double) {
  CleanTestFixture<double> f;
  f.img[0] = 1;
  f.img[1] = 2;
  f.img[4] = 3;
  f.img[7] = 4;
  f.img[15] = 6;
  f.img[14] = 5;
  f.findPeak(3, 5, 0, 5);
  BOOST_CHECK_EQUAL(f.x, 2u);
  BOOST_CHECK_EQUAL(f.y, 4u);
}

BOOST_AUTO_TEST_CASE(findPeakAVX1Single) {
  CleanTestFixture<float> f(36);
  f.img[0] = 1;
  f.findPeak(6, 6, 0, 6);
  BOOST_CHECK_EQUAL(f.x, 0u);
  BOOST_CHECK_EQUAL(f.y, 0u);
}

BOOST_AUTO_TEST_CASE(findPeakAVX2Single) {
  CleanTestFixture<float> f(36);
  f.img[0] = 1;
  f.img[1] = 2;
  f.findPeak(6, 6, 0, 6);
  BOOST_CHECK_EQUAL(f.x, 1u);
  BOOST_CHECK_EQUAL(f.y, 0u);
}

BOOST_AUTO_TEST_CASE(findPeakAVX3Single) {
  CleanTestFixture<float> f(36);
  f.img[0] = 1;
  f.img[1] = 2;
  f.img[6] = 3;
  f.findPeak(6, 6, 0, 6);
  BOOST_CHECK_EQUAL(f.x, 0u);
  BOOST_CHECK_EQUAL(f.y, 1u);
}

BOOST_AUTO_TEST_CASE(findPeakAVX4Single) {
  CleanTestFixture<float> f(36);
  f.img[0] = 1;
  f.img[1] = 2;
  f.img[6] = 3;
  f.img[9] = 4;
  f.findPeak(6, 6, 0, 6);
  BOOST_CHECK_EQUAL(f.x, 3u);
  BOOST_CHECK_EQUAL(f.y, 1u);
}

BOOST_AUTO_TEST_CASE(findPeakAVX5Single) {
  CleanTestFixture<float> f(36);
  f.img[0] = 1;
  f.img[1] = 2;
  f.img[6] = 3;
  f.img[9] = 4;
  f.img[35] = 6;
  f.findPeak(6, 6, 0, 6);
  BOOST_CHECK_EQUAL(f.x, 5u);
  BOOST_CHECK_EQUAL(f.y, 5u);
}

BOOST_AUTO_TEST_CASE(findPeakAVX6Single) {
  CleanTestFixture<float> f(36);
  f.img[0] = 1;
  f.img[1] = 2;
  f.img[6] = 3;
  f.img[9] = 4;
  f.img[35] = 6;
  f.findPeak(2, 18, 0, 18);
  BOOST_CHECK_EQUAL(f.x, 1u);
  BOOST_CHECK_EQUAL(f.y, 17u);
}

BOOST_AUTO_TEST_CASE(findPeakAVX7Single) {
  CleanTestFixture<float> f(38);
  f.img[0] = 1;
  f.img[1] = 2;
  f.img[6] = 3;
  f.img[9] = 4;
  f.img[35] = 6;
  f.img[37] = 7;
  f.findPeak(6, 6, 0, 6);
  BOOST_CHECK_EQUAL(f.x, 5u);
  BOOST_CHECK_EQUAL(f.y, 5u);
}

#endif

BOOST_AUTO_TEST_CASE(findPeakPerformanceDouble) {
  NoiseFixture f;
  for (size_t repeat = 0; repeat != nRepeats; ++repeat) {
    size_t x, y;
    PeakFinder::Find(f.img.data(), f.n, f.n, x, y, true, 0, f.n / 2, 0.0);
  }
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_CASE(findPeakSimplePerformanceDouble) {
  NoiseFixture f;
  for (size_t repeat = 0; repeat != nRepeats; ++repeat) {
    size_t x, y;
    PeakFinder::Simple(f.img.data(), f.n, f.n, x, y, true, 0, f.n / 2, 0, 0);
  }
  BOOST_CHECK(true);
}

#if defined __AVX__ && !defined FORCE_NON_AVX
BOOST_AUTO_TEST_CASE(findPeakAVXPerformanceDouble) {
  NoiseFixture f;
  for (size_t repeat = 0; repeat != nRepeats; ++repeat) {
    size_t x, y;
    PeakFinder::AVX(f.img.data(), f.n, f.n, x, y, true, 0, f.n / 2, 0, 0);
  }
  BOOST_CHECK(true);
}
#endif

BOOST_AUTO_TEST_SUITE_END()
