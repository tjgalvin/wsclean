#include <boost/test/unit_test.hpp>

#include <aocommon/threadpool.h>

#include "../../math/tophatconvolution.h"

using aocommon::Image;

namespace wsclean {

namespace {
constexpr float v = 1.0f / 29.0f;
const Image reference_7x7_radius_3(7, 7,
                                   {
                                       0.0f, 0.0f, 0.0f, v, 0.0f, 0.0f, 0.0f,  //
                                       0.0f, v,    v,    v, v,    v,    0.0f,  //
                                       0.0f, v,    v,    v, v,    v,    0.0f,  //
                                       v,    v,    v,    v, v,    v,    v,     //
                                       0.0f, v,    v,    v, v,    v,    0.0f,  //
                                       0.0f, v,    v,    v, v,    v,    0.0f,  //
                                       0.0f, 0.0f, 0.0f, v, 0.0f, 0.0f, 0.0f   //
                                   });

float ref_value(ssize_t x, ssize_t y) {
  const ssize_t width = reference_7x7_radius_3.Width();
  const ssize_t height = reference_7x7_radius_3.Height();
  if (x < 0 || y < 0 || x >= width || y >= height)
    return 0.0;
  else
    return reference_7x7_radius_3[x + y * width];
}

}  // namespace

BOOST_AUTO_TEST_SUITE(tophat_convolution)

BOOST_AUTO_TEST_CASE(make_tophat_with_large_radius) {
  constexpr size_t kWidth = 10;
  constexpr size_t kHeight = 12;
  constexpr double kRadius = 15.0;
  const Image image = tophat_convolution::MakeTopHatImage(kWidth, kHeight, kRadius);
  for (float value : image) {
    BOOST_CHECK_CLOSE_FRACTION(value, 1.0 / (kWidth * kHeight), 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(make_single_pixel_tophat) {
  constexpr size_t kWidth = 3;
  constexpr size_t kHeight = 3;
  constexpr double kRadius = 0.0;
  const Image image = tophat_convolution::MakeTopHatImage(kWidth, kHeight, kRadius);
  for (size_t i = 0; i != kWidth * kHeight; ++i) {
    const float expected = (i == kWidth * ((kHeight - 1) / 2) + ((kWidth - 1) / 2)) ? 1.0f : 0.0f;
    BOOST_CHECK_CLOSE_FRACTION(image[i], expected, 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(make_medium_tophat) {
  constexpr size_t kWidth = 7;
  constexpr size_t kHeight = 7;
  // Use value slightly larger than 3 to avoid rounding errors for pixels that
  // are exactly a distance of 3 pixels away.
  constexpr double kRadius = 3.01;
  const Image image = tophat_convolution::MakeTopHatImage(kWidth, kHeight, kRadius);
  for (size_t i = 0; i != kWidth * kHeight; ++i) {
    BOOST_CHECK_CLOSE_FRACTION(image[i], reference_7x7_radius_3[i], 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(convolve_single_pixel) {
  constexpr size_t kWidth = 7;
  constexpr size_t kHeight = 7;
  constexpr double kRadius = 3.01;
  Image image(kWidth, kHeight, 0.0);
  image[kWidth * ((kHeight - 1) / 2) + (kWidth - 1) / 2] = 3.0;
  aocommon::ThreadPool::GetInstance().SetNThreads(2);
  tophat_convolution::Convolve(image, kRadius);
  for (size_t i = 0; i != kWidth * kHeight; ++i) {
    // BOOST_CHECK_CLOSE_FRACTION doesn't work when comparing zeros to tiny
    // non-zero values
    BOOST_CHECK_LT(std::fabs(image[i] - 3.0 * reference_7x7_radius_3[i]), 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(convolve_multi_pixel) {
  constexpr size_t kWidth = 35;
  constexpr size_t kHeight = 35;
  constexpr double kRadius = 3.01;

  const ssize_t mid_x = (kWidth - 1) / 2;
  const ssize_t mid_y = (kHeight - 1) / 2;
  const ssize_t ref_mid_x = (reference_7x7_radius_3.Width() - 1) / 2;
  const ssize_t ref_mid_y = (reference_7x7_radius_3.Height() - 1) / 2;

  Image image(kWidth, kHeight, 0.0);
  image[mid_y * kWidth + mid_x] = 3.0;
  image[(mid_y + 1) * kWidth + mid_x + 1] = 1.0;
  image[(mid_y + 2) * kWidth + mid_x + 1] = -2.0;

  aocommon::ThreadPool::GetInstance().SetNThreads(2);
  tophat_convolution::Convolve(image, kRadius);
  for (ssize_t y = 0; y != kHeight; ++y) {
    for (ssize_t x = 0; x != kWidth; ++x) {
      const float value_a = ref_value(x - mid_x + ref_mid_x, y - mid_y + ref_mid_y) * 3.0;
      const float value_b = ref_value(x - 1 - mid_x + ref_mid_x, y - 1 - mid_y + ref_mid_y) * 1.0;
      const float value_c = ref_value(x - 1 - mid_x + ref_mid_x, y - 2 - mid_y + ref_mid_y) * -2.0;
      const float reference = value_a + value_b + value_c;
      BOOST_CHECK_LT(std::fabs(image[y * kWidth + x] - reference), 1e-6);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace wsclean
