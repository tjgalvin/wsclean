#include <boost/test/unit_test.hpp>

#include <schaapcommon/threading/recursivefor.h>
#include <schaapcommon/threading/staticfor.h>
#include <schaapcommon/threading/threadpool.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>

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
  schaapcommon::ThreadPool::GetInstance().SetNThreads(2);
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

  schaapcommon::ThreadPool::GetInstance().SetNThreads(2);
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

BOOST_AUTO_TEST_CASE(mask_noise_threshold_is_below_the_smallest_real_value) {
  // The smallest value that a mask convolved with a tophat kernel can have is
  // one divided by the number of pixels of the kernel, because a single
  // kernel pixel overlapping the edge of the mask is the smallest possible
  // contribution. The noise threshold must be below that value, otherwise
  // real mask data would be removed.
  for (size_t feather_size : {1, 5, 10, 50, 100}) {
    const size_t kernel_size = 2 * feather_size + 1;
    const Image kernel =
        tophat_convolution::MakeTopHatImage(kernel_size, kernel_size, feather_size);
    // The kernel is normalized such that it sums to one, so the smallest
    // non-zero value is one divided by the number of pixels in the kernel.
    float smallest_real_value = std::numeric_limits<float>::max();
    for (float value : kernel) {
      if (value != 0.0f) smallest_real_value = std::min(smallest_real_value, value);
    }
    BOOST_CHECK_LT(tophat_convolution::MaskNoiseThreshold(feather_size), smallest_real_value);
  }
}

BOOST_AUTO_TEST_CASE(remove_mask_noise_removes_the_fft_rounding_noise) {
  constexpr size_t kWidth = 128;
  constexpr size_t kHeight = 128;
  constexpr size_t kRadius = 10;
  Image image(kWidth, kHeight, 0.0f);
  image[kWidth * 32 + 32] = 1.0f;
  tophat_convolution::Convolve(image, kRadius);

  // The convolution is performed with FFTs, which leave rounding noise of
  // either sign in the whole image.
  size_t negative_value_count = 0;
  for (float value : image) {
    if (value < 0.0f) ++negative_value_count;
  }
  BOOST_CHECK_GT(negative_value_count, 0u);

  const float threshold = tophat_convolution::MaskNoiseThreshold(kRadius);
  tophat_convolution::RemoveMaskNoise(image, kRadius);
  for (float value : image) {
    BOOST_CHECK_GE(value, 0.0f);
    if (value != 0.0f) {
      BOOST_CHECK_GE(value, threshold);
    }
  }
}

BOOST_AUTO_TEST_CASE(convolve_from_nested_parallel_loop) {
  constexpr size_t kWidth = 35;
  constexpr size_t kHeight = 35;
  constexpr double kRadius = 3.01;
  constexpr size_t kOuterTasks = 4;

  schaapcommon::ThreadPool::GetInstance().SetNThreads(4);

  std::atomic<bool> all_finite{true};
  BOOST_CHECK_NO_THROW({
    schaapcommon::RecursiveFor recursive_for;
    schaapcommon::StaticFor<size_t> loop;
    loop.ConstrainedRun(0, kOuterTasks, 2, [&](size_t start, size_t end) {
      for (size_t iter = start; iter != end; ++iter) {
        Image image(kWidth, kHeight, 0.0);
        image[(kHeight / 2) * kWidth + (kWidth / 2)] = static_cast<float>(iter + 1);
        tophat_convolution::Convolve(image, kRadius);
        if (!std::all_of(image.begin(), image.end(),
                         [](float value) { return std::isfinite(value); })) {
          all_finite = false;
        }
      }
    });
  });
  BOOST_CHECK(all_finite.load());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace wsclean
