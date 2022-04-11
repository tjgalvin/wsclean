#include <boost/test/unit_test.hpp>

#include "../structures/primarybeamimageset.h"

#include <aocommon/image.h>

#include <radler/component_list.h>

namespace {
const size_t kWidth = 4;
const size_t kHeight = 4;
const size_t kNumScales = 2;
const size_t kNumFreqs = 3;
}  // namespace

BOOST_AUTO_TEST_SUITE(primarybeam_imageset)

BOOST_AUTO_TEST_CASE(correct_component_list) {
  PrimaryBeamImageSet imageSet(kWidth, kHeight);
  imageSet.SetToZero();

  BOOST_CHECK_EQUAL(imageSet.NImages(), 16u);
  BOOST_CHECK_EQUAL(imageSet.Width(), kWidth);
  BOOST_CHECK_EQUAL(imageSet.Height(), kHeight);

  aocommon::Image refImage(kWidth, kHeight, 1.0);
  for (size_t pixel = 0; pixel != kWidth * kHeight; ++pixel) {
    refImage[pixel] *= static_cast<float>(pixel);
  }

  for (size_t i = 0; i != imageSet.NImages(); ++i) {
    // Only fill diagonal, to make sure underlying HMC4x4 matrix is invertible
    if (i == 0 || i == 3 || i == 8 || i == 15) {
      imageSet[i] = refImage;
    }
  }

  for (size_t pixel = 0; pixel != kWidth * kHeight; ++pixel) {
    const size_t x = pixel % kWidth;
    const size_t y = pixel / kWidth;
    const double correctionFactor =
        imageSet.GetUnpolarizedCorrectionFactor(x, y);
    if (pixel != 0) {
      BOOST_CHECK_CLOSE(correctionFactor, 1.0f / static_cast<float>(pixel),
                        1e-5);
    } else {
      BOOST_CHECK_EQUAL(correctionFactor, 0.0f);
    }
  }

  radler::ComponentList list(kWidth, kHeight, kNumScales, kNumFreqs);
  const std::vector<radler::ComponentList::Position> positions = {
      radler::ComponentList::Position(1, 2),
      radler::ComponentList::Position(3, 3)};

  const std::vector<std::vector<float>> values{{1.0, 2.0, 3.0},
                                               {5.0, 6.0, 7.0}};

  // Add components
  for (size_t i = 0; i != kNumScales; ++i) {
    list.Add(positions[i].x, positions[i].y, i, values[i].data());
  }
  list.MergeDuplicates();

  // Correct components for primary beam
  for (size_t channel = 0; channel != list.NFrequencies(); ++channel) {
    imageSet.CorrectComponentList(list, channel);
  }

  for (size_t i = 0; i != list.NScales(); ++i) {
    BOOST_CHECK_EQUAL(list.ComponentCount(i), 1);

    size_t x;
    size_t y;
    aocommon::UVector<float> corrected_values(3);
    list.GetComponent(i, 0, x, y, corrected_values.data());
    const float pixel = static_cast<float>(y * kWidth + x);
    for (size_t channel = 0; channel != list.NFrequencies(); ++channel) {
      BOOST_CHECK_CLOSE(corrected_values[channel], values[i][channel] / pixel,
                        1e-5);
    }
  }
}
BOOST_AUTO_TEST_SUITE_END()