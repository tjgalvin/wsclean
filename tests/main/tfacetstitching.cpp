#include <boost/test/unit_test.hpp>

#include <aocommon/image.h>
#include <aocommon/imagecoordinates.h>

#include <schaapcommon/facets/facet.h>
#include <schaapcommon/facets/facetimage.h>

#include "../../main/facetstitching.h"
#include "../../math/tophatconvolution.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>
#include <vector>

using aocommon::Image;
using schaapcommon::facets::Coord;
using schaapcommon::facets::Facet;
using schaapcommon::facets::FacetImage;

namespace wsclean {

namespace {

// The facets are placed well inside the image, so that the mask of a facet is
// padded before it is convolved. The gap between the facets is smaller than
// twice the feather size, so that the padded boxes of two diagonally
// neighbouring facets overlap. In the overlap only the corners of the boxes
// are covered, which is where the FFT rounding noise of the two facets can
// cancel out.
constexpr int kImageSize = 1024;
constexpr size_t kFeatherSize = 100;
constexpr int kFacetMargin = 200;
constexpr int kFacetGap = 190;
constexpr double kPixelScale = 1e-4;
constexpr double kPhaseCentreRA = 0.5;
constexpr double kPhaseCentreDec = 0.3;

Coord PixelToCoord(const Facet::InitializationData& data, int x, int y) {
  double l, m;
  aocommon::ImageCoordinates::XYToLM(x, y, data.pixel_scale_x, data.pixel_scale_y, data.image_width,
                                     data.image_height, l, m);
  double ra, dec;
  aocommon::ImageCoordinates::LMToRaDec(l, m, data.phase_centre.ra, data.phase_centre.dec, ra, dec);
  return Coord(ra, dec);
}

Facet::InitializationData MakeInitializationData() {
  Facet::InitializationData data(kPixelScale, kPixelScale, kImageSize, kImageSize);
  data.phase_centre.ra = kPhaseCentreRA;
  data.phase_centre.dec = kPhaseCentreDec;
  data.feather_size = kFeatherSize;
  data.align = 1;
  data.make_square = false;
  data.padding = 1.0;
  return data;
}

std::shared_ptr<Facet> MakeRectangularFacet(const Facet::InitializationData& data, int x0, int y0,
                                            int x1, int y1) {
  const std::vector<Coord> corners{PixelToCoord(data, x0, y0), PixelToCoord(data, x1, y0),
                                   PixelToCoord(data, x1, y1), PixelToCoord(data, x0, y1)};
  return std::make_shared<Facet>(data, corners);
}

/**
 * The smallest value that the convolved mask can legitimately have: a single
 * pixel of the tophat kernel that overlaps the mask. The kernel is normalised
 * such that its values sum to one, so this is the value of the smallest pixel
 * of the kernel.
 */
double SmallestMaskValue(size_t feather_size) {
  const Image kernel =
      tophat_convolution::MakeTopHatImage(2 * feather_size + 1, 2 * feather_size + 1, feather_size);
  float smallest = 1.0f;
  for (float value : kernel) {
    if (value != 0.0f) {
      smallest = std::min(smallest, value);
    }
  }
  return smallest;
}

struct StitchResult {
  /// The sum of the masks of the facets; the stitched image is divided by this.
  Image weight;
  /// The stitched image.
  Image image;
  /// The smallest non-zero value of the weight image.
  double min_nonzero_weight = 0.0;
  /// The largest absolute value of the stitched image.
  double max_abs_image = 0.0;
};

/**
 * Stitch the given facets, in the same way as StitchFacets() does: the facet
 * images are added to the full image with their mask as weight, and the full
 * image is divided by the sum of the masks.
 */
StitchResult Stitch(const std::vector<std::shared_ptr<Facet>>& facets, size_t feather_size) {
  FacetImage facet_image(kImageSize, kImageSize, 1);
  StitchResult result;
  result.weight = Image(kImageSize, kImageSize, 0.0f);
  Image numerator(kImageSize, kImageSize, 0.0f);

  // Use deterministic pseudo-random data with a unit standard deviation, so
  // that the noise of the image is of the same order as the noise of a real
  // observation.
  std::mt19937 generator(42);
  std::normal_distribution<float> normal(0.0f, 1.0f);

  for (const std::shared_ptr<Facet>& facet : facets) {
    facet_image.SetFacet(*facet, true);
    const Image mask = CreateFeatheredFacetMask(facet_image, *facet, feather_size);
    float* facet_data = facet_image.Data(0);
    const size_t facet_size = facet_image.Width() * facet_image.Height();
    for (size_t i = 0; i != facet_size; ++i) {
      facet_data[i] = normal(generator);
    }
    facet_image.AddWithMask(numerator.Data(), result.weight.Data(), mask, 0);
  }

  result.image = Image(kImageSize, kImageSize, 0.0f);
  for (size_t i = 0; i != result.weight.Size(); ++i) {
    if (result.weight[i] != 0.0f) {
      result.image[i] = numerator[i] / result.weight[i];
      const double abs_image = std::fabs(result.image[i]);
      result.max_abs_image = std::max(result.max_abs_image, abs_image);
    }
    if (result.weight[i] != 0.0f) {
      const double weight = std::fabs(result.weight[i]);
      if (result.min_nonzero_weight == 0.0 || weight < result.min_nonzero_weight) {
        result.min_nonzero_weight = weight;
      }
    }
  }
  return result;
}

/**
 * Four rectangular facets in a two-by-two grid, with a gap between them.
 */
std::vector<std::shared_ptr<Facet>> MakeFacetGrid() {
  const Facet::InitializationData data = MakeInitializationData();
  const int facet_width = (kImageSize - 2 * kFacetMargin - kFacetGap) / 2;
  std::vector<std::shared_ptr<Facet>> facets;
  for (int grid_y = 0; grid_y != 2; ++grid_y) {
    for (int grid_x = 0; grid_x != 2; ++grid_x) {
      const int x0 = kFacetMargin + grid_x * (facet_width + kFacetGap);
      const int y0 = kFacetMargin + grid_y * (facet_width + kFacetGap);
      const int x1 = (grid_x == 1) ? kImageSize : x0 + facet_width;
      const int y1 = (grid_y == 1) ? kImageSize : y0 + facet_width;
      facets.push_back(MakeRectangularFacet(data, x0, y0, x1, y1));
    }
  }
  return facets;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(facet_stitching)

BOOST_AUTO_TEST_CASE(feathered_mask_contains_no_fft_noise) {
  const Facet::InitializationData data = MakeInitializationData();
  const std::shared_ptr<Facet> facet = MakeRectangularFacet(data, kFacetMargin, kFacetMargin,
                                                            kFacetMargin + 200, kFacetMargin + 200);
  FacetImage facet_image(kImageSize, kImageSize, 1);
  facet_image.SetFacet(*facet, true);
  const Image mask = CreateFeatheredFacetMask(facet_image, *facet, kFeatherSize);

  // The convolution is performed with FFTs, which leave rounding noise of
  // either sign in the whole mask, including in the pixels where the mask
  // should be exactly zero. Every value of the mask should therefore either
  // be exactly zero, or be a real mask value.
  // The convolution has a rounding error of about 1e-7, which is 0.3% of
  // the smallest mask value, so the values are compared with a margin.
  const double smallest_value = SmallestMaskValue(kFeatherSize);
  for (float value : mask) {
    BOOST_CHECK_GE(value, 0.0f);
    if (value != 0.0f) {
      BOOST_CHECK_GE(value, 0.9 * smallest_value);
    }
  }
}

BOOST_AUTO_TEST_CASE(stitching_does_not_divide_by_a_near_zero_weight) {
  const std::vector<std::shared_ptr<Facet>> facets = MakeFacetGrid();

  // The rounding noise of two facets cancels out almost completely at some
  // pixel, which makes the division of the stitched image by the weight blow
  // up.
  const StitchResult result = Stitch(facets, kFeatherSize);
  BOOST_CHECK_GE(result.min_nonzero_weight, 0.9 * SmallestMaskValue(kFeatherSize));
  // The weights are non-negative, so the stitched image is a weighted average
  // of the facet images and cannot be much larger than they are.
  BOOST_CHECK_LT(result.max_abs_image, 10.0);
}

BOOST_AUTO_TEST_CASE(feathered_mask_of_an_edge_facet_is_not_contaminated) {
  // The trimmed bounding box of a facet that touches the edge of the image is
  // clipped by that edge, while its convolution box is not. The tophat kernel
  // of the convolution therefore extends beyond the mask, and the mask has to
  // be padded to keep the FFT convolution from wrapping around. If the mask is
  // not padded, the convolution wraps around and the mask is contaminated with
  // values from the other side of the facet.
  //
  // The mask of a facet does not depend on the position of the facet in the
  // image, so the mask of a facet at the edge of the image has to be equal to
  // the mask of the same facet in the interior of the image.
  const Facet::InitializationData data = MakeInitializationData();
  const std::shared_ptr<Facet> edge_facet = MakeRectangularFacet(data, 0, 0, 300, 300);
  const std::shared_ptr<Facet> inner_facet = MakeRectangularFacet(data, 500, 500, 800, 800);

  FacetImage facet_image(kImageSize, kImageSize, 1);
  facet_image.SetFacet(*edge_facet, true);
  const Image edge_mask = CreateFeatheredFacetMask(facet_image, *edge_facet, kFeatherSize);
  facet_image.SetFacet(*inner_facet, true);
  const Image inner_mask = CreateFeatheredFacetMask(facet_image, *inner_facet, kFeatherSize);

  // The trimmed bounding box of the edge facet is clipped by the edge of the
  // image, so it is smaller than the bounding box of the inner facet. The
  // inner facet is offset by the feather size with respect to the edge facet.
  BOOST_REQUIRE_EQUAL(edge_mask.Width(), 400u);
  BOOST_REQUIRE_EQUAL(edge_mask.Height(), 400u);
  BOOST_REQUIRE_EQUAL(inner_mask.Width(), 500u);
  BOOST_REQUIRE_EQUAL(inner_mask.Height(), 500u);

  for (size_t y = 0; y != edge_mask.Height(); ++y) {
    for (size_t x = 0; x != edge_mask.Width(); ++x) {
      const float edge_value = edge_mask[x + y * edge_mask.Width()];
      const float inner_value =
          inner_mask[(x + kFeatherSize) + (y + kFeatherSize) * inner_mask.Width()];
      // Both masks are zero outside of the feathered facet, and the values
      // inside the feather are only accurate to within the rounding error of
      // the FFTs.
      BOOST_CHECK_SMALL(edge_value - inner_value, 1e-5f);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace wsclean
