#include "../../structures/facetutil.h"

#include <boost/test/unit_test.hpp>

using schaapcommon::facets::BoundingBox;
using schaapcommon::facets::Facet;

BOOST_AUTO_TEST_SUITE(facet_util)

BOOST_AUTO_TEST_CASE(create_grid_single_facet) {
  const size_t kImageWidth = 100;
  const size_t kImageHeight = 142;
  Facet::InitializationData facet_data(0.01, 0.01, kImageWidth, kImageHeight);

  const std::vector<schaapcommon::facets::Facet> facets =
      CreateFacetGrid(facet_data, 1, 1);

  BOOST_TEST_REQUIRE(facets.size() == 1);
  const BoundingBox& box = facets.front().GetTrimmedBoundingBox();
  BOOST_TEST(box.Min().x == 0);
  BOOST_TEST(box.Min().y == 0);
  BOOST_TEST(box.Max().x == kImageWidth);
  BOOST_TEST(box.Max().y == kImageHeight);

  BOOST_TEST(facets.front().DirectionLabel() == "0, 0");
}

BOOST_AUTO_TEST_CASE(create_grid_multiple_facets) {
  const size_t kImageSize = 100;
  const size_t kGridWidth = 4;
  const size_t kGridHeight = 5;
  Facet::InitializationData facet_data(0.01, kImageSize);
  const std::vector<schaapcommon::facets::Facet> facets =
      CreateFacetGrid(facet_data, kGridWidth, kGridHeight);

  BOOST_TEST(facets.size() == kGridWidth * kGridHeight);

  std::set<std::pair<int, int>> grid_cells;
  for (const Facet& facet : facets) {
    const BoundingBox& box = facet.GetTrimmedBoundingBox();
    const int grid_x = box.Centre().x / box.Width();
    const int grid_y = box.Centre().y / box.Height();
    BOOST_TEST(grid_x >= 0);
    BOOST_TEST(grid_y >= 0);
    BOOST_TEST(grid_x < static_cast<int>(kGridWidth));
    BOOST_TEST(grid_y < static_cast<int>(kGridHeight));
    grid_cells.emplace(grid_x, grid_y);

    BOOST_TEST(box.Min().x == grid_x * kImageSize / kGridWidth);
    BOOST_TEST(box.Min().y == grid_y * kImageSize / kGridHeight);
    BOOST_TEST(box.Max().x == (grid_x + 1) * kImageSize / kGridWidth);
    BOOST_TEST(box.Max().y == (grid_y + 1) * kImageSize / kGridHeight);

    BOOST_TEST(facet.DirectionLabel() ==
               std::to_string(grid_x) + ", " + std::to_string(grid_y));
  }

  // Check that the generated grid contains all grid cells.
  BOOST_TEST(grid_cells.size() == kGridWidth * kGridHeight);
}

BOOST_AUTO_TEST_SUITE_END()