#include <boost/test/unit_test.hpp>

#include <numeric>

#include "../../gridding/griddingutils.h"

namespace wsclean {

namespace {
void VerifyTotal(const std::vector<size_t>& sizes, size_t expected) {
  const size_t total = std::accumulate(sizes.begin(), sizes.end(), 0);
  BOOST_CHECK_EQUAL(total, expected);
}
}  // namespace

BOOST_AUTO_TEST_SUITE(CalculateChunkSizesTests)

BOOST_AUTO_TEST_CASE(TwoChunksWithRemainder) {
  const std::vector<size_t> sizes = CalculateChunkSizes(1000, 1000, 1.5f);
  BOOST_CHECK_EQUAL(sizes.size(), 2);
  BOOST_CHECK_EQUAL(sizes[0], 334);
  BOOST_CHECK_EQUAL(sizes[1], 666);
  VerifyTotal(sizes, 1000);
}

BOOST_AUTO_TEST_CASE(MultipleChunks) {
  const std::vector<size_t> sizes = CalculateChunkSizes(10000, 1000, 1.5f);
  BOOST_CHECK_EQUAL(sizes.size(), 16);
  BOOST_CHECK_EQUAL(sizes[0], 10);
  VerifyTotal(sizes, 10000);
}

BOOST_AUTO_TEST_CASE(VariousScenarios) {
  const std::vector<std::tuple<size_t, size_t, float>> test_cases = {
      {1000, 1000, 1.5f}, {5000, 2000, 1.5f}, {1000, 1000, 2.0f}, {500, 10000, 1.5f}};

  for (const auto& [total_rows, max_memory, factor] : test_cases) {
    const std::vector<size_t> sizes = CalculateChunkSizes(total_rows, max_memory, factor);
    BOOST_CHECK_GT(sizes.size(), 0);
    VerifyTotal(sizes, total_rows);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace wsclean
