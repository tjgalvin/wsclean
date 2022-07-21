
#include "../../structures/resources.h"

#include <aocommon/logger.h>

#include <boost/test/unit_test.hpp>

namespace {
const int64_t gb = 1024u * 1024u * 1024u;
}

BOOST_AUTO_TEST_SUITE(resources)

BOOST_AUTO_TEST_CASE(get_part) {
  const Resources r(8, 8 * gb);

  const Resources quarter = r.GetPart(4);
  BOOST_TEST(quarter.NCpus() == 2);
  BOOST_TEST(quarter.Memory() == 2 * gb);

  // Always at least one cpu?
  BOOST_TEST(r.GetPart(16).NCpus() == 1);

  // Round up?
  BOOST_TEST(r.GetPart(3).NCpus() == 3);
}

BOOST_AUTO_TEST_CASE(get_available_memory) {
  aocommon::Logger::SetVerbosity(aocommon::Logger::kQuietVerbosity);
  const int64_t all_memory = GetAvailableMemory(1.0, 0.0);
  BOOST_TEST(all_memory > 0);

  const int64_t half = GetAvailableMemory(0.5, 0.0);
  BOOST_TEST(half > all_memory * 0.4);
  BOOST_TEST(half < all_memory * 0.6);

  const int64_t max_1gb = GetAvailableMemory(1.0, 1.0);
  BOOST_TEST(max_1gb <= all_memory);
  BOOST_TEST(max_1gb <= gb);

  const int64_t max_1gb_half = GetAvailableMemory(0.5, 1.0);
  BOOST_TEST(max_1gb_half < all_memory * 0.6);
  BOOST_TEST(max_1gb_half <= gb);
}

BOOST_AUTO_TEST_SUITE_END()
