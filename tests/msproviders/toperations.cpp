#include "../../msproviders/operations.h"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(ms_provider_operations)

namespace {
constexpr const char* kMsPathWithoutBda = "test_data/MWA_MOCK.ms";
constexpr const char* kTimeBdaMsPath = "test_data/MWA_BDA_MOCK.ms";
}  // namespace

BOOST_AUTO_TEST_CASE(has_frequency_bda) {
  casacore::MeasurementSet ms_without_bda(kMsPathWithoutBda);
  BOOST_CHECK(!wsclean::HasFrequencyBda(ms_without_bda));

  casacore::MeasurementSet time_bda_ms(kTimeBdaMsPath);
  BOOST_CHECK(!wsclean::HasFrequencyBda(time_bda_ms));
}

BOOST_AUTO_TEST_CASE(has_time_bda) {
  casacore::MeasurementSet ms_without_bda(kMsPathWithoutBda);
  BOOST_CHECK(!wsclean::HasTimeBda(ms_without_bda));

  casacore::MeasurementSet time_bda_ms(kTimeBdaMsPath);
  BOOST_CHECK(wsclean::HasTimeBda(time_bda_ms));
}

BOOST_AUTO_TEST_CASE(get_ms_interval) {
  casacore::MeasurementSet ms_without_bda(kMsPathWithoutBda);
  const double interval_without = wsclean::GetMsInterval(ms_without_bda);
  BOOST_CHECK_CLOSE_FRACTION(interval_without, 4.0, 1e-4);

  casacore::MeasurementSet time_bda_ms(kTimeBdaMsPath);
  const double interval_with = wsclean::GetMsInterval(time_bda_ms);
  BOOST_CHECK_CLOSE_FRACTION(interval_with, 8.0, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()
