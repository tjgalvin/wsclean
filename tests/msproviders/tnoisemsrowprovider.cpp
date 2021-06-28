#include "../../msproviders/noisemsrowprovider.h"
#include "../../io/logger.h"

#include <boost/test/unit_test.hpp>

#include <sstream>

BOOST_AUTO_TEST_SUITE(noise_ms_row_provider)

BOOST_AUTO_TEST_CASE(noise_baseline_map) {
  Logger::SetVerbosity(Logger::QuietVerbosity);
  BOOST_CHECK(NoiseMSRowProvider::NoiseMap().Empty());

  const std::string input =
      "0\t0\tnan\n"
      "0\t1\t-3.14\n"
      "0\t2\t42\n"
      "1\t1\t-nan\n"
      "2\t1\t0.987654321\n"
      "2\t2\t0\n";
  std::istringstream str(input);

  NoiseMSRowProvider::NoiseMap map(str);
  BOOST_CHECK(!map.Empty());
  BOOST_CHECK(!std::isfinite(map.GetNoiseValue(0, 0)));
  BOOST_CHECK_CLOSE_FRACTION(map.GetNoiseValue(0, 1), -3.14, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(map.GetNoiseValue(1, 0), -3.14, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(map.GetNoiseValue(0, 2), 42, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(map.GetNoiseValue(2, 0), 42, 1e-5);
  BOOST_CHECK(!std::isfinite(map.GetNoiseValue(1, 1)));
  BOOST_CHECK_CLOSE_FRACTION(map.GetNoiseValue(2, 1), 0.987654321, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(map.GetNoiseValue(1, 2), 0.987654321, 1e-5);
  BOOST_CHECK_EQUAL(map.GetNoiseValue(2, 2), 0);
  BOOST_CHECK_THROW(map.GetNoiseValue(0, 3), std::runtime_error);
  BOOST_CHECK_THROW(map.GetNoiseValue(3, 3), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
