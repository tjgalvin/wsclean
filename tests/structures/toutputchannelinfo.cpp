#include "../../structures/outputchannelinfo.h"

#include <vector>

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(output_channel_info)

BOOST_AUTO_TEST_CASE(smallest_theoretic_beam_size) {
  OutputChannelInfo one;
  one.theoreticBeamSize = 1.0;
  OutputChannelInfo two;
  two.theoreticBeamSize = 2.0;
  OutputChannelInfo nan;
  nan.theoreticBeamSize = std::numeric_limits<double>::quiet_NaN();

  double result = SmallestTheoreticBeamSize({});
  BOOST_CHECK(std::isnan(result));

  result = SmallestTheoreticBeamSize({one});
  BOOST_CHECK_EQUAL(result, 1.0);

  result = SmallestTheoreticBeamSize({nan});
  BOOST_CHECK(std::isnan(result));

  result = SmallestTheoreticBeamSize({one, two});
  BOOST_CHECK_EQUAL(result, 1.0);

  result = SmallestTheoreticBeamSize({two, nan, one});
  BOOST_CHECK_EQUAL(result, 1.0);
}

BOOST_AUTO_TEST_SUITE_END()
