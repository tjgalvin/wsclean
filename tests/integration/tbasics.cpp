#include "../../main/commandline.h"
#include "../../main/wsclean.h"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(basic_integration)

BOOST_AUTO_TEST_CASE(version) {
  WSClean wsclean;
  CommandLine commandLine;
  // wsclean -version
  std::vector<const char*> args = {"wsclean", "-version"};
  BOOST_CHECK(!commandLine.Parse(wsclean, args.size(), args.data(), false));
}

BOOST_AUTO_TEST_CASE(bad_parameter) {
  WSClean wsclean;
  CommandLine commandLine;
  // wsclean -this-is-not-a-valid-parameter
  std::vector<const char*> args = {"wsclean", "-this-is-not-a-valid-parameter"};
  BOOST_CHECK_THROW(commandLine.Parse(wsclean, args.size(), args.data(), false),
                    std::exception);
}

BOOST_AUTO_TEST_SUITE_END()
