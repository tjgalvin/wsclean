#include "../../main/commandline.h"
#include "../../main/wsclean.h"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(basic_integration)

BOOST_AUTO_TEST_CASE(version) {
  // The logger is set to quiet to suppress output by WSClean. Otherwise,
  // WSClean would generate a lot of output to stdout during the tests.
  Logger::SetVerbosity(Logger::QuietVerbosity);
  WSClean wsclean;
  CommandLine commandLine;
  std::vector<const char*> args = {"wsclean", "-version"};
  BOOST_CHECK(!commandLine.Parse(wsclean, args.size(), args.data(), false));
  // Restore logger -- not strictly necessary (since there should be no
  // output), but makes sure there is no side effect in case the test succeeds.
  Logger::SetVerbosity(Logger::NormalVerbosity);
}

BOOST_AUTO_TEST_CASE(bad_parameter) {
  Logger::SetVerbosity(Logger::QuietVerbosity);
  WSClean wsclean;
  CommandLine commandLine;
  std::vector<const char*> args = {"wsclean", "-this-is-not-a-valid-parameter"};
  BOOST_CHECK_THROW(commandLine.Parse(wsclean, args.size(), args.data(), false),
                    std::exception);
  Logger::SetVerbosity(Logger::NormalVerbosity);
}

BOOST_AUTO_TEST_SUITE_END()
