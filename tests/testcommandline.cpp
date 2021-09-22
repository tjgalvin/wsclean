#include <boost/test/unit_test.hpp>

#include "../main/commandline.h"
#include "../main/wsclean.h"
#include "../structures/primarybeam.h"

namespace {
const char* kMWA_MS = "test_data/MWA_MOCK.ms/";

std::vector<const char*> baseArgs() {
  return {"wsclean", "-size",           "1024",   "512",        "-scale",
          "1amin",   "-multiscale",     "-niter", "1000000",    "-mgain",
          "0.8",     "-auto-threshold", "1",      "-auto-mask", "4"};
}
}  // namespace

/**
 * @brief Collection of tests to check the parsing of
 * command line input.
 */
BOOST_AUTO_TEST_SUITE(commandline)

BOOST_AUTO_TEST_CASE(pb_grid_size) {
  BOOST_REQUIRE(boost::filesystem::is_directory(kMWA_MS));

  WSClean wsclean;
  CommandLine commandLine;
  std::vector<const char*> args1 = baseArgs();
  args1.push_back("-apply-primary-beam");
  args1.push_back(kMWA_MS);
  commandLine.Parse(wsclean, args1.size(), args1.data(), false);
  const Settings settings1 = wsclean.GetSettings();
  BOOST_CHECK_EQUAL(settings1.primaryBeamGridSize, 32u);
  PrimaryBeam primaryBeam1(settings1);
  BOOST_CHECK_EQUAL(primaryBeam1.GetUndersamplingFactor(), 16u);
  BOOST_CHECK_EQUAL(primaryBeam1.GetBeamUpdateTime(), 1800u);

  std::vector<const char*> args2 = baseArgs();
  args2.push_back("-apply-primary-beam");
  args2.push_back("-pb-grid-size");
  args2.push_back("64");
  args2.push_back("-beam-aterm-update");
  args2.push_back("900");
  args2.push_back(kMWA_MS);
  commandLine.Parse(wsclean, args2.size(), args2.data(), false);
  const Settings settings2 = wsclean.GetSettings();
  BOOST_CHECK_EQUAL(settings2.primaryBeamGridSize, 64u);
  PrimaryBeam primaryBeam2(settings2);
  BOOST_CHECK_EQUAL(primaryBeam2.GetUndersamplingFactor(), 8u);
  BOOST_CHECK_EQUAL(primaryBeam2.GetBeamUpdateTime(), 900u);
}

BOOST_AUTO_TEST_SUITE_END()