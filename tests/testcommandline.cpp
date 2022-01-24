#include <boost/test/unit_test.hpp>

#include "../main/commandline.h"
#include "../main/wsclean.h"
#include "../structures/primarybeam.h"

namespace {
const char* kMWA_MS = "test_data/MWA_MOCK.ms/";
const char* kMWA_BDA_MS = "test_data/MWA_BDA_MOCK.ms/";

std::vector<const char*> baseArgs() {
  return {"wsclean", "-quiet", "-size",           "1024",   "512",
          "-scale",  "1amin",  "-multiscale",     "-niter", "1000000",
          "-mgain",  "0.8",    "-auto-threshold", "1",      "-auto-mask",
          "4"};
}

std::vector<const char*> baseArgsBda() {
  return {"wsclean", "-quiet", "-size",    "1024",      "1024",
          "-scale",  "1amin",  "-use-idg", "-idg-mode", "cpu"};
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

BOOST_AUTO_TEST_CASE(h5parm) {
  BOOST_REQUIRE(boost::filesystem::is_directory(kMWA_MS));

  CommandLine commandLine;

  // One MSet, one h5parm
  WSClean wsclean;
  std::vector<const char*> args1 = baseArgs();
  args1.push_back("-facet-regions");
  args1.push_back("dummy.reg");
  args1.push_back("-apply-facet-solutions");
  args1.push_back("dummy1.h5");
  args1.push_back("amplitude000,phase000");
  args1.push_back(kMWA_MS);
  commandLine.Parse(wsclean, args1.size(), args1.data(), false);
  const Settings settings1 = wsclean.GetSettings();
  BOOST_CHECK_EQUAL(settings1.facetSolutionFiles.size(), 1u);
  BOOST_CHECK_EQUAL(settings1.facetSolutionFiles[0], std::string("dummy1.h5"));

  // Two MSets, one h5parm
  wsclean.ResetSettings();
  args1.push_back(kMWA_MS);
  commandLine.Parse(wsclean, args1.size(), args1.data(), false);
  const Settings settings2 = wsclean.GetSettings();
  BOOST_CHECK_EQUAL(settings2.filenames.size(), 2u);
  BOOST_CHECK_EQUAL(settings2.facetSolutionFiles.size(), 1u);
  BOOST_CHECK_EQUAL(settings2.facetSolutionFiles[0], std::string("dummy1.h5"));

  // One MSet, two h5parm files
  wsclean.ResetSettings();
  std::vector<const char*> args3 = baseArgs();
  args3.push_back("-facet-regions");
  args3.push_back("dummy.reg");
  args3.push_back("-apply-facet-solutions");
  args3.push_back("dummy1.h5,dummy2.h5");
  args3.push_back("amplitude000,phase000");
  args3.push_back(kMWA_MS);
  BOOST_CHECK_THROW(
      commandLine.Parse(wsclean, args3.size(), args3.data(), false),
      std::runtime_error);

  // Two MSets, two h5parm files
  wsclean.ResetSettings();
  args3.push_back(kMWA_MS);
  commandLine.Parse(wsclean, args3.size(), args3.data(), false);
  const Settings settings3 = wsclean.GetSettings();
  BOOST_CHECK_EQUAL(settings3.filenames.size(), 2u);
  BOOST_CHECK_EQUAL(settings3.facetSolutionFiles.size(), 2u);
  BOOST_CHECK_EQUAL(settings3.facetSolutionFiles[0], std::string("dummy1.h5"));
  BOOST_CHECK_EQUAL(settings3.facetSolutionFiles[1], std::string("dummy2.h5"));
}

BOOST_AUTO_TEST_CASE(idg_bda_averaging) {
  BOOST_REQUIRE(boost::filesystem::is_directory(kMWA_MS));

  WSClean wsclean;
  CommandLine commandLine;

  std::vector<const char*> args = baseArgsBda();

  args.push_back("-baseline-averaging");
  args.push_back("10.0");
  args.push_back(kMWA_MS);

  BOOST_CHECK_THROW(commandLine.Parse(wsclean, args.size(), args.data(), false),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(idg_bda_mset) {
  BOOST_REQUIRE(boost::filesystem::is_directory(kMWA_BDA_MS));

  WSClean wsclean;
  CommandLine commandLine;

  std::vector<const char*> args = baseArgsBda();
  args.push_back(kMWA_BDA_MS);

  BOOST_CHECK_THROW(commandLine.Parse(wsclean, args.size(), args.data(), false),
                    std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
