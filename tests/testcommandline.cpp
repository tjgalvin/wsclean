#include <boost/test/unit_test.hpp>

#include "../main/commandline.h"
#include "../main/wsclean.h"
#include "../structures/primarybeam.h"

namespace {
const char* kMWA_MS = "test_data/MWA_MOCK.ms/";
const char* kMWA_BDA_MS = "test_data/MWA_BDA_MOCK.ms/";
const char* kFacetFile2Facets = FACET_DEFINITION_FILE_2FACETS;
const char* kFacetFile4Facets = FACET_DEFINITION_FILE_4FACETS;
const char* kFacetSolutionFile = "test_data/mock_soltab_2pol.h5";

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

BOOST_AUTO_TEST_CASE(h5parm_basics) {
  BOOST_REQUIRE(boost::filesystem::is_directory(kMWA_MS));
  BOOST_REQUIRE(boost::filesystem::is_regular_file(kFacetFile2Facets));
  BOOST_REQUIRE(boost::filesystem::is_regular_file(kFacetFile4Facets));
  BOOST_REQUIRE(boost::filesystem::is_regular_file(kFacetSolutionFile));
}

BOOST_AUTO_TEST_CASE(h5parm_one_file) {
  CommandLine commandLine;
  WSClean wsclean;
  Settings settings;

  // One solution file, one MSet
  std::vector<const char*> args = baseArgs();
  args.push_back("-facet-regions");
  args.push_back(kFacetFile4Facets);
  args.push_back("-apply-facet-solutions");
  args.push_back(kFacetSolutionFile);
  args.push_back("amplitude000,phase000");
  args.push_back(kMWA_MS);
  commandLine.Parse(wsclean, args.size(), args.data(), false);
  settings = wsclean.GetSettings();
  BOOST_CHECK_EQUAL(settings.facetSolutionFiles.size(), 1u);
  BOOST_CHECK_EQUAL(settings.facetSolutionFiles[0], kFacetSolutionFile);

  // One solution file, two MSets
  wsclean.ResetSettings();
  args.push_back(kMWA_MS);
  commandLine.Parse(wsclean, args.size(), args.data(), false);
  settings = wsclean.GetSettings();
  BOOST_CHECK_EQUAL(settings.filenames.size(), 2u);
  BOOST_CHECK_EQUAL(settings.facetSolutionFiles.size(), 1u);
  BOOST_CHECK_EQUAL(settings.facetSolutionFiles[0], kFacetSolutionFile);
}

BOOST_AUTO_TEST_CASE(h5parm_multiple_files) {
  CommandLine commandLine;
  WSClean wsclean;
  Settings settings;

  // Two solution files, one MSet
  std::vector<const char*> args = baseArgs();
  args.push_back("-facet-regions");
  args.push_back(kFacetFile4Facets);
  args.push_back("-apply-facet-solutions");
  const std::string facet_solutions_str =
      std::string(kFacetSolutionFile) + "," + std::string(kFacetSolutionFile);
  const char* facet_solutions_char = facet_solutions_str.c_str();
  args.push_back(facet_solutions_char);
  args.push_back("amplitude000,phase000");
  args.push_back(kMWA_MS);
  BOOST_CHECK_THROW(commandLine.Parse(wsclean, args.size(), args.data(), false),
                    std::runtime_error);

  // Two solution files, two MSets
  wsclean.ResetSettings();
  args.push_back(kMWA_MS);
  commandLine.Parse(wsclean, args.size(), args.data(), false);
  settings = wsclean.GetSettings();
  BOOST_CHECK_EQUAL(settings.filenames.size(), 2u);
  BOOST_CHECK_EQUAL(settings.facetSolutionFiles.size(), 2u);
  BOOST_CHECK_EQUAL(settings.facetSolutionFiles[0], kFacetSolutionFile);
  BOOST_CHECK_EQUAL(settings.facetSolutionFiles[1], kFacetSolutionFile);
}

BOOST_AUTO_TEST_CASE(h5parm_inconsistent) {
  WSClean wsclean;
  CommandLine commandLine;

  std::vector<const char*> args = baseArgs();
  args.push_back("-facet-regions");
  args.push_back(kFacetFile2Facets);
  args.push_back("-apply-facet-solutions");
  args.push_back(kFacetSolutionFile);
  args.push_back("amplitude000,phase000");
  args.push_back(kMWA_MS);
  BOOST_CHECK_THROW(commandLine.Parse(wsclean, args.size(), args.data(), false),
                    std::runtime_error);
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
