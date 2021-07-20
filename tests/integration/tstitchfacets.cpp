#include "../../main/commandline.h"
#include "../../main/wsclean.h"

#include <boost/test/unit_test.hpp>

namespace {
const char* kMWA_MS = "test_data/MWA_MOCK.ms/";
const char* kFacets = "test_data/ds9_2facets.reg";
}  // namespace

BOOST_AUTO_TEST_SUITE(facet_stitching)

void CheckAndRemoveFileNames(const std::string& prefix, bool isIDG = false) {
  std::vector<std::string> filesExpected;
  if (isIDG) {
    filesExpected.emplace_back(prefix + "-dirty.fits");
    filesExpected.emplace_back(prefix + "-image.fits");
  } else {
    filesExpected.emplace_back(prefix + "-XX-dirty.fits");
    filesExpected.emplace_back(prefix + "-YY-dirty.fits");
    filesExpected.emplace_back(prefix + "-XX-image.fits");
    filesExpected.emplace_back(prefix + "-YY-image.fits");
  }

  for (const std::string& str : filesExpected) {
    BOOST_CHECK(boost::filesystem::is_regular_file(str));
    BOOST_CHECK(boost::filesystem::remove(str));
  }
}

BOOST_AUTO_TEST_CASE(wstacking) {
  BOOST_REQUIRE(boost::filesystem::is_directory(kMWA_MS));
  BOOST_REQUIRE(boost::filesystem::is_regular_file(kFacets));
  const char* prefix = "facet-stitch-wstack";

  WSClean wsclean;
  CommandLine commandLine;
  std::vector<const char*> args = {
      "wsclean", "-quiet", "-size", "256",   "256",
      "-scale",  "4amin",  "-pol",  "XX,YY", "-facet-regions",
      kFacets,   "-name",  prefix,  kMWA_MS};
  commandLine.Parse(wsclean, args.size(), args.data(), false);
  commandLine.Run(wsclean);

  CheckAndRemoveFileNames(std::string(prefix));
}

#ifdef HAVE_WGRIDDER
BOOST_AUTO_TEST_CASE(wgridder) {
  BOOST_REQUIRE(boost::filesystem::is_directory(kMWA_MS));
  BOOST_REQUIRE(boost::filesystem::is_regular_file(kFacets));
  const char* prefix = "facet-stitch-wgridder";

  WSClean wsclean;
  CommandLine commandLine;
  std::vector<const char*> args = {
      "wsclean", "-size",         "256",  "256",   "-scale",
      "4amin",   "-use-wgridder", "-pol", "XX,YY", "-facet-regions",
      kFacets,   "-name",         prefix, kMWA_MS};
  commandLine.Parse(wsclean, args.size(), args.data(), false);
  commandLine.Run(wsclean);

  CheckAndRemoveFileNames(std::string(prefix));
}
#endif

// Since CI does not compile IDG, this test won't be run as part
// of the CI/CD pipeline
#ifdef HAVE_IDG
BOOST_AUTO_TEST_CASE(idg) {
  BOOST_REQUIRE(boost::filesystem::is_directory(kMWA_MS));
  BOOST_REQUIRE(boost::filesystem::is_regular_file(kFacets));
  const char* prefix = "facet-stitch-idg";

  WSClean wsclean;
  CommandLine commandLine;
#ifdef HAVE_EVERYBEAM
  std::vector<const char*> args = {
      "wsclean",  "-quiet",         "-size", "256",   "256",  "-scale", "4amin",
      "-use-idg", "-facet-regions", kFacets, "-name", prefix, kMWA_MS};
#else
  std::vector<const char*> args = {
      "wsclean",  "-size",          "256",   "256",   "-scale", "4amin",
      "-use-idg", "-facet-regions", kFacets, "-name", prefix,   kMWA_MS};
#endif
  commandLine.Parse(wsclean, args.size(), args.data(), false);
  commandLine.Run(wsclean);

  CheckAndRemoveFileNames(std::string(prefix), true);
}
#endif

BOOST_AUTO_TEST_SUITE_END()
