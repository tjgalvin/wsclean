#include "../../main/commandline.h"
#include "../../main/wsclean.h"

#include <boost/test/unit_test.hpp>

namespace {
const char* kMWA_MS = "test_data/MWA_MOCK.ms/";
}

BOOST_AUTO_TEST_SUITE(vela_deconvolution)

BOOST_AUTO_TEST_CASE(wstacking) {
  BOOST_REQUIRE(boost::filesystem::is_directory(kMWA_MS));

  WSClean wsclean;
  CommandLine commandLine;
  // wsclean -niter 1000000 -size 1024 1024 -scale 1amin -parallel-gridding 4
  // -parallel-deconvolution 512 -multiscale -auto-threshold 1 -auto-mask 4
  // -channels-out 8 -join-channels -mgain 0.8 -deconvolution-channels 4
  // -fit-spectral-pol 2 ../MWA_MOCK.ms/
  std::vector<const char*> args = {"wsclean",
                                   "-quiet",
                                   "-size",
                                   "1024",
                                   "1024",
                                   "-scale",
                                   "1amin",
                                   "-parallel-gridding",
                                   "4",
                                   "-multiscale",
                                   "-parallel-deconvolution",
                                   "512",
                                   "-niter",
                                   "1000000",
                                   "-mgain",
                                   "0.8",
                                   "-channels-out",
                                   "8",
                                   "-join-channels",
                                   "-deconvolution-channels",
                                   "4",
                                   "-fit-spectral-pol",
                                   "2",
                                   "-auto-threshold",
                                   "1",
                                   "-auto-mask",
                                   "4",
                                   "-name",
                                   "mwa_test_run",
                                   kMWA_MS};
  commandLine.Parse(wsclean, args.size(), args.data(), false);
  commandLine.Run(wsclean);

  const size_t kNChan = 8;
  std::vector<std::string> filesExpected;
  for (size_t ch = 0; ch != kNChan + 1; ++ch) {
    std::string chStr = (ch == kNChan) ? ("MFS") : ("000" + std::to_string(ch));
    filesExpected.emplace_back("mwa_test_run-" + chStr + "-dirty.fits");
    filesExpected.emplace_back("mwa_test_run-" + chStr + "-image.fits");
    filesExpected.emplace_back("mwa_test_run-" + chStr + "-model.fits");
    filesExpected.emplace_back("mwa_test_run-" + chStr + "-psf.fits");
    filesExpected.emplace_back("mwa_test_run-" + chStr + "-residual.fits");
  }

  for (const std::string& str : filesExpected)
    BOOST_CHECK(boost::filesystem::is_regular_file(str));

  FitsReader dirtyReader("mwa_test_run-MFS-dirty.fits");
  FitsReader residualReader("mwa_test_run-MFS-residual.fits");
  BOOST_CHECK_EQUAL(residualReader.ImageWidth(), 1024u);
  BOOST_CHECK_EQUAL(residualReader.ImageHeight(), 1024u);
  Image residual(residualReader.ImageWidth(), residualReader.ImageHeight()),
      dirty(residualReader.ImageWidth(), residualReader.ImageHeight());
  residualReader.Read(residual.data());
  dirtyReader.Read(dirty.data());

  BOOST_CHECK_LT(residual.RMS(), 0.15);  // RMS should be less than 150 mJy.
  BOOST_CHECK_LT(
      residual.RMS(),
      2 * dirty.RMS());  // RMS of residual should be a lot better than dirty.

  for (const std::string& str : filesExpected)
    BOOST_CHECK(boost::filesystem::remove(str));
}

BOOST_AUTO_TEST_CASE(wgridder) {
  BOOST_REQUIRE(boost::filesystem::is_directory(kMWA_MS));

  WSClean wsclean;
  CommandLine commandLine;
  // wsclean -niter 1000000 -size 1024 1024 -scale 1amin -use-wgridder
  // -parallel-gridding 4 -parallel-deconvolution 512 -multiscale
  // -auto-threshold 1 -auto-mask 4 -channels-out 8 -join-channels -mgain 0.8
  // -deconvolution-channels 4 -fit-spectral-pol 2 ../MWA_MOCK.ms/
  std::vector<const char*> args = {"wsclean",
                                   "-quiet",
                                   "-size",
                                   "1024",
                                   "1024",
                                   "-scale",
                                   "1amin",
                                   "-use-wgridder",
                                   "-parallel-gridding",
                                   "4",
                                   "-multiscale",
                                   "-parallel-deconvolution",
                                   "512",
                                   "-niter",
                                   "1000000",
                                   "-mgain",
                                   "0.8",
                                   "-channels-out",
                                   "8",
                                   "-join-channels",
                                   "-deconvolution-channels",
                                   "4",
                                   "-fit-spectral-pol",
                                   "2",
                                   "-auto-threshold",
                                   "1",
                                   "-auto-mask",
                                   "4",
                                   "-name",
                                   "mwa_test_run",
                                   kMWA_MS};
  commandLine.Parse(wsclean, args.size(), args.data(), false);
  commandLine.Run(wsclean);

  const size_t kNChan = 8;
  std::vector<std::string> filesExpected;
  for (size_t ch = 0; ch != kNChan + 1; ++ch) {
    std::string chStr = (ch == kNChan) ? ("MFS") : ("000" + std::to_string(ch));
    filesExpected.emplace_back("mwa_test_run-" + chStr + "-dirty.fits");
    filesExpected.emplace_back("mwa_test_run-" + chStr + "-image.fits");
    filesExpected.emplace_back("mwa_test_run-" + chStr + "-model.fits");
    filesExpected.emplace_back("mwa_test_run-" + chStr + "-psf.fits");
    filesExpected.emplace_back("mwa_test_run-" + chStr + "-residual.fits");
  }

  for (const std::string& str : filesExpected)
    BOOST_CHECK(boost::filesystem::is_regular_file(str));

  FitsReader dirtyReader("mwa_test_run-MFS-dirty.fits");
  FitsReader residualReader("mwa_test_run-MFS-residual.fits");
  BOOST_CHECK_EQUAL(residualReader.ImageWidth(), 1024u);
  BOOST_CHECK_EQUAL(residualReader.ImageHeight(), 1024u);
  Image residual(residualReader.ImageWidth(), residualReader.ImageHeight()),
      dirty(residualReader.ImageWidth(), residualReader.ImageHeight());
  residualReader.Read(residual.data());
  dirtyReader.Read(dirty.data());

  BOOST_CHECK_LT(residual.RMS(), 0.15);  // RMS should be less than 150 mJy.
  BOOST_CHECK_LT(
      residual.RMS(),
      2 * dirty.RMS());  // RMS of residual should be a lot better than dirty.

  for (const std::string& str : filesExpected)
    BOOST_CHECK(boost::filesystem::remove(str));
}

BOOST_AUTO_TEST_SUITE_END()
