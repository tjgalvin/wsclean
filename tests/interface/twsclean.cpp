#include <complex>
#include <filesystem>
#include <random>

#include <boost/test/unit_test.hpp>

#include <aocommon/fits/fitsreader.h>
#include <aocommon/image.h>

#include "../../interface/wsclean.h"

namespace wsclean {
namespace {

constexpr float kFluxDensity = 13.0;

InMemoryMs MakeInMemorySet() {
  InMemoryMs ms;
  ms.antenna_names = {"Dish1", "Dish2", "Dish3"};
  const std::vector<aocommon::ChannelInfo> kChannels{aocommon::ChannelInfo(120e6, 5e6),
                                                     aocommon::ChannelInfo(130e6, 5e6),
                                                     aocommon::ChannelInfo(140e6, 5e6)};
  aocommon::BandData band(kChannels, 130e6);
  ms.bands.AddBand(band);
  ms.field_name = "Test field";
  ms.has_frequency_bda = false;
  ms.has_time_bda = false;
  ms.interval = 42.0;
  ms.observer = "André";
  ms.phase_centre_dec = M_PI * 0.5;
  ms.phase_centre_ra = 0.0;
  ms.polarizations = {aocommon::PolarizationEnum::XX, aocommon::PolarizationEnum::YY};
  ms.start_time = 61041 * 24.0 * 60.0 * 60.0;  // 1 Jan 2026
  ms.telescope_name = "Test telescope";

  std::mt19937 rng;
  std::uniform_int_distribution<size_t> antenna_distribution(0, 2);
  std::uniform_real_distribution<double> uvw_distribution(-1000.0, 1000.0);
  const size_t kNPolarizations = ms.polarizations.size();

  for (size_t row_index = 0; row_index != 1000; ++row_index) {
    InMemoryRow& row = ms.rows.emplace_back();
    row.antenna1 = antenna_distribution(rng);
    do {
      row.antenna2 = antenna_distribution(rng);
    } while (row.antenna1 == row.antenna2);
    row.data.assign(kChannels.size() * kNPolarizations, kFluxDensity);
    row.data_desc_id = 0;
    row.field_id = 0;
    row.time = ms.start_time;
    row.uvw = {uvw_distribution(rng), uvw_distribution(rng), uvw_distribution(rng)};
    row.weights.assign(kChannels.size() * kNPolarizations, 1.0f);
  }
  return ms;
}

void CheckMetaData(const aocommon::FitsReader& reader) {
  BOOST_REQUIRE_EQUAL(reader.ImageWidth(), 1024);
  BOOST_REQUIRE_EQUAL(reader.ImageHeight(), 1024);
  BOOST_CHECK_CLOSE_FRACTION(reader.PixelSizeX() * 180.0 * 60.0 / M_PI, 1.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(reader.PixelSizeY() * 180.0 * 60.0 / M_PI, 1.0, 1e-6);
}

void CheckCenterPixel(aocommon::FitsReader& reader, float expected_value) {
  aocommon::Image image(reader.ImageWidth(), reader.ImageHeight());
  reader.Read(image.Data());
  float center_pixel = image.Value(reader.ImageWidth() / 2, reader.ImageHeight() / 2);
  BOOST_CHECK_CLOSE_FRACTION(center_pixel, expected_value, 1e-4);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(interface)

BOOST_AUTO_TEST_CASE(image) {
  std::vector<InMemoryMs> data = {MakeInMemorySet()};
  const std::string prefix_name = "interface-test";
  const std::string cmd_line = "-name " + prefix_name +
                               " -size 1024 1024 -quiet -scale 1amin -niter "
                               "1000 -mgain 0.8 -auto-threshold 5";
  Image(cmd_line, std::move(data));

  aocommon::FitsReader reader(prefix_name + "-image.fits");
  CheckMetaData(reader);
  CheckCenterPixel(reader, kFluxDensity);

  reader = aocommon::FitsReader(prefix_name + "-model.fits");
  CheckMetaData(reader);
  CheckCenterPixel(reader, kFluxDensity);

  reader = aocommon::FitsReader(prefix_name + "-psf.fits");
  CheckMetaData(reader);
  CheckCenterPixel(reader, 1.0f);

  reader = aocommon::FitsReader(prefix_name + "-residual.fits");
  CheckMetaData(reader);
  aocommon::Image image(reader.ImageWidth(), reader.ImageHeight());
  reader.Read(image.Data());
  BOOST_CHECK_LT(image.RMS(), 1e-3);

  std::filesystem::remove(prefix_name + "-image.fits");
  std::filesystem::remove(prefix_name + "-model.fits");
  std::filesystem::remove(prefix_name + "-psf.fits");
  std::filesystem::remove(prefix_name + "-residual.fits");
}

BOOST_AUTO_TEST_CASE(overwrite_image_prefix_name) {
  std::vector<InMemoryMs> data = {MakeInMemorySet()};
  const std::string cmd_line =
      "-name unused-prefix -size 128 128 -quiet -scale 1amin "
      "-niter 10 -mgain 0.8 -auto-threshold 5";

  const std::string custom_prefix = "custom-prefix-test";
  Image(cmd_line, std::move(data), custom_prefix);

  // Verify that the FITS files are created with the expected prefix and thus deletable.
  BOOST_TEST(std::filesystem::remove(custom_prefix + "-image.fits"));
  BOOST_TEST(std::filesystem::remove(custom_prefix + "-model.fits"));
  BOOST_TEST(std::filesystem::remove(custom_prefix + "-psf.fits"));
  BOOST_TEST(std::filesystem::remove(custom_prefix + "-residual.fits"));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace wsclean
