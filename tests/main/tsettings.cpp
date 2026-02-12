#include <boost/test/unit_test.hpp>

#include "../../main/settings.h"

namespace wsclean {

namespace {
/// @return A minimal Settings object, which passes Validate().
Settings ValidSettings() {
  Settings settings;
  settings.trimmedImageWidth = 42;
  settings.trimmedImageHeight = 42;
  settings.pixelScaleX = 0.0042;
  settings.pixelScaleY = 0.0042;
  return settings;
}
}  // namespace

BOOST_AUTO_TEST_SUITE(settings)

BOOST_AUTO_TEST_CASE(get_feather_size_default) {
  Settings settings;
  settings.trimmedImageWidth = 10000;
  settings.trimmedImageHeight = 15625;
  BOOST_CHECK_EQUAL(settings.GetFeatherSize(),
                    125);  // 1% of sqrt(10000*15625)
}

BOOST_AUTO_TEST_CASE(get_feather_size_override) {
  Settings settings;
  settings.featherSize = 100;
  settings.trimmedImageWidth = 5000;
  settings.trimmedImageHeight = 7000;
  BOOST_CHECK_EQUAL(settings.GetFeatherSize(), 100);
}

BOOST_AUTO_TEST_CASE(validate_default_settings) {
  Settings settings = ValidSettings();
  BOOST_CHECK_NO_THROW(settings.Validate());
}

BOOST_AUTO_TEST_CASE(validate_mpi_settings) {
  {  // Basic test, with 1 MPI Node and 1 output channel.
    Settings settings = ValidSettings();
    settings.nMpiNodes = 1;
    settings.masterDoesWork = true;
    settings.channelToNode.push_back(0);
    BOOST_CHECK_NO_THROW(settings.Validate());
  }

  {  // Test with 3 MPI nodes, inactive master and 2 output channels.
    Settings settings = ValidSettings();
    settings.channelsOut = 2;
    settings.nMpiNodes = 3;
    settings.masterDoesWork = false;
    settings.channelToNode.push_back(1);
    settings.channelToNode.push_back(2);
    BOOST_CHECK_NO_THROW(settings.Validate());
  }

  {  // Test with inactive master and a single node.
    Settings settings = ValidSettings();
    settings.nMpiNodes = 1;
    settings.masterDoesWork = false;
    settings.channelToNode.push_back(0);
    BOOST_CHECK_THROW(settings.Validate(), std::runtime_error);
  }

  {
    Settings settings = ValidSettings();
    settings.nMpiNodes = 1;
    settings.masterDoesWork = true;
    settings.channelToNode.push_back(42);  // Invalid node index.
    BOOST_CHECK_THROW(settings.Validate(), std::runtime_error);
  }

  {
    Settings settings = ValidSettings();
    settings.nMpiNodes = 1;
    settings.masterDoesWork = true;
    settings.channelToNode.clear();  // Invalid size: Too small.
    BOOST_CHECK_THROW(settings.Validate(), std::runtime_error);
  }

  {
    Settings settings = ValidSettings();
    settings.nMpiNodes = 1;
    settings.masterDoesWork = true;
    settings.channelToNode.push_back(0);  // Invalid size: Too big.
    settings.channelToNode.push_back(0);
    BOOST_CHECK_THROW(settings.Validate(), std::runtime_error);
  }
}

BOOST_AUTO_TEST_CASE(invalid_mpi_message_size) {
  Settings settings = ValidSettings();
  settings.nMpiNodes = 1;  // Enables MPI, including message size check.
  settings.maxMpiMessageSize = 42'000'000'000;
  BOOST_CHECK_THROW(settings.Validate(), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(get_provider_polarization) {
  Settings settings = ValidSettings();
  settings.polarizations = {aocommon::PolarizationEnum::StokesU,
                            aocommon::PolarizationEnum::StokesU};
  settings.facetRegionFilename = "somefile.reg";
  settings.applyFacetBeam = true;
  BOOST_CHECK(
      settings.GetProviderPolarization(aocommon::PolarizationEnum::StokesU) ==
      aocommon::PolarizationEnum::Instrumental);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace wsclean
