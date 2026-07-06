
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "../../gridding/visibilitymodifier.h"

namespace wsclean {

BOOST_AUTO_TEST_SUITE(visibility_modifier)

namespace {
constexpr size_t kNChannels = 4;
constexpr size_t kNStations = 2;
constexpr size_t kAntenna1 = 0;
constexpr size_t kAntenna2 = 1;
constexpr size_t kMsIndex = 0;
constexpr size_t kNPolarizations = 2;
constexpr float kWeights[kNPolarizations * kNChannels] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
constexpr float kImageWeights[kNChannels] = {1.0, 1.0, 1.0, 1.0};
constexpr size_t kTimeOffset = 0;
const std::vector<aocommon::ChannelInfo> kChannels{
    aocommon::ChannelInfo(150e6, 1e6),
    aocommon::ChannelInfo(151e6, 1e6),
    aocommon::ChannelInfo(152e6, 1e6),
    aocommon::ChannelInfo(153e6, 1e6),
};
const aocommon::BandData kFullBand(kChannels, 150e6);

}  // namespace

BOOST_AUTO_TEST_CASE(get_solution_response_for_bda_row) {
  std::vector<aocommon::ChannelInfo> channels1{kChannels[0], kChannels[2], kChannels[3]};
  const aocommon::BandData band1(channels1, kChannels[0].Frequency());
  const aocommon::BandData band2 = kFullBand;
  std::vector<aocommon::ChannelInfo> channels3{kChannels[1], kChannels[3]};
  const aocommon::BandData band3(channels3, kChannels[0].Frequency());

  constexpr size_t kDataDescId1 = 0;
  constexpr size_t kDataDescId2 = 1;
  constexpr size_t kDataDescId3 = 2;

  aocommon::MultiBandData bands;
  bands.SetBand(kDataDescId1, band1);
  bands.SetBand(kDataDescId2, band2);
  bands.SetBand(kDataDescId3, band3);

  const aocommon::JaggedArray<uint32_t> interpolation_table = MakeChannelInterpolationTable(bands);

  BOOST_REQUIRE_EQUAL(interpolation_table[kDataDescId1].size(), 3);
  BOOST_CHECK_EQUAL(interpolation_table[kDataDescId1][0], 0);
  BOOST_CHECK_EQUAL(interpolation_table[kDataDescId1][1], 2);
  BOOST_CHECK_EQUAL(interpolation_table[kDataDescId1][2], 3);

  BOOST_REQUIRE_EQUAL(interpolation_table[kDataDescId2].size(), 4);
  BOOST_CHECK_EQUAL(interpolation_table[kDataDescId2][0], 0);
  BOOST_CHECK_EQUAL(interpolation_table[kDataDescId2][1], 1);
  BOOST_CHECK_EQUAL(interpolation_table[kDataDescId2][2], 2);
  BOOST_CHECK_EQUAL(interpolation_table[kDataDescId2][3], 3);

  BOOST_REQUIRE_EQUAL(interpolation_table[kDataDescId3].size(), 2);
  BOOST_CHECK_EQUAL(interpolation_table[kDataDescId3][0], 1);
  BOOST_CHECK_EQUAL(interpolation_table[kDataDescId3][1], 3);

  const std::vector<aocommon::MC2x2F> beam_response(kNStations * kNChannels,
                                                    aocommon::MC2x2F(1.0, 0.0, 0.0, 1.0));
  std::vector<std::complex<float>> gain_response;
  for (size_t antenna = 0; antenna != kNStations; ++antenna) {
    for (size_t ch = 0; ch != kNChannels; ++ch) {
      gain_response.push_back(ch + 2 + antenna * 3);
      gain_response.push_back(0.0);
      gain_response.push_back(0.0);
      gain_response.push_back(ch + 3 + antenna * 3);
    }
  }
  VisibilityModifier modifier;
  std::vector<schaapcommon::h5parm::GainType> gain_types{
      schaapcommon::h5parm::GainType::kFullJones};
  modifier.SetH5Parm({}, {}, {}, gain_types);
  modifier.InitializeMockResponse(kNStations, bands, beam_response, gain_response);

  const auto [gains_a1, gains_a2] = modifier.GetParmResponseForRow(
      kMsIndex, kDataDescId1, kNStations, kAntenna1, kAntenna2, kTimeOffset);

  // First antenna, channel 0
  BOOST_CHECK_CLOSE_FRACTION(gains_a1[0].real(), 2.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(gains_a1[3].real(), 3.0, 1e-5);
  // First antenna, channel 1 is interpolated from channel 2
  BOOST_CHECK_CLOSE_FRACTION(gains_a1[4].real(), 4.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(gains_a1[7].real(), 5.0, 1e-5);
  // First antenna, channel 2 is interpolated from channel 3
  BOOST_CHECK_CLOSE_FRACTION(gains_a1[8].real(), 5.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(gains_a1[11].real(), 6.0, 1e-5);

  // Second antenna, channel 0
  BOOST_CHECK_CLOSE_FRACTION(gains_a2[0].real(), 5.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(gains_a2[3].real(), 6.0, 1e-5);
  // Second antenna, channel 1 is interpolated from channel 2
  BOOST_CHECK_CLOSE_FRACTION(gains_a2[4].real(), 7.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(gains_a2[7].real(), 8.0, 1e-5);
  // Second antenna, channel 2 is interpolated from channel 3
  BOOST_CHECK_CLOSE_FRACTION(gains_a2[8].real(), 8.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(gains_a2[11].real(), 9.0, 1e-5);

  const auto [gains_b1, gains_b2] = modifier.GetParmResponseForRow(
      kMsIndex, kDataDescId2, kNStations, kAntenna1, kAntenna2, kTimeOffset);

  // kDataDescId2 is the band thas has most channels, and so the interpolation output is the same as
  // the input
  for (size_t ch = 0; ch != kNChannels; ++ch) {
    const size_t index = ch * 4;
    BOOST_CHECK_CLOSE_FRACTION(gains_b1[index + 0].real(), ch + 2 + 0 * 3, 1e-5);
    BOOST_CHECK_CLOSE_FRACTION(gains_b1[index + 3].real(), ch + 3 + 0 * 3, 1e-5);
    BOOST_CHECK_CLOSE_FRACTION(gains_b2[index + 0].real(), ch + 2 + 1 * 3, 1e-5);
    BOOST_CHECK_CLOSE_FRACTION(gains_b2[index + 3].real(), ch + 3 + 1 * 3, 1e-5);
  }

  const auto [gains_c1, gains_c2] = modifier.GetParmResponseForRow(
      kMsIndex, kDataDescId3, kNStations, kAntenna1, kAntenna2, kTimeOffset);

  // First antenna, channel 0 is interpolated from channel 1
  BOOST_CHECK_CLOSE_FRACTION(gains_c1[0].real(), 3.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(gains_c1[3].real(), 4.0, 1e-5);
  // First antenna, channel 1 is interpolated from channel 3
  BOOST_CHECK_CLOSE_FRACTION(gains_c1[4].real(), 5.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(gains_c1[7].real(), 6.0, 1e-5);

  // Second antenna, channel 0 is interpolated from channel 1
  BOOST_CHECK_CLOSE_FRACTION(gains_c2[0].real(), 6.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(gains_c2[3].real(), 7.0, 1e-5);
  // Second antenna, channel 1 is interpolated from channel 3
  BOOST_CHECK_CLOSE_FRACTION(gains_c2[4].real(), 8.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(gains_c2[7].real(), 9.0, 1e-5);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace wsclean
