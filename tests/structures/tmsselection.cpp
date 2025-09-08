#include "../../structures/msselection.h"

#include <boost/test/unit_test.hpp>

#include <aocommon/multibanddata.h>
#include <schaapcommon/reordering/channelrange.h>

using aocommon::BandData;
using aocommon::ChannelInfo;
using aocommon::MultiBandData;

using schaapcommon::reordering::ChannelRange;

namespace wsclean {
namespace {
constexpr size_t kDataDescId = 2;
MultiBandData ExampleBand() {
  const std::vector<ChannelInfo> channels1{ChannelInfo(90.0e6, 10.0e6),
                                           ChannelInfo(100.0e6, 10.0e6),
                                           ChannelInfo(110.0e6, 10.0e6)};
  BandData band1(channels1, 100.0e6);
  MultiBandData bands;
  bands.SetBand(kDataDescId, band1);
  return bands;
}

MultiBandData InvertedBand() {
  const std::vector<ChannelInfo> channels1{ChannelInfo(110.0e6, 10.0e6),
                                           ChannelInfo(100.0e6, 10.0e6),
                                           ChannelInfo(90.0e6, 10.0e6)};
  BandData band1(channels1, 100.0e6);
  MultiBandData bands;
  bands.SetBand(kDataDescId, band1);
  return bands;
}
}  // namespace

BOOST_AUTO_TEST_SUITE(ms_selection)

BOOST_AUTO_TEST_CASE(get_channel_range_empty) {
  const ChannelRange range =
      SelectMsChannels(ExampleBand(), kDataDescId, 102.0e6, 103.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK(range.Empty());
}

BOOST_AUTO_TEST_CASE(get_channel_range_low_range1) {
  const ChannelRange range =
      SelectMsChannels(ExampleBand(), kDataDescId, 80.0e6, 95.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK_EQUAL(range.start, 0);
  BOOST_CHECK_EQUAL(range.end, 1);
}

BOOST_AUTO_TEST_CASE(get_channel_range_low_range2) {
  const ChannelRange range =
      SelectMsChannels(ExampleBand(), kDataDescId, 90.0e6, 95.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK_EQUAL(range.start, 0);
  BOOST_CHECK_EQUAL(range.end, 1);
}

BOOST_AUTO_TEST_CASE(get_channel_range_low_range3) {
  const ChannelRange range =
      SelectMsChannels(ExampleBand(), kDataDescId, 90.0e6, 90.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK_EQUAL(range.start, 0);
  BOOST_CHECK_EQUAL(range.end, 1);
}

BOOST_AUTO_TEST_CASE(get_channel_range_mid_range1) {
  const ChannelRange range =
      SelectMsChannels(ExampleBand(), kDataDescId, 95.0e6, 105.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK_EQUAL(range.start, 1);
  BOOST_CHECK_EQUAL(range.end, 2);
}

BOOST_AUTO_TEST_CASE(get_channel_range_mid_range2) {
  const ChannelRange range =
      SelectMsChannels(ExampleBand(), kDataDescId, 95.0e6, 100.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK_EQUAL(range.start, 1);
  BOOST_CHECK_EQUAL(range.end, 2);
}

BOOST_AUTO_TEST_CASE(get_channel_range_mid_range3) {
  const ChannelRange range =
      SelectMsChannels(ExampleBand(), kDataDescId, 100.0e6, 105.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK_EQUAL(range.start, 1);
  BOOST_CHECK_EQUAL(range.end, 2);
}

BOOST_AUTO_TEST_CASE(get_channel_range_high_range1) {
  const ChannelRange range =
      SelectMsChannels(ExampleBand(), kDataDescId, 105.0e6, 115.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK_EQUAL(range.start, 2);
  BOOST_CHECK_EQUAL(range.end, 3);
}

BOOST_AUTO_TEST_CASE(get_channel_range_high_range2) {
  const ChannelRange range =
      SelectMsChannels(ExampleBand(), kDataDescId, 110.0e6, 115.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK_EQUAL(range.start, 2);
  BOOST_CHECK_EQUAL(range.end, 3);
}

BOOST_AUTO_TEST_CASE(get_channel_range_high_range3) {
  const ChannelRange range =
      SelectMsChannels(ExampleBand(), kDataDescId, 105.0e6, 110.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK_EQUAL(range.start, 2);
  BOOST_CHECK_EQUAL(range.end, 3);
}

BOOST_AUTO_TEST_CASE(get_channel_range_full_range1) {
  const ChannelRange range =
      SelectMsChannels(ExampleBand(), kDataDescId, 90.0e6, 110.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK_EQUAL(range.start, 0);
  BOOST_CHECK_EQUAL(range.end, 3);
}

BOOST_AUTO_TEST_CASE(get_channel_range_full_range2) {
  const ChannelRange range =
      SelectMsChannels(ExampleBand(), kDataDescId, 0.0e6, 200.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK_EQUAL(range.start, 0);
  BOOST_CHECK_EQUAL(range.end, 3);
}

BOOST_AUTO_TEST_CASE(get_channel_range_inverted_low_range1) {
  const ChannelRange range =
      SelectMsChannels(InvertedBand(), kDataDescId, 80.0e6, 95.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK_EQUAL(range.start, 2);
  BOOST_CHECK_EQUAL(range.end, 3);
}

BOOST_AUTO_TEST_CASE(get_channel_range_inverted_low_range2) {
  const ChannelRange range =
      SelectMsChannels(InvertedBand(), kDataDescId, 90.0e6, 95.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK_EQUAL(range.start, 2);
  BOOST_CHECK_EQUAL(range.end, 3);
}

BOOST_AUTO_TEST_CASE(get_channel_range_inverted_low_range3) {
  const ChannelRange range =
      SelectMsChannels(InvertedBand(), kDataDescId, 90.0e6, 90.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK_EQUAL(range.start, 2);
  BOOST_CHECK_EQUAL(range.end, 3);
}

BOOST_AUTO_TEST_CASE(get_channel_range_inverted_mid_range1) {
  const ChannelRange range =
      SelectMsChannels(InvertedBand(), kDataDescId, 95.0e6, 105.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK_EQUAL(range.start, 1);
  BOOST_CHECK_EQUAL(range.end, 2);
}

BOOST_AUTO_TEST_CASE(get_channel_range_inverted_mid_range2) {
  const ChannelRange range =
      SelectMsChannels(InvertedBand(), kDataDescId, 95.0e6, 100.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK_EQUAL(range.start, 1);
  BOOST_CHECK_EQUAL(range.end, 2);
}

BOOST_AUTO_TEST_CASE(get_channel_range_inverted_mid_range3) {
  const ChannelRange range =
      SelectMsChannels(InvertedBand(), kDataDescId, 100.0e6, 105.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK_EQUAL(range.start, 1);
  BOOST_CHECK_EQUAL(range.end, 2);
}

BOOST_AUTO_TEST_CASE(get_channel_range_inverted_high_range1) {
  const ChannelRange range =
      SelectMsChannels(InvertedBand(), kDataDescId, 105.0e6, 115.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK_EQUAL(range.start, 0);
  BOOST_CHECK_EQUAL(range.end, 1);
}

BOOST_AUTO_TEST_CASE(get_channel_range_inverted_high_range2) {
  const ChannelRange range =
      SelectMsChannels(InvertedBand(), kDataDescId, 110.0e6, 115.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK_EQUAL(range.start, 0);
  BOOST_CHECK_EQUAL(range.end, 1);
}

BOOST_AUTO_TEST_CASE(get_channel_range_inverted_high_range3) {
  const ChannelRange range =
      SelectMsChannels(InvertedBand(), kDataDescId, 105.0e6, 110.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK_EQUAL(range.start, 0);
  BOOST_CHECK_EQUAL(range.end, 1);
}

BOOST_AUTO_TEST_CASE(get_channel_range_inverted_full_range1) {
  const ChannelRange range =
      SelectMsChannels(InvertedBand(), kDataDescId, 90.0e6, 110.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK_EQUAL(range.start, 0);
  BOOST_CHECK_EQUAL(range.end, 3);
}

BOOST_AUTO_TEST_CASE(get_channel_range_inverted_full_range2) {
  const ChannelRange range =
      SelectMsChannels(InvertedBand(), kDataDescId, 0.0e6, 200.0e6);
  BOOST_CHECK_EQUAL(range.data_desc_id, kDataDescId);
  BOOST_CHECK_EQUAL(range.start, 0);
  BOOST_CHECK_EQUAL(range.end, 3);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace wsclean
