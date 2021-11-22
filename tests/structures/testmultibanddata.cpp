#include <boost/test/unit_test.hpp>

#include "../../structures/multibanddata.h"

#include <vector>

using aocommon::ChannelInfo;

BOOST_AUTO_TEST_SUITE(multi_band_data)

BOOST_AUTO_TEST_CASE(empty) {
  MultiBandData multiBand;
  BOOST_CHECK_EQUAL(multiBand.BandCount(), 0u);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.BandStart(), 0.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.BandEnd(), 0.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.Bandwidth(), 0.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.CentreFrequency(), 0.0, 1e-6);
  BOOST_CHECK_EQUAL(multiBand.DataDescCount(), 0u);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.HighestFrequency(), 0.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.LowestFrequency(), 0.0, 1e-6);

  MultiBandData moved(std::move(multiBand));
  BOOST_CHECK_EQUAL(moved.BandCount(), 0u);

  MultiBandData copied(multiBand);
  BOOST_CHECK_EQUAL(copied.BandCount(), 0u);
}

BOOST_AUTO_TEST_CASE(irregular_bands) {
  // Band 1 has (purposely):
  // - A higher frequency than Band 2
  // - Fewer channels than Band 2
  // - Have a different channel width
  // MultiBandData should be able to handle this.
  std::vector<ChannelInfo> channels1{ChannelInfo(180e6, 10e6),
                                     ChannelInfo(190e6, 10e6)},
      channels2{ChannelInfo(140e6, 5e6), ChannelInfo(145e6, 5e6),
                ChannelInfo(150e6, 5e6)};
  MultiBandData multiBand;
  const size_t dataDescId1 = multiBand.AddBand(BandData(channels1, 185e6));
  const size_t dataDescId2 = multiBand.AddBand(BandData(channels2, 145e6));

  BOOST_CHECK_EQUAL(multiBand.BandCount(), 2u);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.BandStart(), 137.5e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.BandEnd(), 195e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.Bandwidth(), 195e6 - 137.5e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.CentreFrequency(),
                             0.5 * (195e6 + 137.5e6), 1e-6);
  BOOST_CHECK_EQUAL(multiBand.DataDescCount(), 2u);
  BOOST_CHECK_EQUAL(multiBand.GetBandIndex(dataDescId1), 0u);
  BOOST_CHECK_EQUAL(multiBand.GetBandIndex(dataDescId2), 1u);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.HighestFrequency(), 190e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.LowestFrequency(), 140e6, 1e-6);
  BOOST_CHECK_EQUAL(multiBand[dataDescId1].ChannelCount(), 2u);
  BOOST_CHECK_EQUAL(multiBand[dataDescId2].ChannelCount(), 3u);

  MultiBandData partialBandA(multiBand, 1, 2);
  BOOST_CHECK_EQUAL(partialBandA.BandCount(), 2u);
  BOOST_CHECK_CLOSE_FRACTION(partialBandA.BandStart(), 142.5e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(partialBandA.BandEnd(), 195e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(partialBandA.Bandwidth(), 195e6 - 142.5e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(partialBandA.CentreFrequency(),
                             0.5 * (195e6 + 142.5e6), 1e-6);
  BOOST_CHECK_EQUAL(partialBandA.DataDescCount(), 2u);
  BOOST_CHECK_EQUAL(partialBandA.GetBandIndex(dataDescId1), 0u);
  BOOST_CHECK_EQUAL(partialBandA.GetBandIndex(dataDescId2), 1u);
  BOOST_CHECK_CLOSE_FRACTION(partialBandA.HighestFrequency(), 190e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(partialBandA.LowestFrequency(), 145e6, 1e-6);
  BOOST_CHECK_EQUAL(partialBandA[dataDescId1].ChannelCount(), 1u);
  BOOST_CHECK_CLOSE_FRACTION(partialBandA[dataDescId1].ChannelFrequency(0),
                             190e6, 1e-6);
  BOOST_CHECK_EQUAL(partialBandA[dataDescId2].ChannelCount(), 1u);
  BOOST_CHECK_CLOSE_FRACTION(partialBandA[dataDescId2].ChannelFrequency(0),
                             145e6, 1e-6);

  MultiBandData partialBandB(multiBand, 1, 3);
  BOOST_CHECK_EQUAL(partialBandB.BandCount(), 2u);
  BOOST_CHECK_CLOSE_FRACTION(partialBandB.BandStart(), 142.5e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(partialBandB.BandEnd(), 195e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(partialBandB.Bandwidth(), 195e6 - 142.5e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(partialBandB.CentreFrequency(),
                             0.5 * (195e6 + 142.5e6), 1e-6);
  BOOST_CHECK_EQUAL(partialBandB.DataDescCount(), 2u);
  BOOST_CHECK_EQUAL(partialBandB.GetBandIndex(dataDescId1), 0u);
  BOOST_CHECK_EQUAL(partialBandB.GetBandIndex(dataDescId2), 1u);
  BOOST_CHECK_CLOSE_FRACTION(partialBandB.HighestFrequency(), 190e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(partialBandB.LowestFrequency(), 145e6, 1e-6);
  BOOST_CHECK_EQUAL(partialBandB[dataDescId1].ChannelCount(), 1u);
  BOOST_CHECK_CLOSE_FRACTION(partialBandB[dataDescId1].ChannelFrequency(0),
                             190e6, 1e-6);
  BOOST_CHECK_EQUAL(partialBandB[dataDescId2].ChannelCount(), 2u);
  BOOST_CHECK_CLOSE_FRACTION(partialBandB[dataDescId2].ChannelFrequency(0),
                             145e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(partialBandB[dataDescId2].ChannelFrequency(1),
                             150e6, 1e-6);

  MultiBandData copy(multiBand);
  BOOST_CHECK_EQUAL(copy.BandCount(), 2u);

  MultiBandData moved(std::move(multiBand));
  BOOST_CHECK_EQUAL(moved.BandCount(), 2u);
}

BOOST_AUTO_TEST_SUITE_END()
