#include <vector>

#include <boost/test/unit_test.hpp>

#include "../../main/mshelper.h"
#include "../../main/settings.h"

using aocommon::BandData;
using aocommon::ChannelInfo;
using aocommon::MultiBandData;
using schaapcommon::reordering::MSSelection;

namespace wsclean {

BOOST_AUTO_TEST_SUITE(ms_helper)

BOOST_AUTO_TEST_CASE(initialize_ms_list) {
  Settings settings;
  settings.filenames = {"dummy_a.ms", "dummy_b.ms"};
  settings.doReorder = false;

  MSSelection global_selection;
  std::vector<aocommon::MultiBandData> bands(2);

  const std::vector<ChannelInfo> channels1{ChannelInfo(145.0e6, 10.0e6),
                                           ChannelInfo(155.0e6, 10.0e6)};
  const BandData band1(channels1, 150.0e6);
  bands[0].AddBand(band1);

  const std::vector<ChannelInfo> channels2{ChannelInfo(245.0e6, 10.0e6),
                                           ChannelInfo(255.0e6, 10.0e6)};
  const BandData band2(channels2, 250.0e6);
  bands[1].AddBand(band2);

  MsHelper ms_helper(settings, global_selection, bands);

  ImagingTableEntry entry;
  entry.lowestFrequency = 200.0e6;
  entry.highestFrequency = 300.0e6;
  entry.polarization = aocommon::PolarizationEnum::StokesI;
  const std::vector<MsListItem> result = ms_helper.InitializeMsList(entry);
  // Since the entry only covers the second band, the result should
  // have only one item (instead of two), with ms_index 1.
  BOOST_REQUIRE_EQUAL(result.size(), 1);
  BOOST_CHECK_EQUAL(result[0].ms_index, 1);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace wsclean
