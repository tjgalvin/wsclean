#include "../../structures/primarybeam.h"

#include <boost/test/unit_test.hpp>

#include "../../msproviders/contiguousms.h"

using schaapcommon::reordering::MSSelection;
using schaapcommon::reordering::StorageManagerType;

namespace wsclean {

BOOST_AUTO_TEST_SUITE(primary_beam)

BOOST_AUTO_TEST_CASE(get_beam_intervals) {
  const std::string ms_path = "test_data/MWA_MOCK.ms";
  ContiguousMS ms(ms_path, "DATA", "MODEL_DATA", StorageManagerType::Default,
                  MSSelection(), aocommon::PolarizationEnum::StokesI, 0, false);
  const std::vector<BeamInterval> intervals = GetBeamIntervals(ms, 1);
  BOOST_REQUIRE_EQUAL(intervals.size(), 1);
  BeamInterval first = intervals.front();
  BOOST_CHECK_EQUAL(first.start_row, 0);
  BOOST_CHECK_EQUAL(first.end_row, 2555);
  BOOST_CHECK_CLOSE_FRACTION(first.central_time, 4875418129.0, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace wsclean
