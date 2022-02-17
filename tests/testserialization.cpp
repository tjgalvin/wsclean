#include <boost/test/unit_test.hpp>

#include "../scheduling/griddingtask.h"
#include "../idg/averagebeam.h"

#include <aocommon/image.h>
#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

using aocommon::SerialIStream;
using aocommon::SerialOStream;

BOOST_AUTO_TEST_SUITE(serialization)

BOOST_AUTO_TEST_CASE(empty_gridding_task) {
  GriddingTask a, b;
  a.subtractModel = true;
  b.subtractModel = false;

  SerialOStream ostr;
  a.Serialize(ostr);
  BOOST_CHECK_NE(ostr.size(), 0u);

  SerialIStream istr(std::move(ostr));
  b.Unserialize(istr);
  BOOST_CHECK_EQUAL(a.subtractModel, b.subtractModel);
}

BOOST_AUTO_TEST_CASE(image_weights) {
  ImageWeights weightsA(WeightMode::Briggs(0.5), 1024, 2048, 0.01, 0.01, false,
                        1.0);
  weightsA.SetAllValues(3.14);

  SerialOStream ostr;
  weightsA.Serialize(ostr);
  BOOST_CHECK_NE(ostr.size(), 0u);

  ImageWeights weightsB;
  SerialIStream istr(std::move(ostr));
  weightsB.Unserialize(istr);
  BOOST_CHECK_EQUAL(weightsA.Width(), weightsB.Width());
  BOOST_CHECK_EQUAL(weightsA.Height(), weightsB.Height());
  BOOST_CHECK(weightsA.GetWeightMode() == weightsB.GetWeightMode());

  std::vector<double> gridA(weightsA.Width() * weightsA.Height()),
      gridB(weightsB.Width() * weightsB.Height());
  weightsA.GetGrid(gridA.data());
  weightsB.GetGrid(gridB.data());
  BOOST_CHECK_EQUAL_COLLECTIONS(gridA.begin(), gridA.end(), gridB.begin(),
                                gridB.end());
}

BOOST_AUTO_TEST_CASE(msselection) {
  MSSelection a;
  a.SetBandId(3);
  a.SetChannelRange(4, 5);
  a.SetEvenOrOddTimesteps(MSSelection::EvenTimesteps);
  a.SetFieldIds(std::vector<size_t>{6, 7});
  a.SetInterval(8, 9);
  a.SetMaxUVWInM(11);
  a.SetMinUVWInM(10);

  SerialOStream ostr;
  a.Serialize(ostr);
  BOOST_CHECK_NE(ostr.size(), 0u);

  MSSelection b;
  SerialIStream istr(std::move(ostr));
  b.Unserialize(istr);
  BOOST_CHECK_EQUAL_COLLECTIONS(a.FieldIds().begin(), a.FieldIds().end(),
                                b.FieldIds().begin(), b.FieldIds().end());
  BOOST_CHECK_EQUAL(a.BandId(), b.BandId());
  BOOST_CHECK_EQUAL(a.ChannelRangeStart(), b.ChannelRangeStart());
  BOOST_CHECK_EQUAL(a.ChannelRangeEnd(), b.ChannelRangeEnd());
  BOOST_CHECK_EQUAL(a.EvenOrOddTimesteps(), b.EvenOrOddTimesteps());
  BOOST_CHECK_EQUAL(a.IntervalStart(), b.IntervalStart());
  BOOST_CHECK_EQUAL(a.IntervalEnd(), b.IntervalEnd());
  BOOST_CHECK_EQUAL(a.MaxUVWInM(), b.MaxUVWInM());
  BOOST_CHECK_EQUAL(a.MinUVWInM(), b.MinUVWInM());
}

BOOST_AUTO_TEST_CASE(image) {
  aocommon::Image a(12, 13);
  for (size_t i = 0; i != 12 * 13; ++i) a[i] = i + 1;

  SerialOStream ostr;
  a.Serialize(ostr);
  BOOST_CHECK_NE(ostr.size(), 0u);

  aocommon::Image b;
  SerialIStream istr(std::move(ostr));
  b.Unserialize(istr);
  BOOST_CHECK_EQUAL(a.Width(), b.Width());
  BOOST_CHECK_EQUAL(a.Height(), b.Height());
  BOOST_CHECK_EQUAL_COLLECTIONS(a.begin(), a.end(), b.begin(), b.end());
}

BOOST_AUTO_TEST_SUITE_END()
