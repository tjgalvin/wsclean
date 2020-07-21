#include <boost/test/unit_test.hpp>

#include "../scheduling/griddingtask.h"
#include "../idg/averagebeam.h"

#include "../serialostream.h"
#include "../serialistream.h"

BOOST_AUTO_TEST_SUITE(serialization)

BOOST_AUTO_TEST_CASE(basic) {
  SerialOStream ostr;
  ostr.Bool(true)
      .UInt8(80)
      .UInt16(160)
      .UInt32(320)
      .UInt64(640)
      .Float(1.5)
      .Double(3.14)
      .LDouble(2.71)
      .String("hi!");

  SerialIStream istr(std::move(ostr));

  BOOST_CHECK_EQUAL(istr.Bool(), true);
  BOOST_CHECK_EQUAL(istr.UInt8(), 80);
  BOOST_CHECK_EQUAL(istr.UInt16(), 160);
  BOOST_CHECK_EQUAL(istr.UInt32(), 320);
  BOOST_CHECK_EQUAL(istr.UInt64(), 640);
  BOOST_CHECK_EQUAL(istr.Float(), 1.5);
  BOOST_CHECK_EQUAL(istr.Double(), 3.14);
  BOOST_CHECK_EQUAL(istr.LDouble(), 2.71);
  BOOST_CHECK_EQUAL(istr.String(), "hi!");
}

BOOST_AUTO_TEST_CASE(vector64) {
  std::vector<int32_t> int32vecA, int32vecB{12, 13, 14};
  std::vector<uint32_t> uint32vecA, uint32vecB{15, 16, 17};
  std::vector<uint64_t> int64vecA, int64vecB{18, 19, 20};

  SerialOStream ostr;
  ostr.VectorUInt64(int32vecA)
      .VectorUInt64(int32vecB)
      .VectorUInt64(uint32vecA)
      .VectorUInt64(uint32vecB)
      .VectorUInt64(int64vecA)
      .VectorUInt64(int64vecB);

  SerialIStream istr(std::move(ostr));

  std::vector<int32_t> out_int32vecA, out_int32vecB;
  std::vector<uint32_t> out_uint32vecA, out_uint32vecB;
  std::vector<uint64_t> out_int64vecA, out_int64vecB;
  istr.VectorUInt64(out_int32vecA)
      .VectorUInt64(out_int32vecB)
      .VectorUInt64(out_uint32vecA)
      .VectorUInt64(out_uint32vecB)
      .VectorUInt64(out_int64vecA)
      .VectorUInt64(out_int64vecB);
  BOOST_CHECK_EQUAL_COLLECTIONS(int32vecA.begin(), int32vecA.end(),
                                out_int32vecA.begin(), out_int32vecA.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(int32vecB.begin(), int32vecB.end(),
                                out_int32vecB.begin(), out_int32vecB.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(uint32vecA.begin(), uint32vecA.end(),
                                out_uint32vecA.begin(), out_uint32vecA.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(uint32vecB.begin(), uint32vecB.end(),
                                out_uint32vecB.begin(), out_uint32vecB.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(int64vecA.begin(), int64vecA.end(),
                                out_int64vecA.begin(), out_int64vecA.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(int64vecB.begin(), int64vecB.end(),
                                out_int64vecB.begin(), out_int64vecB.end());
}

BOOST_AUTO_TEST_CASE(empty_gridding_task) {
  GriddingTask a, b;
  a.addToModel = true;
  b.addToModel = false;

  SerialOStream ostr;
  a.Serialize(ostr);
  BOOST_CHECK_NE(ostr.size(), 0);

  SerialIStream istr(std::move(ostr));
  b.Unserialize(istr);
  BOOST_CHECK_EQUAL(a.addToModel, b.addToModel);
}

BOOST_AUTO_TEST_CASE(image_weights) {
  ImageWeights weightsA(WeightMode::Briggs(0.5), 1024, 2048, 0.01, 0.01, false,
                        1.0);
  weightsA.SetAllValues(3.14);

  SerialOStream ostr;
  weightsA.Serialize(ostr);
  BOOST_CHECK_NE(ostr.size(), 0);

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
  BOOST_CHECK_NE(ostr.size(), 0);

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
  Image a(12, 13);
  for (size_t i = 0; i != 12 * 13; ++i) a[i] = i + 1;

  SerialOStream ostr;
  a.Serialize(ostr);
  BOOST_CHECK_NE(ostr.size(), 0);

  Image b;
  SerialIStream istr(std::move(ostr));
  b.Unserialize(istr);
  BOOST_CHECK_EQUAL(a.Width(), b.Width());
  BOOST_CHECK_EQUAL(a.Height(), b.Height());
  BOOST_CHECK_EQUAL_COLLECTIONS(a.begin(), a.end(), b.begin(), b.end());
}

BOOST_AUTO_TEST_CASE(average_beam_empty) {
  AverageBeam a, b;

  SerialOStream ostr;
  a.Serialize(ostr);
  BOOST_CHECK_NE(ostr.size(), 0);

  b.SetMatrixInverseBeam(std::shared_ptr<std::vector<std::complex<float>>>(
      new std::vector<std::complex<float>>(12, 3)));
  b.SetScalarBeam(
      std::shared_ptr<std::vector<float>>(new std::vector<float>(11, 4)));

  SerialIStream istr(std::move(ostr));
  b.Unserialize(istr);
  BOOST_CHECK(!b.MatrixInverseBeam());
  BOOST_CHECK(!b.ScalarBeam());
}

BOOST_AUTO_TEST_CASE(average_beam_filled) {
  AverageBeam a, b;
  a.SetMatrixInverseBeam(std::shared_ptr<std::vector<std::complex<float>>>(
      new std::vector<std::complex<float>>(12, std::complex<float>(3, 0))));
  a.SetScalarBeam(
      std::shared_ptr<std::vector<float>>(new std::vector<float>(11, 4)));

  SerialOStream ostr;
  a.Serialize(ostr);
  BOOST_CHECK_NE(ostr.size(), 0);

  SerialIStream istr(std::move(ostr));
  b.Unserialize(istr);
  BOOST_CHECK_EQUAL(b.MatrixInverseBeam()->size(), 12);
  BOOST_CHECK_EQUAL(b.ScalarBeam()->size(), 11);
  BOOST_CHECK_EQUAL(b.MatrixInverseBeam()->at(8), std::complex<float>(3, 0));
  BOOST_CHECK_EQUAL(b.ScalarBeam()->at(7), 4);
}

BOOST_AUTO_TEST_SUITE_END()
