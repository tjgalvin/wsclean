#include <boost/test/unit_test.hpp>

#include "../scheduling/griddingtask.h"
#include "../idg/averagebeam.h"

BOOST_AUTO_TEST_SUITE(serialization)

BOOST_AUTO_TEST_CASE( basic )
{
	std::ostringstream str;
	Serializable::SerializeToBool(str, true);
	Serializable::SerializeToDouble(str, 3.14);
	Serializable::SerializeToFloat(str, 1.5);
	Serializable::SerializeToLDouble(str, 2.71);
	Serializable::SerializeToString(str, "hi!");
	Serializable::SerializeToUInt16(str, 160);
	Serializable::SerializeToUInt32(str, 320);
	Serializable::SerializeToUInt64(str, 640);
	Serializable::SerializeToUInt8(str, 80);
	
	std::istringstream istr(str.str());
	
	BOOST_CHECK_EQUAL(Serializable::UnserializeBool(istr), true);
	BOOST_CHECK_EQUAL(Serializable::UnserializeDouble(istr), 3.14);
	BOOST_CHECK_EQUAL(Serializable::UnserializeFloat(istr), 1.5);
	BOOST_CHECK_EQUAL(Serializable::UnserializeLDouble(istr), 2.71);
	BOOST_CHECK_EQUAL(Serializable::UnserializeString(istr), "hi!");
	BOOST_CHECK_EQUAL(Serializable::UnserializeUInt16(istr), 160);
	BOOST_CHECK_EQUAL(Serializable::UnserializeUInt32(istr), 320);
	BOOST_CHECK_EQUAL(Serializable::UnserializeUInt64(istr), 640);
	BOOST_CHECK_EQUAL(Serializable::UnserializeUInt8(istr), 80);
}

BOOST_AUTO_TEST_CASE( empty_gridding_task )
{
	GriddingTask a, b;
	a.addToModel = true;
	b.addToModel = false;
	
	std::ostringstream ostr;
	a.Serialize(ostr);
	std::string s(ostr.str());
	BOOST_CHECK_NE(s.size(),0);
	
	std::istringstream istr(s);
	b.Unserialize(istr);
	BOOST_CHECK_EQUAL(a.addToModel, b.addToModel);
}

BOOST_AUTO_TEST_CASE( image_weights )
{
	ImageWeights weightsA(
		WeightMode::Briggs(0.5),
		1024, 2048,
		0.01, 0.01,
		false,
		1.0
	);
	weightsA.SetAllValues(3.14);

	std::ostringstream ostr;
	weightsA.Serialize(ostr);
	std::string s(ostr.str());
	BOOST_CHECK_NE(s.size(),0);
	
	ImageWeights weightsB;
	std::istringstream istr(s);
	weightsB.Unserialize(istr);
	BOOST_CHECK_EQUAL(weightsA.Width(), weightsB.Width());
	BOOST_CHECK_EQUAL(weightsA.Height(), weightsB.Height());
	BOOST_CHECK(weightsA.GetWeightMode() == weightsB.GetWeightMode());
	
	std::vector<double>
		gridA(weightsA.Width()*weightsA.Height()),
		gridB(weightsB.Width()*weightsB.Height());
	weightsA.GetGrid(gridA.data());
	weightsB.GetGrid(gridB.data());
	BOOST_CHECK_EQUAL_COLLECTIONS(gridA.begin(), gridA.end(), gridB.begin(), gridB.end());
}

BOOST_AUTO_TEST_CASE( msselection )
{
	MSSelection a;
	a.SetBandId(3);
	a.SetChannelRange(4, 5);
	a.SetEvenOrOddTimesteps(MSSelection::EvenTimesteps);
	a.SetFieldIds(std::vector<size_t>{6, 7});
	a.SetInterval(8, 9);
	a.SetMaxUVWInM(11);
	a.SetMinUVWInM(10);

	std::ostringstream ostr;
	a.Serialize(ostr);
	std::string s(ostr.str());
	BOOST_CHECK_NE(s.size(),0);
	
	MSSelection b;
	std::istringstream istr(s);
	b.Unserialize(istr);
	BOOST_CHECK_EQUAL_COLLECTIONS(a.FieldIds().begin(), a.FieldIds().end(), b.FieldIds().begin(), b.FieldIds().end());
	BOOST_CHECK_EQUAL(a.BandId(), b.BandId());
	BOOST_CHECK_EQUAL(a.ChannelRangeStart(), b.ChannelRangeStart());
	BOOST_CHECK_EQUAL(a.ChannelRangeEnd(), b.ChannelRangeEnd());
	BOOST_CHECK_EQUAL(a.EvenOrOddTimesteps(), b.EvenOrOddTimesteps());
	BOOST_CHECK_EQUAL(a.IntervalStart(), b.IntervalStart());
	BOOST_CHECK_EQUAL(a.IntervalEnd(), b.IntervalEnd());
	BOOST_CHECK_EQUAL(a.MaxUVWInM(), b.MaxUVWInM());
	BOOST_CHECK_EQUAL(a.MinUVWInM(), b.MinUVWInM());
}

BOOST_AUTO_TEST_CASE( image )
{
	Image a(12, 13);
	for(size_t i=0; i!=12*13; ++i)
		a[i] = i+1;

	std::ostringstream ostr;
	a.Serialize(ostr);
	std::string s(ostr.str());
	BOOST_CHECK_NE(s.size(), 0);
	
	Image b;
	std::istringstream istr(s);
	b.Unserialize(istr);
	BOOST_CHECK_EQUAL(a.Width(), b.Width());
	BOOST_CHECK_EQUAL(a.Height(), b.Height());
	BOOST_CHECK_EQUAL_COLLECTIONS(a.begin(), a.end(), b.begin(), b.end());
}

BOOST_AUTO_TEST_CASE( average_beam_empty )
{
	AverageBeam a, b;
	
	std::ostringstream ostr;
	a.Serialize(ostr);
	std::string s(ostr.str());
	BOOST_CHECK_NE(s.size(), 0);
	
	b.SetMatrixInverseBeam(std::shared_ptr<std::vector<std::complex<float>>>(new std::vector<std::complex<float>>(12, 3)));
	b.SetScalarBeam(std::shared_ptr<std::vector<float>>(new std::vector<float>(11, 4)));
	
	std::istringstream istr(s);
	b.Unserialize(istr);
	BOOST_CHECK(!b.MatrixInverseBeam());
	BOOST_CHECK(!b.ScalarBeam());
}

BOOST_AUTO_TEST_CASE( average_beam_filled )
{
	AverageBeam a, b;
	a.SetMatrixInverseBeam(std::shared_ptr<std::vector<std::complex<float>>>(new std::vector<std::complex<float>>(12, std::complex<float>(3, 0))));
	a.SetScalarBeam(std::shared_ptr<std::vector<float>>(new std::vector<float>(11, 4)));
	
	std::ostringstream ostr;
	a.Serialize(ostr);
	std::string s(ostr.str());
	BOOST_CHECK_NE(s.size(), 0);
	
	std::istringstream istr(s);
	b.Unserialize(istr);
	BOOST_CHECK_EQUAL(b.MatrixInverseBeam()->size(), 12);
	BOOST_CHECK_EQUAL(b.ScalarBeam()->size(), 11);
	BOOST_CHECK_EQUAL(b.MatrixInverseBeam()->at(8), std::complex<float>(3, 0));
	BOOST_CHECK_EQUAL(b.ScalarBeam()->at(7), 4);
}

BOOST_AUTO_TEST_SUITE_END()
