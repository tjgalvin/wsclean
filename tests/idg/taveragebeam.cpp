#include <boost/test/unit_test.hpp>

#include "../../idg/averagebeam.h"
#include "../../io/cachedimageset.h"

#include <aocommon/image.h>
#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

using aocommon::SerialIStream;
using aocommon::SerialOStream;

namespace {
constexpr size_t kNPol = 4;
constexpr size_t kMatrixWidth = 3;
constexpr size_t kMatrixHeight = 4;
constexpr std::complex<float> kMatrixValue(3.0, 11.0);
constexpr size_t kScalarWidth = 5;
constexpr size_t kScalarHeight = 2;
constexpr float kScalarValue = 4.0;

void CheckBeam(const AverageBeam& average_beam) {
  BOOST_REQUIRE(average_beam.MatrixInverseBeam());
  BOOST_CHECK_EQUAL(average_beam.MatrixInverseBeam()->size(),
                    kMatrixWidth * kMatrixHeight * kNPol * kNPol);
  BOOST_CHECK_EQUAL(average_beam.MatrixWidth(), kMatrixWidth);
  BOOST_CHECK_EQUAL(average_beam.MatrixHeight(), kMatrixHeight);
  BOOST_CHECK_EQUAL(average_beam.MatrixInverseBeam()->at(8), kMatrixValue);
  BOOST_REQUIRE(average_beam.ScalarBeam());
  BOOST_CHECK_EQUAL(average_beam.ScalarBeam()->size(),
                    kScalarWidth * kScalarHeight);
  BOOST_CHECK_EQUAL(average_beam.ScalarWidth(), kScalarWidth);
  BOOST_CHECK_EQUAL(average_beam.ScalarHeight(), kScalarHeight);
  BOOST_CHECK_EQUAL(average_beam.ScalarBeam()->at(7), kScalarValue);
}

AverageBeam FilledBeam() {
  AverageBeam result;
  auto scalar_beam = std::make_shared<std::vector<float>>(
      kScalarWidth * kScalarHeight, kScalarValue);
  auto matrix_beam = std::make_shared<std::vector<std::complex<float>>>(
      kMatrixWidth * kMatrixHeight * kNPol * kNPol, kMatrixValue);
  result.SetScalarBeam(scalar_beam, kScalarWidth, kScalarHeight);
  result.SetMatrixInverseBeam(matrix_beam, kMatrixWidth, kMatrixHeight);
  return result;
}
}  // namespace

BOOST_AUTO_TEST_SUITE(average_beam)

BOOST_AUTO_TEST_CASE(empty_serialization) {
  AverageBeam a;

  SerialOStream ostr;
  a.Serialize(ostr);
  BOOST_CHECK_NE(ostr.size(), 0u);

  SerialIStream istr(std::move(ostr));
  AverageBeam b = FilledBeam();
  b.Unserialize(istr);
  BOOST_CHECK(!b.ScalarBeam());
  BOOST_CHECK(!b.MatrixInverseBeam());
}

BOOST_AUTO_TEST_CASE(filled_serialization) {
  AverageBeam a = FilledBeam();

  SerialOStream ostr;
  a.Serialize(ostr);
  BOOST_CHECK_NE(ostr.size(), 0u);

  SerialIStream istr(std::move(ostr));
  AverageBeam b;
  b.Unserialize(istr);
  CheckBeam(b);
}

BOOST_AUTO_TEST_CASE(empty_store_load) {
  CachedImageSet scalar_cache;
  CachedImageSet matrix_cache;
  FitsWriter writer;
  scalar_cache.Initialize(writer, 2, 1, 0, "test_scalar");
  matrix_cache.Initialize(writer, 2, 1, 0, "test_matrix");

  const size_t frequency_index = 0;
  AverageBeam a;
  a.Store(scalar_cache, matrix_cache, frequency_index);
  BOOST_CHECK(scalar_cache.Empty());
  BOOST_CHECK(matrix_cache.Empty());

  std::unique_ptr<AverageBeam> b =
      AverageBeam::Load(scalar_cache, matrix_cache, frequency_index);
  BOOST_CHECK(b == nullptr);
}

BOOST_AUTO_TEST_CASE(filled_store_load) {
  AverageBeam a = FilledBeam();

  CachedImageSet scalar_cache;
  CachedImageSet matrix_cache;
  FitsWriter writer;
  scalar_cache.Initialize(writer, 2, 1, 0, "test_scalar");
  matrix_cache.Initialize(writer, 2, 1, 0, "test_matrix");

  const size_t frequency_index = 0;
  a.Store(scalar_cache, matrix_cache, frequency_index);
  BOOST_CHECK(!scalar_cache.Empty());
  BOOST_CHECK(!matrix_cache.Empty());

  std::unique_ptr<AverageBeam> b =
      AverageBeam::Load(scalar_cache, matrix_cache, frequency_index);
  CheckBeam(*b);
}

BOOST_AUTO_TEST_SUITE_END()
