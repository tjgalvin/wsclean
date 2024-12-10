#include <boost/test/unit_test.hpp>

#include <aocommon/hmatrix4x4.h>

#include "../../gridding/averagecorrection.h"

#include <cmath>
#include <complex>
#include <iostream>

using aocommon::HMC4x4;
using aocommon::MC2x2;

namespace wsclean {

void CheckMatrix(const HMC4x4& m, const HMC4x4& reference,
                 const std::string& msg = {}) {
  for (size_t i = 0; i != 16; ++i) {
    std::ostringstream str;
    str << "index " << i << ", diff "
        << std::round(std::abs(m[i] - reference[i]) * 100.0) / 100.0 << ", got "
        << m[i] << ", expected " << reference[i];
    if (!msg.empty()) str << "(" << msg << ")";
    BOOST_TEST(std::abs(m[i] - reference[i]) < 1e-2, str.str());
  }
}

BOOST_AUTO_TEST_SUITE(average_correction)

BOOST_AUTO_TEST_CASE(square_root_unit) {
  const HMC4x4 result = PrincipalSquareRoot(HMC4x4::Unit());
  CheckMatrix(result, HMC4x4::Unit());
}

BOOST_AUTO_TEST_CASE(square_root_scalar) {
  const HMC4x4 result = PrincipalSquareRoot(HMC4x4::Unit() * 4.0);
  CheckMatrix(result, HMC4x4::Unit() * 2.0);
}

BOOST_AUTO_TEST_CASE(square_root_diagonal) {
  const HMC4x4 input{1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 4, 0, 0, 0, 0, 16};
  const HMC4x4 result = PrincipalSquareRoot(input.Square());
  CheckMatrix(result, input);
}

BOOST_AUTO_TEST_CASE(square_root_zero_eigenvalue) {
  const std::complex<double> j(0, 1);
  const HMC4x4 input{5, 2.0 + j, 0, 0, 2.0 - j, 5, 0, 0,
                     0, 0,       0, 0, 0,       0, 0, 0};
  const HMC4x4 result = PrincipalSquareRoot(input);
  CheckMatrix(result.Square(), input);
}

BOOST_AUTO_TEST_CASE(square_root_complex) {
  const std::complex<double> j(0, 1);
  const HMC4x4 complex_matrix{2,  j,  j, j, -j, 2,  j, j,
                              -j, -j, 2, 1, -j, -j, 1, 2};
  const HMC4x4 result = PrincipalSquareRoot(complex_matrix);
  CheckMatrix(result.Square(), complex_matrix);
}

BOOST_AUTO_TEST_CASE(square_root_non_positive_semi_definite) {
  const std::complex<double> j(0, 1);
  const HMC4x4 complex_matrix{2,  j,  j, j, -j, 2,  j,  j,
                              -j, -j, 2, j, -j, -j, -j, 2};
  const HMC4x4 result = PrincipalSquareRoot(complex_matrix);
  BOOST_CHECK(!std::isfinite(result[0].real()));
  BOOST_CHECK(!std::isfinite(result[5].real()));
  BOOST_CHECK(!std::isfinite(result[10].real()));
  BOOST_CHECK(!std::isfinite(result[15].real()));
}

void CheckConjugateCorrection(const MC2x2& a, const MC2x2& b) {
  const HMC4x4 partial_a = aocommon::HMC4x4::KroneckerProduct(
      a.HermitianSquare().Transpose(), b.HermitianSquare());
  const HMC4x4 partial_b = aocommon::HMC4x4::KroneckerProduct(
      b.HermitianSquare().Transpose(), a.HermitianSquare());

  const HMC4x4 reference_full = (partial_a + partial_b) * 0.5;
  std::stringstream str_full;
  str_full << "CheckConjugateCorrection(), full: a=" << a << ", b=" << b
           << ", \nleft=" << partial_a.String()
           << "\nright=" << partial_b.String();
  CheckMatrix(AddConjugateCorrectionPart(partial_a), reference_full,
              str_full.str());
}

BOOST_AUTO_TEST_CASE(add_conjugate_correction_part) {
  CheckConjugateCorrection(MC2x2::Zero(), MC2x2::Zero());
  CheckConjugateCorrection(MC2x2::Zero(), MC2x2::Unity());
  CheckConjugateCorrection(MC2x2::Unity(), MC2x2::Unity());
  CheckConjugateCorrection(MC2x2::Unity() * 2.0, MC2x2::Unity() * 3.0);
  const MC2x2 ones{1.0, 1.0, 1.0, 1.0};
  CheckConjugateCorrection(ones, ones);
  CheckConjugateCorrection(ones * -3.0, ones * 2.0);
  CheckConjugateCorrection(ones * std::complex{1.0, 2.0},
                           ones * std::complex{3.0, 4.0});
  const MC2x2 a{1.0, 2.0, 2.0, 1.0};
  CheckConjugateCorrection(a, ones);
  const MC2x2 b{3.0, 4.0, 4.0, 3.0};
  CheckConjugateCorrection(a, b);

  const MC2x2 p{{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0}};
  CheckConjugateCorrection(a, p);
  CheckConjugateCorrection(p, a);

  const MC2x2 q{{16.0, 15.0}, {14.0, 13.0}, {12.0, 11.0}, {10.0, 9.0}};
  CheckConjugateCorrection(p, q);
}

BOOST_AUTO_TEST_CASE(kronecker_square_benchmark_old,
                     *boost::unit_test::disabled()) {
  MC2x2 a{{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0}};
  MC2x2 b{{0.8, 0.7}, {0.6, 0.5}, {0.4, 0.3}, {0.2, 0.1}};
  HMC4x4 r = HMC4x4::Zero();
  for (size_t i = 0; i != 100000000; ++i) {
    // Be aware that we found out these equations miss a transpose. However,
    // the benchmark is still useful to compare how much we gained.
    const aocommon::Matrix4x4 p = aocommon::Matrix4x4::KroneckerProduct(a, b);
    const aocommon::Matrix4x4 q = aocommon::Matrix4x4::KroneckerProduct(b, a);
    r += (p.HermitianSquare() + q.HermitianSquare()) * 0.5;
  }
  BOOST_CHECK(r.Norm() != 0.0);
}

BOOST_AUTO_TEST_CASE(kronecker_square_benchmark_new,
                     *boost::unit_test::disabled()) {
  MC2x2 a{{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0}};
  MC2x2 b{{0.8, 0.7}, {0.6, 0.5}, {0.4, 0.3}, {0.2, 0.1}};
  HMC4x4 r = HMC4x4::Zero();
  for (size_t i = 0; i != 100000000; ++i) {
    r += aocommon::HMC4x4::KroneckerProduct(a.HermitianSquare().Transpose(),
                                            b.HermitianSquare());
  }
  BOOST_CHECK(r.Norm() != 0.0);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace wsclean
