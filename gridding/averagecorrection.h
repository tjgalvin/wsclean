#ifndef GRIDDING_AVERAGE_CORRECTION_H_
#define GRIDDING_AVERAGE_CORRECTION_H_

#include <aocommon/hmatrix4x4.h>
#include <aocommon/io/serialistream.h>
#include <aocommon/io/serialostream.h>

#include "gainmode.h"

namespace wsclean {

/**
 * Returns the principal square root of an Hermitian matrix.
 * Given matrix A, the returned matrix B is such that A = BB.
 * Because B is guaranteed to be Hermitian, A = BB^H also holds.
 *
 * The specified matrix should be positive (semi)definite. If
 * not, the eigen values may be negative, leading to values
 * on the diagonal of the square root matrix that are complex,
 * and thus the matrix is no longer Hermitian.
 *
 * The principal square root of a Hermitian positive semi-definite
 * matrix is unique, and is itself also a Hermitian positive
 * semi-definite matrix.
 */
aocommon::HMC4x4 PrincipalSquareRoot(const aocommon::HMC4x4& matrix);

/**
 * Given a matrix that was formed from KroneckerProduct(a^T, b), this
 * functions constructs 0.5 * [KroneckerProduct(a^T, b) + KroneckerProduct(b^T,
 * a)].
 */
aocommon::HMC4x4 AddConjugateCorrectionPart(const aocommon::HMC4x4& m);

namespace internal {
/**
 * Given a diagonal matrix that is to be used in a KroneckerProduct.
 * This function computes the Hermitian square of the matrix with
 * all results that won't be used in the KroneckerProduct optimised out.
 * In other words:
 *   Skip the transpose and all multiplications involving the zero diagonals.
 *   Skip computing the imaginary part of the final complex multiplication.
 *   Only the real part is needed for the optimised KroneckerProduct.
 */
inline std::array<double, 2> PartialHermitianSquareForDiagonalKronecker(
    const aocommon::MC2x2FDiag& matrix) {
  const aocommon::MC2x2FDiag conjugate = matrix.Conjugate();
  return {(std::real(conjugate.Get(0)) * double{std::real(matrix.Get(0))}) -
              (std::imag(conjugate.Get(0)) * double{std::imag(matrix.Get(0))}),
          (std::real(conjugate.Get(1)) * double{std::real(matrix.Get(1))}) -
              (std::imag(conjugate.Get(1)) * double{std::imag(matrix.Get(1))})};
}
}  // namespace internal

/**
 * This class is used to collect the Jones corrections that are applied
 * to visibilities while gridding. The Kronecker product of the Jones
 * matrices is taken to form a Mueller matrix. The Mueller matrices are
 * squared and averaged.
 */
class AverageCorrection {
 public:
  template <GainMode Mode>
  void Add(const aocommon::MC2x2F& gain1, const aocommon::MC2x2F& gain2,
           float visibility_weight) {
    static_assert(!AllowScalarCorrection(Mode),
                  "Use MC2x2Diag for scalar corrections");
    // Add: w * [ (g1^H g1)^T (x) (g2^H g2) ].
    // The conjugate part is added later (see AddConjugateCorrectionPart()).
    const aocommon::MC2x2 g1(gain1);
    const aocommon::MC2x2 g2(gain2);
    matrix_ += aocommon::HMC4x4::KroneckerProduct(
                   g1.HermitianSquare().Transpose(), g2.HermitianSquare()) *
               visibility_weight;
  }

  /**
   * Optimised specialisation of @ref Add() for diagonal matrices
   */
  template <GainMode Mode>
  void Add(const aocommon::MC2x2FDiag& gain1, const aocommon::MC2x2FDiag& gain2,
           double visibility_weight) {
    using internal::PartialHermitianSquareForDiagonalKronecker;
    if constexpr (AllowScalarCorrection(Mode)) {
      const std::complex<float> g = GetGainElement<Mode>(gain1, gain2);
      sum_ += std::norm(g) * visibility_weight;
    } else {
      // Compute the hermitian square of gain1 and gain2.
      // In an optimised way that skips computation of parts we don't need.
      // Skip transpose of gain1 Hermitian square, this is a no-op on a diagonal
      // matrix.
      std::array<double, 2> g1_herm_square_real =
          PartialHermitianSquareForDiagonalKronecker(gain1);
      std::array<double, 2> g2_herm_square_real =
          PartialHermitianSquareForDiagonalKronecker(gain2);

      // Compute the KroneckerProduct of gain1 and gain2:
      //   * Only compute the real values for full matrix elements:
      //     {0,0}, {0,3}, {3,0} and {3,3}
      //   * These have indexes 0, 3, 8 and 15 in the HMatrix class.
      //   * All other values are guaranteed to be zero.
      matrix_ += aocommon::HMatrix4x4::FromData({
          g1_herm_square_real[0] * g2_herm_square_real[0] * visibility_weight,
          0.0f,
          0.0f,
          g1_herm_square_real[0] * g2_herm_square_real[1] * visibility_weight,
          0.0f,
          0.0f,
          0.0f,
          0.0f,
          g1_herm_square_real[1] * g2_herm_square_real[0] * visibility_weight,
          0.0f,
          0.0f,
          0.0f,
          0.0f,
          0.0f,
          0.0f,
          g1_herm_square_real[1] * g2_herm_square_real[1] * visibility_weight,
      });
    }
  }

  AverageCorrection operator/(double denominator) const {
    AverageCorrection result;
    result.sum_ = sum_ / denominator;
    result.matrix_ = matrix_ * (1.0 / denominator);
    return result;
  }

  void Serialize(aocommon::SerialOStream& stream) const {
    stream.LDouble(sum_).Object(matrix_);
  }

  void Unserialize(aocommon::SerialIStream& stream) {
    stream.LDouble(sum_).Object(matrix_);
  }

  /**
   * True if the correction is completely zero.
   */
  constexpr bool IsZero() const {
    return sum_ == 0.0 && matrix_ == aocommon::HMC4x4::Zero();
  }

  /**
   * True if the correction is a scalar correction.
   */
  constexpr bool IsScalar() const {
    return matrix_ == aocommon::HMC4x4::Zero();
  }

  constexpr long double GetScalarValue() const {
    assert(matrix_ == aocommon::HMC4x4::Zero());
    return sum_;
  }

  const aocommon::HMC4x4 GetMatrixValue() const {
    assert(sum_ == 0.0);
    return AddConjugateCorrectionPart(matrix_);
  }

  double GetStokesIValue() const {
    if (matrix_ == aocommon::HMC4x4::Zero()) {
      return sum_;
    } else {
      // The matrix hasn't been finalized with
      // AddConjugateCorrectionPart(matrix_) yet. However, doing so doesn't
      // change the result of the sum and can therefore be skipped. Proof: if r
      // is the finalized matrix and m is matrix_, then:
      // - r_00 = 0.5 * (m_00 + m_00))
      // - r_30 = 0.5 * (m_30 + conj(m_30))
      // - r_33 = 0.5 * (m_33 + m_33)
      // (Equations are from AddConjugateCorrectionPart()).
      // Hence, r_00 and r_33 are not changed, and the real part of r_30 is also
      // not changed (and is the only used part).
      return 0.5 * (matrix_.Data(0) + 2.0 * matrix_.Data(9) + matrix_.Data(15));
    }
  }

 private:
  /**
   * @brief Compute the gain from the given solution matrices.
   *
   * @tparam Mode Which entry or entries from the gain matrices should be
   * taken into account? Must be a mode for which NVisibilities(Mode) == 1.
   * See @ref GainMode for further documentation.
   */
  template <GainMode Mode>
  std::complex<float> GetGainElement(const aocommon::MC2x2FDiag& gain1,
                                     const aocommon::MC2x2FDiag& gain2) {
    assert(GetNVisibilities(Mode) == 1);
    if constexpr (Mode == GainMode::kXX)
      return gain2.Get(0) * std::conj(gain1.Get(0));
    else if constexpr (Mode == GainMode::kYY)
      return gain2.Get(1) * std::conj(gain1.Get(1));
    else  // Mode == GainMode::kTrace
      return 0.5f * (gain2.Get(0) * std::conj(gain1.Get(0)) +
                     gain2.Get(1) * std::conj(gain1.Get(1)));
  }

  long double sum_ = 0.0L;
  aocommon::HMC4x4 matrix_ = aocommon::HMC4x4::Zero();
};

std::string ToString(const AverageCorrection& average_correction);

}  // namespace wsclean

#endif
