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
    if constexpr (AllowScalarCorrection(Mode)) {
      const std::complex<float> g = GetGainElement<Mode>(gain1, gain2);
      sum_ += std::norm(g) * visibility_weight;
    } else {
      // Add: w * ((g1^H g1)^T (x) g2^H g2 + (g2^H g2)^T (x) g1^H g1)
      const aocommon::MC2x2 g1(gain1);
      const aocommon::MC2x2 g2(gain2);
      const aocommon::HMatrix4x4 a = aocommon::HMatrix4x4::KroneckerProduct(
          g1.HermitianSquare().Transpose(), g2.HermitianSquare());
      const aocommon::HMatrix4x4 b = aocommon::HMatrix4x4::KroneckerProduct(
          g2.HermitianSquare().Transpose(), g1.HermitianSquare());
      matrix_ += (a + b) * (0.5 * visibility_weight);
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
    return sum_ == 0.0 && matrix_ == aocommon::HMatrix4x4::Zero();
  }

  /**
   * True if the correction is a scalar correction.
   */
  constexpr bool IsScalar() const {
    return matrix_ == aocommon::HMatrix4x4::Zero();
  }

  constexpr long double GetScalarValue() const {
    assert(matrix_ == aocommon::HMatrix4x4::Zero());
    return sum_;
  }

  constexpr const aocommon::HMatrix4x4& GetMatrixValue() const {
    assert(sum_ == 0.0);
    return matrix_;
  }

  double GetStokesIValue() const {
    if (matrix_ == aocommon::HMatrix4x4::Zero()) {
      return sum_;
    } else {
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
  std::complex<float> GetGainElement(const aocommon::MC2x2F& gain1,
                                     const aocommon::MC2x2F& gain2) {
    assert(GetNVisibilities(Mode) == 1);
    if constexpr (Mode == GainMode::kXX)
      return gain2.Get(0) * std::conj(gain1.Get(0));
    else if constexpr (Mode == GainMode::kYY)
      return gain2.Get(3) * std::conj(gain1.Get(3));
    else  // Mode == GainMode::kTrace
      return 0.5f * (gain2.Get(0) * std::conj(gain1.Get(0)) +
                     gain2.Get(3) * std::conj(gain1.Get(3)));
  }

  long double sum_ = 0.0L;
  aocommon::HMatrix4x4 matrix_;
};

std::string ToString(const AverageCorrection& average_correction);

}  // namespace wsclean

#endif
