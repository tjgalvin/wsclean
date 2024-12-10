#include "averagecorrection.h"

#include <complex>

extern "C" void zheev_(const char* jobz, const char* uplo, const int* n,
                       std::complex<double>* a, const int* lda, double* w,
                       std::complex<double>* work, const int* lwork,
                       double* rwork, int* info);

namespace wsclean {

aocommon::HMC4x4 PrincipalSquareRoot(const aocommon::HMC4x4& matrix) {
  double ev[4];
  constexpr int n = 4;
  constexpr int work_size = 20;
  std::complex<double> work[work_size];
  double rwork[std::max(1, 3 * n - 2)];  // size as required by zheev
  int info = 0;
  std::complex<double> a[16];
  for (size_t col = 0; col != 4; ++col) {
    for (size_t row = 0; row != 4; ++row) {
      // LAPACK uses col-first, HMC4x4 uses row first, so transpose:
      a[col * 4 + row] = matrix[row * 4 + col];
    }
  }
  const char job_mode = 'V';        // Get eigenvalues and eigenvectors
  const char upper_or_lower = 'L';  // Lower triangle of A is stored
  // ZHEEV computes all eigenvalues and, optionally, eigenvectors of a
  // complex Hermitian matrix.
  zheev_(&job_mode, &upper_or_lower,
         &n,  // Order of A
         a,
         &n,  // leading dimension of the array A
         ev, work, &work_size, rwork, &info);
  if (info == 0) {
    bool is_positive_semi_definite = true;
    for (size_t i = 0; i != 4; ++i) {
      if (ev[i] < 0.0) {
        is_positive_semi_definite = false;
        break;
      }
      ev[i] = std::sqrt(ev[i]);
    }

    if (is_positive_semi_definite) {
      const auto c = [](std::complex<double> z) -> std::complex<double> {
        return std::conj(z);
      };
      // Recompose the full matrix: Calculate M = A Λ^1/2 A^H
      // where A is the 4 × 4 matrix whose ith column is the eigenvector qi of
      // the matrix. Λ is the diagonal matrix whose diagonal elements are the
      // corresponding eigenvalues, Λii = λi.
      //
      // Note that LAPACK uses row-first ordering, HMC4x4 uses col first.
      // Therefore, the LAPACK results are transposed Because the last term
      // (A^H) was already transposed, it is transposed twice and thus does not
      // need to be transposed.
      return aocommon::HMC4x4{
          // Row 1
          std::norm(a[0]) * ev[0] + std::norm(a[4]) * ev[1] +
              std::norm(a[8]) * ev[2] + std::norm(a[12]) * ev[3],
          0, 0, 0,
          // Row 2
          a[1] * ev[0] * c(a[0]) + a[5] * ev[1] * c(a[4]) +
              a[9] * ev[2] * c(a[8]) + a[13] * ev[3] * c(a[12]),
          std::norm(a[1]) * ev[0] + std::norm(a[5]) * ev[1] +
              std::norm(a[9]) * ev[2] + std::norm(a[13]) * ev[3],
          0, 0,
          // Row 3
          a[2] * ev[0] * c(a[0]) + a[6] * ev[1] * c(a[4]) +
              a[10] * ev[2] * c(a[8]) + a[14] * ev[3] * c(a[12]),
          a[2] * ev[0] * c(a[1]) + a[6] * ev[1] * c(a[5]) +
              a[10] * ev[2] * c(a[9]) + a[14] * ev[3] * c(a[13]),
          std::norm(a[2]) * ev[0] + std::norm(a[6]) * ev[1] +
              std::norm(a[10]) * ev[2] + std::norm(a[14]) * ev[3],
          0,
          // Row 4
          a[3] * ev[0] * c(a[0]) + a[7] * ev[1] * c(a[4]) +
              a[11] * ev[2] * c(a[8]) + a[15] * ev[3] * c(a[12]),
          a[3] * ev[0] * c(a[1]) + a[7] * ev[1] * c(a[5]) +
              a[11] * ev[2] * c(a[9]) + a[15] * ev[3] * c(a[13]),
          a[3] * ev[0] * c(a[2]) + a[7] * ev[1] * c(a[6]) +
              a[11] * ev[2] * c(a[10]) + a[15] * ev[3] * c(a[14]),
          std::norm(a[3]) * ev[0] + std::norm(a[7]) * ev[1] +
              std::norm(a[11]) * ev[2] + std::norm(a[15]) * ev[3]};
    }
  }
  // Return a matrix with NaNs on the diagonal.
  return aocommon::HMC4x4::Unit() * std::numeric_limits<double>::quiet_NaN();
}

aocommon::HMC4x4 AddConjugateCorrectionPart(const aocommon::HMC4x4& m) {
  using T = const std::complex<double>;
  using RT = const double;
  using std::conj;
  // See HMC4x4::KroneckerProduct() for the kronecker terms
  // r00 = q00 * p00 = m00;
  RT r00 = m[0].real();
  // r10 = q00 * p10 = q00 * conj(p01) = conj(m20)
  T r10 = conj(m[8]);
  // r11 = q00 * p11 = m22
  RT r11 = m[10].real();
  // r20 = q01 * p00 = conj(m10)
  T r20 = conj(m[4]);
  // r21 = q01 * conj(p10) = q01 * p01 = conj(q10) * p01 = m21
  T r21 = m[9];
  // r22 = q11 * p00 = m11
  RT r22 = m[5].real();
  // r30 = q01 * p10 = conj(q10) * conj(p01) = conj(m30)
  T r30 = conj(m[12]);
  // r31 = q01 * p11 = conj(q10) * p11 = conj(m32)
  T r31 = conj(m[14]);
  // r32 = q11 * p10 = q11 * conj(p01) = conj(m31)
  T r32 = conj(m[13]);
  // r33 = q11 * p11 = m33;
  RT r33 = m[15].real();

  return (m + aocommon::HMC4x4::FromData(
                  {r00, r10.real(), r10.imag(), r11, r20.real(), r20.imag(),
                   r21.real(), r21.imag(), r22, r30.real(), r30.imag(),
                   r31.real(), r31.imag(), r32.real(), r32.imag(), r33})) *
         0.5;
}

std::string ToString(const AverageCorrection& average_correction) {
  std::ostringstream str;
  str << "AverageCorrection, ";
  if (average_correction.IsScalar()) {
    str << "scalar: " << average_correction.GetScalarValue();
  } else {
    str << "matrix:\n" << average_correction.GetMatrixValue().String();
  }
  return str.str();
}

}  // namespace wsclean
