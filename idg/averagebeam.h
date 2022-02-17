#ifndef AVERAGE_BEAM_H
#define AVERAGE_BEAM_H

#include <aocommon/io/serialstreamfwd.h>

#include <complex>
#include <memory>
#include <vector>

class CachedImageSet;

class AverageBeam {
 public:
  AverageBeam()
      : scalar_width_(0),
        scalar_height_(0),
        matrix_width_(0),
        matrix_height_(0) {}
  bool Empty() const { return (!scalar_beam_ || !matrix_inverse_beam_); }

  void SetScalarBeam(const std::shared_ptr<std::vector<float>>& scalar_beam,
                     size_t width, size_t height);

  void SetMatrixInverseBeam(
      const std::shared_ptr<std::vector<std::complex<float>>>&
          matrix_inverse_beam,
      size_t width, size_t height);

  /**
   * The image resulting from IDG gridding is multiplied by the scalar beam. It
   * is the result of this multiplication that is returned to WSClean. The
   * scalar beam is chosen such that the Stokes I image is flat noise. There is
   * no guarantee that Stokes Q, U, V will be flat noise, nor that there is no
   * correlation between the noise in I,Q,U,V, but for all practical purposes
   * they can be treated as such.
   *
   * This image has the size of the (trimmed) full gridded image. The size of
   * this image can also be obtained with @ref ScalarWidth() and
   * @ref ScalarHeight().
   */
  std::shared_ptr<std::vector<float>>& ScalarBeam() { return scalar_beam_; }
  std::shared_ptr<const std::vector<float>> ScalarBeam() const {
    return scalar_beam_;
  }
  size_t ScalarWidth() const { return scalar_width_; }
  size_t ScalarHeight() const { return scalar_height_; }

  /**
   * The matrix inverse beam is applied while gridding. It is the inverse of the
   * mean square matrix beam.
   *
   * The returned vector contains sub_grid_size^2 matrices, where
   * each matrix is 4x4 complex numbers.
   */
  std::shared_ptr<std::vector<std::complex<float>>>& MatrixInverseBeam() {
    return matrix_inverse_beam_;
  }
  std::shared_ptr<const std::vector<std::complex<float>>> MatrixInverseBeam()
      const {
    return matrix_inverse_beam_;
  }
  size_t MatrixWidth() const { return matrix_width_; }
  size_t MatrixHeight() const { return matrix_height_; }

  static std::unique_ptr<AverageBeam> Load(const CachedImageSet& scalar_cache,
                                           const CachedImageSet& matrix_cache,
                                           size_t frequency_index);
  void Store(CachedImageSet& scalar_cache, CachedImageSet& matrix_cache,
             size_t frequency_index) const;

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);

 private:
  std::shared_ptr<std::vector<float>> scalar_beam_;
  size_t scalar_width_;
  size_t scalar_height_;
  std::shared_ptr<std::vector<std::complex<float>>> matrix_inverse_beam_;
  size_t matrix_width_;
  size_t matrix_height_;
};

#endif
