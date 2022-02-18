#include "averagebeam.h"

#include "../io/cachedimageset.h"

#include <aocommon/io/serialistream.h>
#include <aocommon/io/serialostream.h>
#include <aocommon/logger.h>

#include <cassert>

namespace {
constexpr size_t kNPolarizations = 4;
}  // namespace

void AverageBeam::SetMatrixInverseBeam(
    const std::shared_ptr<std::vector<std::complex<float>>>&
        matrix_inverse_beam,
    size_t width, size_t height) {
  assert(matrix_inverse_beam->size() ==
         width * height * kNPolarizations * kNPolarizations);
  matrix_inverse_beam_ = matrix_inverse_beam;
  matrix_width_ = width;
  matrix_height_ = height;
}

void AverageBeam::SetScalarBeam(
    const std::shared_ptr<std::vector<float>>& scalar_beam, size_t width,
    size_t height) {
  assert(scalar_beam->size() == width * height);
  scalar_beam_ = scalar_beam;
  scalar_width_ = width;
  scalar_height_ = height;
}

void AverageBeam::Serialize(aocommon::SerialOStream& stream) const {
  stream.Bool(static_cast<bool>(scalar_beam_));
  if (scalar_beam_) {
    stream.Vector(*scalar_beam_).UInt64(scalar_width_).UInt64(scalar_height_);
  }

  stream.Bool(static_cast<bool>(matrix_inverse_beam_));
  if (matrix_inverse_beam_) {
    stream.Vector(*matrix_inverse_beam_)
        .UInt64(matrix_width_)
        .UInt64(matrix_height_);
  }
}

void AverageBeam::Unserialize(aocommon::SerialIStream& stream) {
  const bool has_scalar = stream.Bool();
  if (has_scalar) {
    scalar_beam_.reset(new std::vector<float>());
    stream.Vector(*scalar_beam_).UInt64(scalar_width_).UInt64(scalar_height_);
  } else
    scalar_beam_.reset();

  const bool has_matrix_inverse = stream.Bool();
  if (has_matrix_inverse) {
    matrix_inverse_beam_.reset(new std::vector<std::complex<float>>());
    stream.Vector(*matrix_inverse_beam_)
        .UInt64(matrix_width_)
        .UInt64(matrix_height_);
  } else
    matrix_inverse_beam_.reset();
}

std::unique_ptr<AverageBeam> AverageBeam::Load(
    const CachedImageSet& scalar_cache, const CachedImageSet& matrix_cache,
    size_t frequency_index) {
  std::unique_ptr<AverageBeam> result;
  if (!scalar_cache.Empty()) {
    aocommon::Logger::Debug << "Loading average beam from cache.\n";
    result.reset(new AverageBeam());

    // Scalar beam
    result->scalar_width_ = scalar_cache.Writer().Width();
    result->scalar_height_ = scalar_cache.Writer().Height();
    const size_t n_scalar_pixels =
        result->scalar_width_ * result->scalar_height_;
    result->scalar_beam_ =
        std::make_shared<std::vector<float>>(n_scalar_pixels);
    scalar_cache.Load(result->scalar_beam_->data(),
                      aocommon::PolarizationEnum::StokesI, frequency_index,
                      false);

    // Inverse matrix beam
    // It is stored as one real and one imaginary fits image. The 16 elements
    // are consecutively stored along the X axis, so its width is 16x its
    // height. This makes these images not easy to interpret, but it avoids
    // having to reshuffle the data, and they are only for temporary storage,
    // not for interpretation.
    assert(!matrix_cache.Empty());
    result->matrix_width_ =
        matrix_cache.Writer().Width() / (kNPolarizations * kNPolarizations * 2);
    result->matrix_height_ = matrix_cache.Writer().Height();
    const size_t n_matrix_pixels = result->matrix_width_ *
                                   result->matrix_height_ * kNPolarizations *
                                   kNPolarizations;
    result->matrix_inverse_beam_ =
        std::make_shared<std::vector<std::complex<float>>>(n_matrix_pixels);
    matrix_cache.Load(
        reinterpret_cast<float*>(result->matrix_inverse_beam_->data()),
        aocommon::PolarizationEnum::StokesI, frequency_index, false);
  }
  return result;
}

void AverageBeam::Store(CachedImageSet& scalar_cache,
                        CachedImageSet& matrix_cache,
                        size_t frequency_index) const {
  if (!Empty()) {
    // Scalar beam
    scalar_cache.Writer().SetImageDimensions(scalar_width_, scalar_height_);
    scalar_cache.Store(scalar_beam_->data(),
                       aocommon::PolarizationEnum::StokesI, frequency_index,
                       false);

    // Matrix beam
    matrix_cache.Writer().SetImageDimensions(
        matrix_width_ * kNPolarizations * kNPolarizations * 2, matrix_height_);
    matrix_cache.Store(reinterpret_cast<float*>(matrix_inverse_beam_->data()),
                       aocommon::PolarizationEnum::StokesI, frequency_index,
                       false);
  }
}
