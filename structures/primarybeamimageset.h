#ifndef PRIMARY_BEAM_IMAGE_SET
#define PRIMARY_BEAM_IMAGE_SET

#include <array>

#include <boost/filesystem/operations.hpp>

#include <aocommon/hmatrix4x4.h>
#include <aocommon/image.h>
#include <aocommon/matrix2x2.h>
#include <aocommon/polarization.h>
#include <aocommon/fits/fitsreader.h>

#include <radler/component_list.h>

namespace wsclean {

/**
 * A set of images that describe the beam correction. The image set
 * has maximally 16 images that present the 16 elements of a Hermitian Mueller
 * matrix. In the case not all elements are required (e.g. for Stokes I
 * correction), some of the images may remain empty (i.e. unallocated),
 * indicating these elements should be assumed to be zero.
 *
 * All images should have the same size at all time, unless they are empty.
 */
class PrimaryBeamImageSet {
 public:
  PrimaryBeamImageSet() {}

  PrimaryBeamImageSet(size_t width, size_t height) {
    for (size_t i = 0; i != _beamImages.size(); ++i)
      _beamImages[i] = aocommon::Image(width, height);
  }

  void SetToZero() {
    for (aocommon::Image& img : _beamImages) img = 0.0f;
  }

  double GetUnpolarizedCorrectionFactor(size_t x, size_t y) const {
    const double stokes_i = StokesIValue(x + y * Width());
    return stokes_i != 0.0 ? (1.0 / stokes_i) : 0.0;
  }

  const aocommon::Image& operator[](size_t index) const {
    return _beamImages[index];
  }
  aocommon::Image& operator[](size_t index) { return _beamImages[index]; }

  PrimaryBeamImageSet& operator+=(const PrimaryBeamImageSet& rhs) {
    if (_beamImages.size() != rhs._beamImages.size())
      throw std::runtime_error("Primary beam image sets don't match");
    for (size_t i = 0; i != _beamImages.size(); ++i) {
      if (_beamImages[i].Empty())
        _beamImages[i] = rhs._beamImages[i];
      else if (!rhs._beamImages[i].Empty())
        _beamImages[i] += rhs._beamImages[i];
    }
    return *this;
  }

  PrimaryBeamImageSet& operator*=(double factor) {
    for (aocommon::Image& image : _beamImages) image *= factor;
    return *this;
  }

  void ApplyStokesI(float* stokes_i, double beam_limit) const {
    assert(beam_limit >= 0);
    const size_t size = _beamImages[0].Size();
    for (size_t j = 0; j != size; ++j) {
      const double beam_value = StokesIValue(j);
      if (beam_value <= beam_limit) {
        stokes_i[j] = std::numeric_limits<float>::quiet_NaN();
      } else {
        const double inverted_beam = 1.0 / beam_value;
        stokes_i[j] *= inverted_beam;
      }
    }
  }

  void ApplyDiagonal(float* images[2], double beamLimit) const {
    // The beam will be compared to a squared quantity (matrix norm), so square
    // it:
    beamLimit = beamLimit * beamLimit;
    const size_t size = _beamImages[0].Size();
    for (size_t j = 0; j != size; ++j) {
      // xx' = xx * m_00 + yy * m_03
      // yy' = xx * m_30 + yy * m_33
      // Simplify the Mueller matrix to a 2x2 matrix so it can be
      // properly inverted.
      aocommon::MC2x2 beam = Value2x2(j);
      if (Norm(beam) > beamLimit && beam.Invert()) {
        const double xx = images[0][j];
        const double yy = images[1][j];
        images[0][j] = xx * beam.Get(0).real() + yy * beam.Get(1).real();
        images[1][j] = xx * beam.Get(2).real() + yy * beam.Get(3).real();
      } else {
        for (size_t p = 0; p != 2; ++p)
          images[p][j] = std::numeric_limits<float>::quiet_NaN();
      }
    }
  };

  void ApplyFullStokes(float* images[4], double beamLimit) const {
    // The beam will be compared to a squared quantity (matrix norm), so square
    // it:
    beamLimit = beamLimit * beamLimit;
    const size_t size = _beamImages[0].Size();
    for (size_t j = 0; j != size; ++j) {
      aocommon::HMC4x4 beam = Value(j);
      if (beam.Norm() > beamLimit) {
        if (!beam.Invert()) beam = aocommon::HMC4x4::Zero();
        double stokesVal[4] = {images[0][j], images[1][j], images[2][j],
                               images[3][j]};
        aocommon::Vector4 v;
        aocommon::Polarization::StokesToLinear(stokesVal, v.data());
        v = beam * v;
        aocommon::Polarization::LinearToStokes(v.data(), stokesVal);
        for (size_t p = 0; p != 4; ++p) images[p][j] = stokesVal[p];
      } else {
        for (size_t p = 0; p != 4; ++p)
          images[p][j] = std::numeric_limits<float>::quiet_NaN();
      }
    }
  }

  /**
   * @brief Correct component list for primary beam given a (output)
   * channel index.
   */
  void CorrectComponentList(radler::ComponentList& componentList,
                            size_t channel) {
    componentList.MergeDuplicates();

    for (size_t i = 0; i != componentList.NScales(); ++i) {
      const aocommon::UVector<radler::ComponentList::Position>& positions =
          componentList.GetPositions(i);
      for (size_t j = 0; j != positions.size(); ++j) {
        const double beamCorrectionFactor =
            GetUnpolarizedCorrectionFactor(positions[j].x, positions[j].y);
        componentList.MultiplyScaleComponent(i, j, channel,
                                             beamCorrectionFactor);
      }
    }
  }

  static constexpr size_t NImages() { return kNImages; }
  size_t Width() const { return _beamImages[0].Width(); }
  size_t Height() const { return _beamImages[0].Height(); }

 private:
  aocommon::HMC4x4 Value(size_t x, size_t y) const {
    return Value(x + y * Width());
  }

  aocommon::HMC4x4 Value(size_t pixel_index) const {
    aocommon::HMC4x4 beam_values;
    for (size_t i = 0; i != kNImages; ++i) {
      if (!_beamImages[i].Empty()) {
        beam_values.Data(i) = _beamImages[i][pixel_index];
      }
    }
    return beam_values;
  }

  /**
   * Returns a 2x2 matrix for diagonal correction.
   */
  aocommon::MC2x2 Value2x2(size_t pixel_index) const {
    const double xx_to_xx = _beamImages[0][pixel_index];   // m_00
    const double yy_to_yy = _beamImages[15][pixel_index];  // m_33
    const double xx_to_yy_r =
        _beamImages[9].Empty() ? 0.0 : _beamImages[9][pixel_index];
    const double xx_to_yy_i =
        _beamImages[10].Empty() ? 0.0 : _beamImages[10][pixel_index];
    return aocommon::MC2x2{xx_to_xx,
                           {xx_to_yy_r, -xx_to_yy_i},
                           {xx_to_yy_r, xx_to_yy_i},
                           yy_to_yy};
  }

  double StokesIValue(size_t pixel_index) const {
    const double xx_to_xx = _beamImages[0][pixel_index];   // m_00
    const double yy_to_yy = _beamImages[15][pixel_index];  // m_33
    const double xx_to_yy = _beamImages[9].Empty()
                                ? 0.0
                                : _beamImages[9][pixel_index];  // real(m_30)
    // (xx_to_xx + yy_to_yy + m_30 + m_03) / 2
    // = 0.5 * (xx_to_xx + yy_to_yy + m_30 + conj(m_30))
    // = 0.5 * (xx_to_xx + yy_to_yy) + real(m_30)
    return 0.5 * (xx_to_xx + yy_to_yy) + xx_to_yy;
  }

  static constexpr size_t kNImages = 16;
  /**
   * Some images of this array may be left unset if they are not relevant for
   * the beam correction (e.g. for Stokes I correction). The first and last
   * images are always set, unless the entire PrimaryBeamImageSet object is
   * empty.
   */
  std::array<aocommon::Image, kNImages> _beamImages;
};

}  // namespace wsclean

#endif
