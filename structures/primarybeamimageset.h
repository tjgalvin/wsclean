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

class PrimaryBeamImageSet {
 public:
  PrimaryBeamImageSet() : _width(0), _height(0) {}

  PrimaryBeamImageSet(size_t width, size_t height)
      : _width(width), _height(height) {
    for (size_t i = 0; i != _beamImages.size(); ++i)
      _beamImages[i] = aocommon::Image(width, height);
  }

  void SetToZero() {
    for (aocommon::Image& img : _beamImages) img = 0.0f;
  }

  double GetUnpolarizedCorrectionFactor(size_t x, size_t y) const {
    const aocommon::HMC4x4 beam = Value(x, y);
    const double sum = beam.Data(0) + beam.Data(15);
    return sum != 0.0 ? (2.0 / sum) : 0.0;
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
    // The beam will be compared to a squared quantity (matrix norm), so square
    // it:
    beam_limit = beam_limit * beam_limit;
    const size_t size = _width * _height;
    for (size_t j = 0; j != size; ++j) {
      const aocommon::HMC4x4 beam = Value(j);
      const std::array<double, 4> diagonal = beam.DiagonalValues();
      if (beam.Norm() <= beam_limit || diagonal[0] + diagonal[3] == 0.0) {
        stokes_i[j] = std::numeric_limits<float>::quiet_NaN();
      } else {
        const double inverted_beam = 2.0 / (diagonal[0] + diagonal[3]);
        stokes_i[j] *= inverted_beam;
      }
    }
  }

  void ApplyDiagonal(float* images[2], double beamLimit) const {
    // The beam will be compared to a squared quantity (matrix norm), so square
    // it:
    beamLimit = beamLimit * beamLimit;
    const size_t size = _width * _height;
    for (size_t j = 0; j != size; ++j) {
      aocommon::HMC4x4 beam = Value(j);
      if (beam.Norm() > beamLimit) {
        if (!beam.Invert()) beam = aocommon::HMC4x4::Zero();
        aocommon::Vector4 v{images[0][j], 0., 0., images[1][j]};
        // Implementation could be made more efficient - at the expense of
        // explicitness - since only a few entries in the beam matrix are
        // relevant
        v = beam * v;
        images[0][j] = v[0].real();
        images[1][j] = v[3].real();
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
    const size_t size = _width * _height;
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
  size_t Width() const { return _width; }
  size_t Height() const { return _height; }

 private:
  aocommon::HMC4x4 Value(size_t x, size_t y) const {
    return Value(x + y * _width);
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

  static constexpr size_t kNImages = 16;
  std::array<aocommon::Image, kNImages> _beamImages;
  size_t _width, _height;
};

#endif
