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
  PrimaryBeamImageSet() : _beamImages(), _width(0), _height(0) {}

  PrimaryBeamImageSet(size_t width, size_t height)
      : _width(width), _height(height) {
    for (size_t i = 0; i != _beamImages.size(); ++i)
      _beamImages[i] = aocommon::Image(width, height);
  }

  void SetToZero() {
    for (aocommon::Image& img : _beamImages)
      std::fill_n(img.Data(), _width * _height, 0.0);
  }

  double GetUnpolarizedCorrectionFactor(size_t x, size_t y) {
    const size_t j = y * _width + x;
    aocommon::HMC4x4 beam = aocommon::HMC4x4::FromData(
        {_beamImages[0][j], _beamImages[1][j], _beamImages[2][j],
         _beamImages[3][j], _beamImages[4][j], _beamImages[5][j],
         _beamImages[6][j], _beamImages[7][j], _beamImages[8][j],
         _beamImages[9][j], _beamImages[10][j], _beamImages[11][j],
         _beamImages[12][j], _beamImages[13][j], _beamImages[14][j],
         _beamImages[15][j]});
    if (!beam.Invert()) beam = aocommon::HMC4x4::Zero();
    aocommon::Vector4 v{0.5, 0.0, 0.0, 0.5};
    v = beam * v;
    return v[0].real() + v[3].real();
  }

  const aocommon::Image& operator[](size_t index) const {
    return _beamImages[index];
  }
  aocommon::Image& operator[](size_t index) { return _beamImages[index]; }

  PrimaryBeamImageSet& operator+=(const PrimaryBeamImageSet& rhs) {
    if (_beamImages.size() != rhs._beamImages.size())
      throw std::runtime_error("Primary beam image sets don't match");
    for (size_t i = 0; i != _beamImages.size(); ++i) {
      for (size_t j = 0; j != _width * _height; ++j)
        _beamImages[i][j] += rhs._beamImages[i][j];
    }
    return *this;
  }

  PrimaryBeamImageSet& operator*=(double factor) {
    for (size_t i = 0; i != _beamImages.size(); ++i) {
      for (size_t j = 0; j != _width * _height; ++j)
        _beamImages[i][j] *= factor;
    }
    return *this;
  }

  void ApplyScalarStokesI(float* stokes_i, double beam_limit) const {
    // The beam will be compared to a squared quantity (matrix norm), so square
    // it:
    beam_limit = beam_limit * beam_limit;
    const size_t size = _width * _height;
    for (size_t j = 0; j != size; ++j) {
      const aocommon::HMC4x4 beam = aocommon::HMC4x4::FromData(
          {_beamImages[0][j], _beamImages[1][j], _beamImages[2][j],
           _beamImages[3][j], _beamImages[4][j], _beamImages[5][j],
           _beamImages[6][j], _beamImages[7][j], _beamImages[8][j],
           _beamImages[9][j], _beamImages[10][j], _beamImages[11][j],
           _beamImages[12][j], _beamImages[13][j], _beamImages[14][j],
           _beamImages[15][j]});
      const std::array<double, 4> diagonal = beam.DiagonalValues();
      if (beam.Norm() <= beam_limit || diagonal[0] + diagonal[3] == 0.0) {
        stokes_i[j] = std::numeric_limits<float>::quiet_NaN();
      } else {
        const double inverted_beam = 2.0 / (diagonal[0] + diagonal[3]);
        stokes_i[j] *= inverted_beam;
      }
    }
  }

  void ApplyStokesI(float* stokesI, double beamLimit) const {
    // The beam will be compared to a squared quantity (matrix norm), so square
    // it:
    beamLimit = beamLimit * beamLimit;
    const size_t size = _width * _height;
    for (size_t j = 0; j != size; ++j) {
      aocommon::HMC4x4 beam = aocommon::HMC4x4::FromData(
          {_beamImages[0][j], _beamImages[1][j], _beamImages[2][j],
           _beamImages[3][j], _beamImages[4][j], _beamImages[5][j],
           _beamImages[6][j], _beamImages[7][j], _beamImages[8][j],
           _beamImages[9][j], _beamImages[10][j], _beamImages[11][j],
           _beamImages[12][j], _beamImages[13][j], _beamImages[14][j],
           _beamImages[15][j]});
      if (beam.Norm() > beamLimit) {
        if (!beam.Invert()) beam = aocommon::HMC4x4::Zero();
        const float factor = stokesI[j] * 0.5;
        aocommon::Vector4 v{factor, 0.0, 0.0, factor};
        v = beam * v;
        stokesI[j] = v[0].real() + v[3].real();
      } else {
        stokesI[j] = std::numeric_limits<float>::quiet_NaN();
      }
    }
  }

  void ApplyDiagonal(float* images[2], double beamLimit) {
    // The beam will be compared to a squared quantity (matrix norm), so square
    // it:
    beamLimit = beamLimit * beamLimit;
    const size_t size = _width * _height;
    for (size_t j = 0; j != size; ++j) {
      aocommon::HMC4x4 beam = aocommon::HMC4x4::FromData(
          {_beamImages[0][j], _beamImages[1][j], _beamImages[2][j],
           _beamImages[3][j], _beamImages[4][j], _beamImages[5][j],
           _beamImages[6][j], _beamImages[7][j], _beamImages[8][j],
           _beamImages[9][j], _beamImages[10][j], _beamImages[11][j],
           _beamImages[12][j], _beamImages[13][j], _beamImages[14][j],
           _beamImages[15][j]});
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
      aocommon::HMC4x4 beam = aocommon::HMC4x4::FromData(
          {_beamImages[0][j], _beamImages[1][j], _beamImages[2][j],
           _beamImages[3][j], _beamImages[4][j], _beamImages[5][j],
           _beamImages[6][j], _beamImages[7][j], _beamImages[8][j],
           _beamImages[9][j], _beamImages[10][j], _beamImages[11][j],
           _beamImages[12][j], _beamImages[13][j], _beamImages[14][j],
           _beamImages[15][j]});
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

  size_t NImages() const { return _beamImages.size(); }
  size_t Width() const { return _width; }
  size_t Height() const { return _height; }

 private:
  std::array<aocommon::Image, 16> _beamImages;
  size_t _width, _height;
};

#endif
