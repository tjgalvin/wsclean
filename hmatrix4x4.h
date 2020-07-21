#ifndef HMATRIX_4X4_H
#define HMATRIX_4X4_H

#include <complex>
#include <string>
#include <sstream>
#include <stdexcept>

#include <aocommon/matrix2x2.h>
#include <aocommon/matrix4x4.h>

class HMatrix4x4 {
 public:
  HMatrix4x4() { std::fill_n(_data, 16, 0.0); }

  explicit HMatrix4x4(const aocommon::Matrix4x4& src) {
    // row 0
    _data[0] = src[0].real();

    // row 1
    _data[1] = src[4].real();
    _data[2] = src[4].imag();
    _data[3] = src[5].real();

    // row 2
    _data[4] = src[8].real();
    _data[5] = src[8].imag();
    _data[6] = src[9].real();
    _data[7] = src[9].imag();
    _data[8] = src[10].real();

    // row 3
    _data[9] = src[12].real();
    _data[10] = src[12].imag();
    _data[11] = src[13].real();
    _data[12] = src[13].imag();
    _data[13] = src[14].real();
    _data[14] = src[14].imag();
    _data[15] = src[15].real();
  }

  HMatrix4x4(std::initializer_list<std::complex<double>> list) {
    *this = HMatrix4x4(aocommon::Matrix4x4(list));
  }

  static HMatrix4x4 FromData(std::initializer_list<double> list) {
    if (list.size() != 16)
      throw std::runtime_error("FromData parameter should have 16 elements");
    HMatrix4x4 m;
    double* ptr = m._data;
    for (double e : list) {
      *ptr = e;
      ++ptr;
    }
    return m;
  }

  static HMatrix4x4 Zero() { return HMatrix4x4(); }

  static HMatrix4x4 Unit() {
    HMatrix4x4 unit;
    unit._data[0] = 1.0;
    unit._data[3] = 1.0;
    unit._data[8] = 1.0;
    unit._data[15] = 1.0;
    return unit;
  }

  HMatrix4x4& operator+=(const HMatrix4x4& rhs) {
    for (size_t i = 0; i != 16; ++i) _data[i] += rhs._data[i];
    return *this;
  }

  HMatrix4x4 operator*(double rhs) const {
    HMatrix4x4 m;
    for (size_t i = 0; i != 16; ++i) m._data[i] = _data[i] * rhs;
    return m;
  }

  aocommon::Vector4 operator*(const aocommon::Vector4& rhs) const {
    aocommon::Vector4 v(_data[0] * rhs[0], (*this)[4] * rhs[0],
                        (*this)[8] * rhs[0], (*this)[12] * rhs[0]);
    v[0] += (*this)[1] * rhs[1];
    v[1] += _data[3] * rhs[1];
    v[2] += (*this)[1 + 8] * rhs[1];
    v[3] += (*this)[1 + 12] * rhs[1];

    v[0] += (*this)[2] * rhs[2];
    v[1] += (*this)[2 + 4] * rhs[2];
    v[2] += _data[8] * rhs[2];
    v[3] += (*this)[2 + 12] * rhs[2];

    v[0] += (*this)[3] * rhs[3];
    v[1] += (*this)[3 + 4] * rhs[3];
    v[2] += (*this)[3 + 8] * rhs[3];
    v[3] += _data[15] * rhs[3];

    return v;
  }

  HMatrix4x4& operator*=(double rhs) {
    for (size_t i = 0; i != 16; ++i) _data[i] *= rhs;
    return *this;
  }

  bool Invert() {
    aocommon::Matrix4x4 inv = ToMatrix();
    if (!inv.Invert())
      return false;
    else {
      *this = HMatrix4x4(inv);
      return true;
    }
  }

  std::complex<double> operator[](size_t i) const {
    const size_t lookup[16] = {32, 17, 20, 25, 1, 35, 22, 27,
                               4,  6,  40, 29, 9, 11, 13, 47};
    const size_t l = lookup[i];
    return ((l & 32) == 0)
               ? (((l & 16) == 0)
                      ? std::complex<double>(_data[l], _data[l + 1])
                      : std::complex<double>(_data[l & (~16)],
                                             -_data[(l & (~16)) + 1]))
               : (_data[l & (~32)]);
  }

  aocommon::Matrix4x4 ToMatrix() const {
    aocommon::Matrix4x4 m;
    for (size_t i = 0; i != 16; ++i) {
      m[i] = (*this)[i];
    }
    return m;
  }

  double Norm() const {
    return
        // diagonal
        _data[0] * _data[0] + _data[3] * _data[3] + _data[8] * _data[8] +
        _data[15] * _data[15] +
        // lower half x 2
        2.0 * (std::norm(toComplex(1)) + std::norm(toComplex(4)) +
               std::norm(toComplex(6)) + std::norm(toComplex(9)) +
               std::norm(toComplex(11)) + std::norm(toComplex(13)));
  }

  std::string String() const {
    std::ostringstream str;
    for (size_t y = 0; y != 4; ++y) {
      for (size_t x = 0; x != 3; ++x) {
        str << (*this)[x + y * 4] << '\t';
      }
      str << (*this)[3 + y * 4] << '\n';
    }
    return str.str();
  }

  static HMatrix4x4 KroneckerProduct(const aocommon::MC2x2& hma,
                                     const aocommon::MC2x2& hmb) {
    HMatrix4x4 result;

    // top left submatrix
    result._data[0] = (hma[0] * hmb[0]).real();
    result.setComplex(1, hma[0] * hmb[2]);
    result._data[3] = (hma[0] * hmb[3]).real();

    // bottom left submatrix
    result.setComplex(4, hma[2] * hmb[0]);
    result.setComplex(6, hma[2] * hmb[1]);
    result.setComplex(9, hma[2] * hmb[2]);
    result.setComplex(11, hma[2] * hmb[3]);

    // bottom right submatrix
    result._data[8] = (hma[3] * hmb[0]).real();
    result.setComplex(13, hma[3] * hmb[2]);
    result._data[15] = (hma[3] * hmb[3]).real();

    return result;
  }

  const double& Data(size_t index) const { return _data[index]; }
  double& Data(size_t index) { return _data[index]; }

 private:
  std::complex<double> toComplex(size_t singleIndex) const {
    return std::complex<double>(_data[singleIndex], _data[singleIndex + 1]);
  }
  void setComplex(size_t singleIndex, std::complex<double> val) {
    _data[singleIndex] = val.real();
    _data[singleIndex + 1] = val.imag();
  }

  /**
   * Internally, the data is stored as 16 doubles. The diagonal
   * is real, and only the lower (complex) half is stored, in
   * column-first order. The elements can therefore be index in
   * the following way:
   *  0
   *  1  3
   *  4  6 8
   *  9 11 13 15
   */
  double _data[16];
};

typedef HMatrix4x4 HMC4x4;

#endif
