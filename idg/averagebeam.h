#ifndef AVERAGE_BEAM_H
#define AVERAGE_BEAM_H

#include <complex>
#include <memory>
#include <vector>

class AverageBeam {
 public:
  AverageBeam() {}
  bool Empty() { return (!_scalarBeam || !_matrixInverseBeam); }
  void SetScalarBeam(const std::shared_ptr<std::vector<float>>& scalarBeam) {
    _scalarBeam = scalarBeam;
  }
  void SetMatrixInverseBeam(
      const std::shared_ptr<std::vector<std::complex<float>>>&
          matrixInverseBeam) {
    _matrixInverseBeam = matrixInverseBeam;
  }

  /**
   * The image resulting from IDG gridding is multiplied by the scalar beam. It
   * is the result of this multiplication that is returned to WSClean. The
   * scalar beam is chosen such that the Stokes I image is flat noise. There is
   * no guarantee that Stokes Q, U, V will be flat noise, nor that there is no
   * correlation between the noise in I,Q,U,V, but for all practical purposes
   * they can be treated as such.
   */
  std::shared_ptr<std::vector<float>>& ScalarBeam() { return _scalarBeam; }

  /**
   * The matrix inverse beam is applied while gridding. It is the inverse of the
   * mean square matrix beam.
   */
  std::shared_ptr<std::vector<std::complex<float>>>& MatrixInverseBeam() {
    return _matrixInverseBeam;
  }

  void Serialize(class SerialOStream& stream) const;
  void Unserialize(class SerialIStream& stream);

 private:
  std::shared_ptr<std::vector<float>> _scalarBeam;
  std::shared_ptr<std::vector<std::complex<float>>> _matrixInverseBeam;
};

#endif
