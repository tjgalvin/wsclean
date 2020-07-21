#include "averagebeam.h"

#include "../serialostream.h"
#include "../serialistream.h"

void AverageBeam::Serialize(SerialOStream& stream) const {
  if (_scalarBeam) {
    stream.Bool(true);
    stream.VectorFloat(*_scalarBeam);
  } else {
    stream.Bool(false);
  }

  if (_matrixInverseBeam) {
    stream.Bool(true);
    stream.VectorCFloat(*_matrixInverseBeam);
  } else {
    stream.Bool(false);
  }
}

void AverageBeam::Unserialize(SerialIStream& stream) {
  bool hasScalar = stream.Bool();
  if (hasScalar) {
    _scalarBeam.reset(new std::vector<float>());
    stream.VectorFloat(*_scalarBeam);
  } else
    _scalarBeam.reset();

  bool hasMatrixInverse = stream.Bool();
  if (hasMatrixInverse) {
    _matrixInverseBeam.reset(new std::vector<std::complex<float>>());
    stream.VectorCFloat(*_matrixInverseBeam);
  } else
    _matrixInverseBeam.reset();
}
