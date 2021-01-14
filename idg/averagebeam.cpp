#include "averagebeam.h"

#include <aocommon/io/serialistream.h>
#include <aocommon/io/serialostream.h>

void AverageBeam::Serialize(aocommon::SerialOStream& stream) const {
  if (_scalarBeam) {
    stream.Bool(true);
    stream.Vector(*_scalarBeam);
  } else {
    stream.Bool(false);
  }

  if (_matrixInverseBeam) {
    stream.Bool(true);
    stream.Vector(*_matrixInverseBeam);
  } else {
    stream.Bool(false);
  }
}

void AverageBeam::Unserialize(aocommon::SerialIStream& stream) {
  bool hasScalar = stream.Bool();
  if (hasScalar) {
    _scalarBeam.reset(new std::vector<float>());
    stream.Vector(*_scalarBeam);
  } else
    _scalarBeam.reset();

  bool hasMatrixInverse = stream.Bool();
  if (hasMatrixInverse) {
    _matrixInverseBeam.reset(new std::vector<std::complex<float>>());
    stream.Vector(*_matrixInverseBeam);
  } else
    _matrixInverseBeam.reset();
}
