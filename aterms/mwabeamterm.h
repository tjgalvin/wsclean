#ifndef MWA_BEAM_TERM_H
#define MWA_BEAM_TERM_H

#include <thread>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MPosition.h>

#include <complex>

#include "../mwa/tilebeambase.h"
#include "../mwa/tilebeam2016.h"

#include <aocommon/lane.h>
#include <aocommon/matrix2x2.h>

#include "atermbeam.h"
#include "atermstub.h"

class MWABeamTerm : public ATermBeam {
 public:
  MWABeamTerm(casacore::MeasurementSet& ms,
              const CoordinateSystem& coordinateSystem);

  void SetSearchPath(const std::string& searchPath) {
    _searchPath = searchPath;
  }

 private:
  virtual bool calculateBeam(std::complex<float>* buffer, double time,
                             double frequency, size_t fieldId) final override;

  size_t _width, _height, _nStations;
  double _phaseCentreRA, _phaseCentreDec, _dl, _dm, _phaseCentreDL,
      _phaseCentreDM;
  casacore::MPosition _arrayPos;
  double _delays[16];
  bool _frequencyInterpolation;
  std::string _searchPath;
  std::unique_ptr<TileBeamBase<TileBeam2016>> _tileBeam;
};

#endif
