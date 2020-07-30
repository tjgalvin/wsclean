#ifndef MWA_BEAM_TERM_H
#define MWA_BEAM_TERM_H

#include <complex>
#include <thread>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MPosition.h>

#include <aocommon/lane.h>
#include <aocommon/matrix2x2.h>

#include "atermbeam.h"
#include "atermstub.h"

#ifdef HAVE_EVERYBEAM

#include <EveryBeam/load.h>
#include <EveryBeam/coords/coordutils.h>

class MWABeamTerm : public ATermBeam {
 public:
  MWABeamTerm(casacore::MeasurementSet& ms,
              const CoordinateSystem& coordinateSystem,
              const std::string& search_path);

 private:
  virtual bool calculateBeam(std::complex<float>* buffer, double time,
                             double frequency, size_t fieldId) final override;

  std::unique_ptr<everybeam::telescope::Telescope> _telescope;
  everybeam::coords::CoordinateSystem _coordinate_system;
  // TODO: can be removed, will be stored in telescope::Options anyway
  bool _frequencyInterpolation;
  std::string _coeff_path;
};
#else
using MWABeamTerm = ATermStub;
#endif  // HAVE_EVERYBEAM

#endif
