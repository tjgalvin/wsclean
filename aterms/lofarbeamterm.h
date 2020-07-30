#ifndef LOFAR_BEAM_TERM_H
#define LOFAR_BEAM_TERM_H

#include <thread>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <complex>

#include "atermstub.h"
#include "atermbeam.h"

#ifdef HAVE_EVERYBEAM

#include <EveryBeam/load.h>
#include <EveryBeam/coords/coordutils.h>

class LofarBeamTerm : public ATermBeam {
 public:
  LofarBeamTerm(casacore::MeasurementSet& ms,
                const CoordinateSystem& coordinateSystem,
                const std::string& dataColumnName);

  void SetUseDifferentialBeam(bool useDiffBeam) {
    _useDifferentialBeam = useDiffBeam;
  }

  void SetUseChannelFrequency(bool useChannelFrequency) {
    _useChannelFrequency = useChannelFrequency;
  }

 private:
  bool calculateBeam(std::complex<float>* buffer, double time, double frequency,
                     size_t fieldId) final override;

  std::unique_ptr<everybeam::telescope::Telescope> _telescope;
  everybeam::coords::CoordinateSystem _coordinate_system;
  // TODO: can be removed, will be stored in telescope::Options anyway
  bool _useDifferentialBeam, _useChannelFrequency;
};

#else
using LofarBeamTerm = ATermStub;
#endif  // HAVE_EVERYBEAM

#endif
