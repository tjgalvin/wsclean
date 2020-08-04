#ifndef WSCLEAN_ATERMS_EVERYBEAMATERMS_H_
#define WSCLEAN_ATERMS_EVERYBEAMATERMS_H_

#include <complex>

#include "atermbeam.h"
#include "atermstub.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#ifdef HAVE_EVERYBEAM
#include <EveryBeam/load.h>
#include <EveryBeam/options.h>
#include <EveryBeam/coords/coordutils.h>

/**
 * @brief Wraps the EveryBeam::Telescope classes
 * for computing the gridded beam response
 *
 */
class EveryBeamATerm final : public ATermBeam {
 public:
  EveryBeamATerm(casacore::MeasurementSet& ms,
                 const CoordinateSystem& coordinateSystem,
                 const everybeam::Options& settings);

 private:
  /**
   * @brief Calculate the gridded response for the \param _telescope
   *
   * @param buffer Buffer
   * @param time Time MJD (s)
   * @param frequency Frequency (Hz)
   * @param fieldId Field id
   */
  bool calculateBeam(std::complex<float>* buffer, double time, double frequency,
                     size_t fieldId) final override;

  std::unique_ptr<everybeam::telescope::Telescope> _telescope;
  everybeam::coords::CoordinateSystem _coordinate_system;

  size_t _cachedFieldId;
  double _cachedFrequency;
};
#else
using EveryBeamATerm = ATermStub;
#endif  // HAVE_EVERYBEAM

#endif