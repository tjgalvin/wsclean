#ifndef PAF_BEAM_TERM_H
#define PAF_BEAM_TERM_H

#include <complex>
#include <string>
#include <vector>

#include "atermbase.h"
#include "atermresampler.h"

/**
 * Class for reading fits files as they are used in "Phased-array feed"
 * such as Apertif.
 */
class PAFBeamTerm final : public ATermBase {
 public:
  PAFBeamTerm(const CoordinateSystem& coordinateSystem);

  bool Calculate(std::complex<float>* buffer, double time, double frequency,
                 size_t fieldId, const double* uvwInM) override;

  void SetUpdateInterval(double updateInterval) {
    _updateInterval = updateInterval;
  }

  double AverageUpdateTime() const override { return _updateInterval; }

  void Open(const std::string& filenameTemplate,
            const std::vector<std::string>& antennaMap,
            const std::string& beamName, double beamRA, double beamDec);

  void SetTukeyWindow(double padding) { _resampler.SetTukeyWindow(padding); }

  void SetWindow(WindowFunction::Type window) { _resampler.SetWindow(window); }

  void SetDownSample(bool downsample) { _resampler.SetDownSample(downsample); }

  void SetCorrectForFrequencyOffset(bool correct) {
    _correctForFrequencyOffset = correct;
  }
  void SetReferenceFrequency(double refFrequency) {
    _refFrequency = refFrequency;
  }

 private:
  std::vector<FitsReader> _readers;
  const CoordinateSystem _coordinateSystem;
  ATermResampler _resampler;
  size_t _nAntenna, _nFrequencies;
  double _freq0, _dfreq, _beamRA, _beamDec;
  double _updateInterval, _previousTime;
  bool _correctForFrequencyOffset;
  double _refFrequency;
};

#endif
