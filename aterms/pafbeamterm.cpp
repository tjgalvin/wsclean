#include "pafbeamterm.h"

#include <boost/algorithm/string.hpp>

#include "../wsclean/logger.h"

#include "../fitsreader.h"

PAFBeamTerm::PAFBeamTerm(const ATermBase::CoordinateSystem& coordinateSystem)
    : _coordinateSystem(coordinateSystem),
      _resampler(coordinateSystem),
      _freq0(0.0),
      _dfreq(0.0),
      _beamRA(0.0),
      _beamDec(0.0),
      _updateInterval(3600),
      _previousTime(std::numeric_limits<double>::lowest()),
      _correctForFrequencyOffset(true),
      _refFrequency(0.0) {}

bool PAFBeamTerm::Calculate(std::complex<float>* buffer, double time,
                            double frequency, size_t /*fieldId*/,
                            const double* /*uvwInM*/) {
  bool outdated = std::fabs(time - _previousTime) > _updateInterval;
  if (!outdated) return false;
  _previousTime = time;

  const size_t width = _coordinateSystem.width,
               height = _coordinateSystem.height,
               freqIndex = std::min<size_t>(
                   _nFrequencies - 1,
                   std::max(0.0, round((frequency - _freq0) / _dfreq)));

  double refFrequency;
  if (_refFrequency != 0.0)
    refFrequency = _refFrequency;
  else
    refFrequency = double(freqIndex) * _dfreq + _freq0;
  const double frequencyFactor =
      _correctForFrequencyOffset ? frequency / refFrequency : 1.0;

  aocommon::UVector<double> scratch(_resampler.ScratchASize());
  aocommon::UVector<double> output(_resampler.ScratchBSize(_readers[0]));
  for (size_t antenna = 0; antenna != _nAntenna; ++antenna) {
    _resampler.OverrideFitsPhaseCentre(_beamRA, _beamDec);
    _resampler.ReadAndResample(_readers[antenna], freqIndex, scratch, output,
                               frequencyFactor);
    std::complex<float>* imgPtr = buffer + antenna * width * height * 4;
    for (size_t i = 0; i != width * height; ++i) {
      imgPtr[i * 4] = output[i];
      imgPtr[i * 4 + 1] = 0.0;
      imgPtr[i * 4 + 2] = 0.0;
      imgPtr[i * 4 + 3] = output[i];
    }
  }

  return true;
}

void PAFBeamTerm::Open(const std::string& filenameTemplate,
                       const std::vector<std::string>& antennaMap,
                       const std::string& beamName, double beamRA,
                       double beamDec) {
  _nAntenna = antennaMap.size();
  _readers.clear();
  _beamRA = beamRA;
  _beamDec = beamDec;
  for (size_t antenna = 0; antenna != _nAntenna; ++antenna) {
    std::string filename = boost::algorithm::replace_all_copy(
        filenameTemplate, "$ANT", antennaMap[antenna]);
    boost::algorithm::replace_all(filename, "$BEAM", beamName);
    _readers.emplace_back(filename, false, true);
    if (antenna == 0) {
      _freq0 = _readers.back().FrequencyDimensionStart();
      _dfreq = _readers.back().FrequencyDimensionIncr();
      _nFrequencies = _readers.back().NFrequencies();
    } else {
      if (_freq0 != _readers.back().FrequencyDimensionStart() ||
          _dfreq != _readers.back().FrequencyDimensionIncr() ||
          _nFrequencies != _readers.back().NFrequencies())
        throw std::runtime_error(
            "All antenna fits files should have the same frequency axis");
    }
  }
}
