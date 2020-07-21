#ifndef ATERM_BEAM_H
#define ATERM_BEAM_H

#include "atermbase.h"

class ATermBeam : public ATermBase {
 public:
  ATermBeam()
      : _updateInterval(0),
        _lastATermUpdate(0),
        _lastFrequency(0.0),
        _lastFieldId(0) {}

  virtual bool Calculate(std::complex<float>* buffer, double time,
                         double frequency, size_t fieldId,
                         const double*) final override {
    if (time - _lastATermUpdate > _updateInterval || fieldId != _lastFieldId ||
        frequency != _lastFrequency) {
      _lastATermUpdate = time;
      _lastFieldId = fieldId;
      _lastFrequency = frequency;
      return calculateBeam(buffer, time + _updateInterval * 0.5, frequency,
                           fieldId);
    } else {
      return false;
    }
  }

  void SetUpdateInterval(double updateInterval) {
    _updateInterval = updateInterval;
    _lastATermUpdate = -_updateInterval - 1;
    _lastFieldId = 0;
  }

  virtual double AverageUpdateTime() const override { return _updateInterval; }

 protected:
  virtual bool calculateBeam(std::complex<float>* buffer, double time,
                             double frequency, size_t fieldId) = 0;

 private:
  double _updateInterval, _lastATermUpdate, _lastFrequency;
  size_t _lastFieldId;
};

#endif
