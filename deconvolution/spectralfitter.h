#ifndef SPECTRAL_FITTER_H
#define SPECTRAL_FITTER_H

#include <aocommon/uvector.h>

enum SpectralFittingMode {
  NoSpectralFitting,
  PolynomialSpectralFitting,
  LogPolynomialSpectralFitting
};

class SpectralFitter {
 public:
  SpectralFitter(SpectralFittingMode mode, size_t nTerms)
      : _mode(mode), _nTerms(nTerms) {}

  SpectralFittingMode Mode() const { return _mode; }

  void SetMode(SpectralFittingMode mode, size_t nTerms) {
    _mode = mode;
    _nTerms = nTerms;
  }

  void FitAndEvaluate(double* values) const;

  void Fit(aocommon::UVector<double>& terms, const double* values) const;

  void Evaluate(double* values, const aocommon::UVector<double>& terms) const;

  double Evaluate(const aocommon::UVector<double>& terms,
                  double frequency) const;

  void SetFrequencies(const double* frequencies, const double* weights,
                      size_t n) {
    _frequencies.assign(frequencies, frequencies + n);
    _weights.assign(weights, weights + n);
    double weightSum = 0.0;
    _referenceFrequency = 0.0;
    for (size_t i = 0; i != n; ++i) {
      _referenceFrequency += _frequencies[i] * _weights[i];
      weightSum += _weights[i];
    }
    if (weightSum != 0.0)
      _referenceFrequency /= weightSum;
    else
      _referenceFrequency = 150e6;
  }

  double Frequency(size_t index) const { return _frequencies[index]; }

  double Weight(size_t index) const { return _weights[index]; }

  size_t NTerms() const { return _nTerms; }

  size_t NFrequencies() const { return _frequencies.size(); }

  double ReferenceFrequency() const { return _referenceFrequency; }

 private:
  enum SpectralFittingMode _mode;
  size_t _nTerms;
  aocommon::UVector<double> _frequencies, _weights;
  double _referenceFrequency;
};

#endif
