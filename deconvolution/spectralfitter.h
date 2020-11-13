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
  typedef float num_t;

  SpectralFitter(SpectralFittingMode mode, size_t nTerms)
      : _mode(mode), _nTerms(nTerms) {}

  SpectralFittingMode Mode() const { return _mode; }

  void SetMode(SpectralFittingMode mode, size_t nTerms) {
    _mode = mode;
    _nTerms = nTerms;
  }

  void FitAndEvaluate(num_t* values) const;

  void Fit(aocommon::UVector<num_t>& terms, const num_t* values) const;

  void Evaluate(num_t* values, const aocommon::UVector<num_t>& terms) const;

  num_t Evaluate(const aocommon::UVector<num_t>& terms, double frequency) const;

  void SetFrequencies(const double* frequencies, const num_t* weights,
                      size_t n) {
    _frequencies.assign(frequencies, frequencies + n);
    _weights.assign(weights, weights + n);
    num_t weightSum = 0.0;
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

  num_t Weight(size_t index) const { return _weights[index]; }

  size_t NTerms() const { return _nTerms; }

  size_t NFrequencies() const { return _frequencies.size(); }

  double ReferenceFrequency() const { return _referenceFrequency; }

 private:
  enum SpectralFittingMode _mode;
  size_t _nTerms;
  aocommon::UVector<double> _frequencies;
  aocommon::UVector<num_t> _weights;
  double _referenceFrequency;
};

#endif
