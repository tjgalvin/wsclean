#include "spectralfitter.h"

#include "../math/polynomialfitter.h"
#include "../math/nlplfitter.h"

#include "../io/logger.h"

#include <limits>

void SpectralFitter::Fit(aocommon::UVector<num_t>& terms, const num_t* values,
                         size_t x, size_t y) const {
  if (IsForced()) {
    terms.resize(_nTerms);
    // We need to find a such that
    // y[i] = a f(x[i], terms), with f the shape.
    // The least-squares fit is:
    // a = sum (y[i] w[i] f[i]) / sum (w[i] f[i]^2)
    terms[0] = 1.0;
    for (size_t term = 1; term != _nTerms; ++term) {
      const Image& termImage = _forcedTerms[term - 1];
      terms[term] = termImage[x + y * termImage.Width()];
    }
    float aNumerator = 0.0;
    float aDivisor = 0.0;
    for (size_t i = 0; i != _frequencies.size(); ++i) {
      const float w = _weights[i];
      const float f = NonLinearPowerLawFitter::Evaluate(_frequencies[i], terms,
                                                        ReferenceFrequency());
      aNumerator += w * f * values[i];
      aDivisor += w * f * f;
    }
    terms[0] = aNumerator / aDivisor;
  } else {
    switch (_mode) {
      default:
      case SpectralFittingMode::NoFitting:
        break;

      case SpectralFittingMode::Polynomial: {
        PolynomialFitter fitter;
        double refFreq = ReferenceFrequency();
        for (size_t i = 0; i != _frequencies.size(); ++i) {
          fitter.AddDataPoint(_frequencies[i] / refFreq - 1.0, values[i],
                              _weights[i]);
        }

        fitter.Fit(terms, _nTerms);
      } break;

      case SpectralFittingMode::LogPolynomial: {
        NonLinearPowerLawFitter fitter;
        double refFreq = ReferenceFrequency();
        for (size_t i = 0; i != _frequencies.size(); ++i)
          fitter.AddDataPoint(_frequencies[i] / refFreq, values[i]);

        fitter.Fit(terms, _nTerms);
      } break;
    }
  }
}

void SpectralFitter::Evaluate(num_t* values,
                              const aocommon::UVector<num_t>& terms) const {
  switch (_mode) {
    default:
    case SpectralFittingMode::NoFitting:
      break;

    case SpectralFittingMode::Polynomial: {
      const double refFreq = ReferenceFrequency();
      for (size_t i = 0; i != _frequencies.size(); ++i) {
        num_t newValue =
            PolynomialFitter::Evaluate(_frequencies[i] / refFreq - 1.0, terms);
        values[i] = newValue;
      }
    } break;

    case SpectralFittingMode::LogPolynomial: {
      const double refFreq = ReferenceFrequency();
      for (size_t i = 0; i != _frequencies.size(); ++i) {
        num_t newValue =
            NonLinearPowerLawFitter::Evaluate(_frequencies[i], terms, refFreq);
        values[i] = newValue;
      }
    } break;
  }
}

SpectralFitter::num_t SpectralFitter::Evaluate(
    const aocommon::UVector<num_t>& terms, double frequency) const {
  switch (_mode) {
    default:
    case SpectralFittingMode::NoFitting:
      throw std::runtime_error(
          "Something is inconsistent: can't evaluate terms at frequency "
          "without fitting");

    case SpectralFittingMode::Polynomial:
      return PolynomialFitter::Evaluate(frequency / ReferenceFrequency() - 1.0,
                                        terms);

    case SpectralFittingMode::LogPolynomial:
      return NonLinearPowerLawFitter::Evaluate(frequency, terms,
                                               ReferenceFrequency());
  }
}
