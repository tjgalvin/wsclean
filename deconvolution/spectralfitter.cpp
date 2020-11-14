#include "spectralfitter.h"

#include "../math/polynomialfitter.h"
#include "../math/nlplfitter.h"

#include "../io/logger.h"

#include <limits>

void SpectralFitter::FitAndEvaluate(num_t* values) const {
  aocommon::UVector<num_t> terms;
  Fit(terms, values);
  Evaluate(values, terms);
}

void SpectralFitter::Fit(aocommon::UVector<num_t>& terms,
                         const num_t* values) const {
  switch (_mode) {
    default:
    case NoSpectralFitting:
      break;

    case PolynomialSpectralFitting: {
      PolynomialFitter fitter;
      double refFreq = ReferenceFrequency();
      for (size_t i = 0; i != _frequencies.size(); ++i) {
        fitter.AddDataPoint(_frequencies[i] / refFreq - 1.0, values[i],
                            _weights[i]);
      }

      fitter.Fit(terms, _nTerms);
    } break;

    case LogPolynomialSpectralFitting: {
      NonLinearPowerLawFitter fitter;
      double refFreq = ReferenceFrequency();
      for (size_t i = 0; i != _frequencies.size(); ++i)
        fitter.AddDataPoint(_frequencies[i] / refFreq, values[i]);

      fitter.Fit(terms, _nTerms);
    } break;
  }
}

void SpectralFitter::Evaluate(num_t* values,
                              const aocommon::UVector<num_t>& terms) const {
  switch (_mode) {
    default:
    case NoSpectralFitting:
      break;

    case PolynomialSpectralFitting: {
      double refFreq = ReferenceFrequency();
      for (size_t i = 0; i != _frequencies.size(); ++i) {
        num_t newValue =
            PolynomialFitter::Evaluate(_frequencies[i] / refFreq - 1.0, terms);
        // std::cout << values[i] << "->" << newValue << ' ';
        values[i] = newValue;
      }
      // std::cout << '\n';

    } break;

    case LogPolynomialSpectralFitting: {
      double refFreq = ReferenceFrequency();
      for (size_t i = 0; i != _frequencies.size(); ++i) {
        num_t newValue =
            NonLinearPowerLawFitter::Evaluate(_frequencies[i], terms, refFreq);
        // std::cout << values[i] << "->" << newValue << ' ';
        values[i] = newValue;
      }
      // std::cout << '\n';

    } break;
  }
}

SpectralFitter::num_t SpectralFitter::Evaluate(
    const aocommon::UVector<num_t>& terms, double frequency) const {
  switch (_mode) {
    default:
    case NoSpectralFitting:
      throw std::runtime_error(
          "Something is inconsistent: can't evaluate terms at frequency "
          "without fitting");

    case PolynomialSpectralFitting:
      return PolynomialFitter::Evaluate(frequency / ReferenceFrequency() - 1.0,
                                        terms);

    case LogPolynomialSpectralFitting:
      return NonLinearPowerLawFitter::Evaluate(frequency, terms,
                                               ReferenceFrequency());
  }
}
