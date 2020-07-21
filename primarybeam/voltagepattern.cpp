#include "voltagepattern.h"

#include "../wsclean/logger.h"
#include "../wsclean/primarybeamimageset.h"

#include <aocommon/imagecoordinates.h>

#include <cmath>

using namespace aocommon;

void VoltagePattern::EvaluatePolynomial(
    const aocommon::UVector<double>& coefficients, bool invert) {
  // This comes from casa's: void PBMath1DIPoly::fillPBArray(), wideband case
  size_t nSamples = 10000;
  size_t nfreq = frequencies.size();
  size_t ncoef = coefficients.size() / nfreq;
  values.resize(nSamples * nfreq);
  inverseIncrementRadius = double(nSamples - 1) / maximumRadiusArcMin;
  double* output = values.data();
  for (size_t n = 0; n != nfreq; n++) {
    const double* freqcoefficients = &coefficients[n * ncoef];
    for (size_t i = 0; i < nSamples; i++) {
      double taper = 0.0;
      double x2 = double(i) / inverseIncrementRadius;
      x2 = x2 * x2;
      double y = 1.0;

      for (size_t j = 0; j < ncoef; j++) {
        taper += y * freqcoefficients[j];
        y *= x2;
      }
      if (taper >= 0.0) {
        if (invert && taper != 0.0) {
          taper = 1.0 / std::sqrt(taper);
        } else {
          taper = std::sqrt(taper);
        }
      } else {
        taper = 0.0;
      }
      *output = taper;
      ++output;
    }
  }
};

aocommon::UVector<double> VoltagePattern::interpolateValues(double freq) const {
  aocommon::UVector<double> result;
  size_t ifit = 0;
  size_t nFreq = frequencies.size();
  for (ifit = 0; ifit != nFreq; ifit++) {
    if (freq <= frequencies[ifit]) break;
  }
  if (ifit == 0) {
    Logger::Info << "Using voltage pattern coefficients for "
                 << frequencies[0] * 1e-6 << " MHz instead of requested "
                 << freq * 1e-6 << '\n';
    result.assign(values.begin(), values.begin() + NSamples());
  } else if (ifit == nFreq) {
    Logger::Info << "Using voltage pattern coefficients for "
                 << frequencies[nFreq - 1] * 1e-6
                 << " MHz instead of requested " << freq * 1e-6 << '\n';
    result.assign(values.begin() + (nFreq - 1) * NSamples(), values.end());
  } else {
    Logger::Info << "Interpolating voltage pattern coefficients from "
                 << frequencies[ifit - 1] * 1e-6 << " and "
                 << frequencies[ifit] * 1e-6 << " MHz to " << freq * 1e-6
                 << " MHz.\n";
    size_t n = NSamples();
    double l = (freq - frequencies[ifit - 1]) /
               (frequencies[ifit] - frequencies[ifit - 1]);
    const double* vpA = FreqIndexValues(ifit - 1);
    const double* vpB = FreqIndexValues(ifit);
    result.resize(n);
    for (size_t i = 0; i != n; ++i) {
      result[i] = vpA[i] * (1.0 - l) + vpB[i] * l;
    }
  }
  return result;
}

const double* VoltagePattern::interpolateValues(
    double frequencyHz, aocommon::UVector<double>& interpolatedVals) const {
  if (frequencies.size() > 1) {
    interpolatedVals = interpolateValues(frequencyHz);
    return interpolatedVals.data();
  } else {
    return FreqIndexValues(0);
  }
}

double VoltagePattern::lmMaxSquared(double frequencyHz) const {
  double factor =
      (180.0 / M_PI) * 60.0 * frequencyHz * 1.0e-9;  // arcminutes * GHz
  double rmax = maximumRadiusArcMin / factor;
  return rmax * rmax;
}

void VoltagePattern::Render(PrimaryBeamImageSet& beamImages, double pixelScaleX,
                            double pixelScaleY, double phaseCentreRA,
                            double phaseCentreDec, double pointingRA,
                            double pointingDec, double phaseCentreDL,
                            double phaseCentreDM, double frequencyHz) const {
  size_t width = beamImages.Width(), height = beamImages.Height();
  double lmMaxSq = lmMaxSquared(frequencyHz);

  aocommon::UVector<double> interpolatedVals;
  const double* vp = interpolateValues(frequencyHz, interpolatedVals);

  double factor =
      (180.0 / M_PI) * 60.0 * frequencyHz * 1.0e-9;  // arcminutes * GHz
  double l0, m0;
  ImageCoordinates::RaDecToLM(pointingRA, pointingDec, phaseCentreRA,
                              phaseCentreDec, l0, m0);
  l0 += phaseCentreDL;
  m0 += phaseCentreDM;
  size_t imgIndex = 0;
  Logger::Debug << "Interpolating 1D voltage pattern to output image...\n";
  for (size_t iy = 0; iy != height; ++iy) {
    for (size_t ix = 0; ix != width; ++ix) {
      double l, m, ra, dec;
      ImageCoordinates::XYToLM(ix, iy, pixelScaleX, pixelScaleY, width, height,
                               l, m);
      l += phaseCentreDL;
      m += m0;
      ImageCoordinates::LMToRaDec(l, m, phaseCentreRA, phaseCentreDec, ra, dec);
      ImageCoordinates::RaDecToLM(ra, dec, pointingRA, pointingDec, l, m);
      l -= l0;
      m -= m0;
      double r2 = l * l + m * m;
      double out;
      if (r2 > lmMaxSq) {
        out = 0.0;
      } else {
        double r = std::sqrt(r2) * factor;
        int indx = int(r * inverseIncrementRadius);
        out = vp[indx];
      }

      beamImages[0][imgIndex] = out;
      beamImages[1][imgIndex] = 0.0;
      beamImages[2][imgIndex] = 0.0;
      beamImages[3][imgIndex] = 0.0;
      beamImages[4][imgIndex] = 0.0;
      beamImages[5][imgIndex] = 0.0;
      beamImages[6][imgIndex] = out;
      beamImages[7][imgIndex] = 0.0;
      ++imgIndex;
    }
  }
}

void VoltagePattern::Render(std::complex<float>* aterm, size_t width,
                            size_t height, double pixelScaleX,
                            double pixelScaleY, double phaseCentreRA,
                            double phaseCentreDec, double pointingRA,
                            double pointingDec, double phaseCentreDL,
                            double phaseCentreDM, double frequencyHz) const {
  double lmMaxSq = lmMaxSquared(frequencyHz);

  aocommon::UVector<double> interpolatedVals;
  const double* vp = interpolateValues(frequencyHz, interpolatedVals);

  double factor =
      (180.0 / M_PI) * 60.0 * frequencyHz * 1.0e-9;  // arcminutes * GHz
  double l0, m0;
  ImageCoordinates::RaDecToLM(pointingRA, pointingDec, phaseCentreRA,
                              phaseCentreDec, l0, m0);
  l0 += phaseCentreDL;
  m0 += phaseCentreDM;
  for (size_t iy = 0; iy != height; ++iy) {
    std::complex<float>* row = aterm + iy * width * 4;
    for (size_t ix = 0; ix != width; ++ix) {
      double l, m, ra, dec;
      ImageCoordinates::XYToLM(ix, iy, pixelScaleX, pixelScaleY, width, height,
                               l, m);
      l += phaseCentreDL;
      m += m0;
      ImageCoordinates::LMToRaDec(l, m, phaseCentreRA, phaseCentreDec, ra, dec);
      ImageCoordinates::RaDecToLM(ra, dec, pointingRA, pointingDec, l, m);
      l -= l0;
      m -= m0;
      double r2 = l * l + m * m;
      double out;
      if (r2 > lmMaxSq) {
        out = 1e-4;
      } else {
        double r = std::sqrt(r2) * factor;
        int indx = int(r * inverseIncrementRadius);
        out = vp[indx] * (1.0 - 1e-4) + 1e-4;
      }

      std::complex<float>* ptr = row + ix * 4;
      ptr[0] = out;
      ptr[1] = 0.0;
      ptr[2] = 0.0;
      ptr[3] = out;
    }
  }
}
