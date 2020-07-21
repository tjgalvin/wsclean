#include "atcabeam.h"

#include <cmath>

#include "../wsclean/logger.h"
#include "../wsclean/primarybeamimageset.h"

enum ATCABeam::Band ATCABeam::GetBand(double freqGHz) {
  if (freqGHz >= 7.0 && freqGHz < 12.5) {
    return ATCA_X;
  } else if (freqGHz > 4.0 && freqGHz < 7.0) {
    return ATCA_C;
  } else if (freqGHz > 1.0 && freqGHz < 3.5) {
    return ATCA_16;
  } else if (freqGHz > 15.0 && freqGHz < 26.0) {
    return ATCA_K;
  } else if (freqGHz > 30.0 && freqGHz < 50.0) {
    return ATCA_Q;
  } else if (freqGHz > 80 && freqGHz < 120.0) {
    return ATCA_W;
  } else {
    return ATCA_L1;  // ATCA_Unknown ?
  }
}

enum ATCABeam::Band ATCABeam::GetBand(const std::string& str) {
  if (str == "ATCA_L1") {
    return ATCA_L1;
  } else if (str == "ATCA_L2") {
    return ATCA_L2;
  } else if (str == "ATCA_L3") {
    return ATCA_L3;
  } else if (str == "ATCA_16") {
    return ATCA_16;
  } else if (str == "ATCA_S") {
    return ATCA_S;
  } else if (str == "ATCA_C") {
    return ATCA_C;
  } else if (str == "ATCA_C_RI") {
    return ATCA_C_RI;
  } else if (str == "ATCA") {
    return ATCA_Unknown;
  } else if (str == "ATCA_X") {
    return ATCA_X;
  } else if (str == "ATCA_K") {
    return ATCA_K;
  } else if (str == "ATCA_Q") {
    return ATCA_Q;
  } else if (str == "ATCA_W") {
    return ATCA_W;
  }
  return ATCA_Unknown;
}

VoltagePattern ATCABeam::CalculateVoltagePattern(enum Band band) {
  VoltagePattern vPattern;
  switch (band) {
    case ATCA_16: {
      Logger::Debug << "Using voltage pattern for ATCA_16 band\n";
      // coef x nfreq
      aocommon::UVector<double> coefficients(
          {1.0, 1.06274e-03, 1.32342e-06, -8.72013e-10, 1.08020e-12,
           1.0, 9.80817e-04, 1.17898e-06, -7.83160e-10, 8.66199e-13,
           1.0, 9.53553e-04, 9.33233e-07, -4.26759e-10, 5.63667e-13,
           1.0, 9.78268e-04, 6.63231e-07, 4.18235e-11,  2.62297e-13,
           1.0, 1.02424e-03, 6.12726e-07, 2.25733e-10,  2.04834e-13,
           1.0, 1.05818e-03, 5.37473e-07, 4.22386e-10,  1.17530e-13,
           1.0, 1.10650e-03, 5.11574e-07, 5.89732e-10,  8.13628e-14});
      vPattern.frequencies.resize(7);
      for (size_t i = 0; i != 7; ++i) {
        vPattern.frequencies[i] = (1332 + i * 256) * 1.e6;
      }
      vPattern.maximumRadiusArcMin = 53.0;
      vPattern.EvaluatePolynomial(coefficients, true);
    } break;
    default:
      throw std::runtime_error(
          "ATCA PB for given spectral band not implemented");
  }
  return vPattern;
}
