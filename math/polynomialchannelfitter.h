#ifndef POLYNOMIAL_CHANNEL_FITTER_H
#define POLYNOMIAL_CHANNEL_FITTER_H

#include <aocommon/uvector.h>

class PolynomialChannelFitter {
 public:
  void Clear() {
    _channels.clear();
    _dataPoints.clear();
  }

  void AddChannel(double startFrequency, double endFrequency) {
    _channels.push_back(std::make_pair(startFrequency, endFrequency));
  }

  void AddDataPoint(size_t channel, double y) {
    _dataPoints.push_back(std::make_pair(channel, y));
  }

  void Fit(aocommon::UVector<double>& terms, size_t nTerms);

  static double Evaluate(double x, const aocommon::UVector<double>& terms) {
    double val = terms[0];
    double f = 1.0;
    for (size_t i = 1; i != terms.size(); ++i) {
      f *= x;
      val += f * terms[i];
    }
    return val;
  }

 private:
  /**
   * Start and end frequencies of the channels
   */
  aocommon::UVector<std::pair<double, double>> _channels;
  aocommon::UVector<std::pair<size_t, double>> _dataPoints;
};

#endif
