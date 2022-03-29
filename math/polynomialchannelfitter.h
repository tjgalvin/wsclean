#ifndef POLYNOMIAL_CHANNEL_FITTER_H
#define POLYNOMIAL_CHANNEL_FITTER_H

#include <cstddef>
#include <vector>

class PolynomialChannelFitter {
 public:
  void Clear() {
    channels_.clear();
    data_points_.clear();
  }

  void AddChannel(double start_frequency, double end_frequency) {
    channels_.emplace_back(start_frequency, end_frequency);
  }

  void AddDataPoint(size_t channel, double y) {
    data_points_.emplace_back(channel, y);
  }

  void Fit(std::vector<double>& terms, size_t n_terms);

  static double Evaluate(double x, const std::vector<double>& terms);

 private:
  /**
   * Start and end frequencies of the channels
   */
  std::vector<std::pair<double, double>> channels_;
  std::vector<std::pair<size_t, double>> data_points_;
};

#endif
