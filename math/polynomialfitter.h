#ifndef POLYNOMIAL_FITTER_H
#define POLYNOMIAL_FITTER_H

#include <array>
#include <vector>

class PolynomialFitter {
 public:
  typedef float num_t;

  void Clear() { data_points_.clear(); }

  void AddDataPoint(num_t x, num_t y, num_t w) {
    data_points_.emplace_back(std::array<num_t, 3>{{x, y, w}});
  }

  /**
   * @param [out] terms The resulting terms.
   * Using a pre-allocated vector instead of a return value avoids
   * memory allocations in this performance-critical function.
   */
  void Fit(std::vector<num_t>& terms, size_t nTerms);

  static num_t Evaluate(num_t x, const std::vector<num_t>& terms) {
    num_t val = terms[0];
    num_t f = 1.0;
    for (size_t i = 1; i != terms.size(); ++i) {
      f *= x;
      val += f * terms[i];
    }
    return val;
  }

  size_t Size() const { return data_points_.size(); }

 private:
  std::vector<std::array<num_t, 3>> data_points_;
};

#endif
