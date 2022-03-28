#include "polynomialchannelfitter.h"

#include <cassert>
#include <stdexcept>

#include "gsl/gsl_multifit.h"

void PolynomialChannelFitter::Fit(std::vector<double>& terms, size_t n_terms) {
  const size_t n_points = data_points_.size();

  gsl_multifit_linear_workspace* work =
      gsl_multifit_linear_alloc(n_points, n_terms);

  gsl_matrix* x = gsl_matrix_alloc(n_points, n_terms);
  gsl_vector* y = gsl_vector_alloc(n_points);
  gsl_vector* c = gsl_vector_calloc(n_terms);
  gsl_matrix* cov = gsl_matrix_calloc(n_terms, n_terms);
  double chisq;

  for (size_t i = 0; i != n_points; ++i) {
    size_t channel_index = data_points_[i].first;
    const double nu_begin = channels_[channel_index].first;
    const double nu_end = channels_[channel_index].second;
    double nu_begin_term = nu_begin;
    double nu_end_term = nu_end;
    gsl_matrix_set(x, i, 0, 1.0);
    for (size_t j = 1; j != n_terms; ++j) {
      nu_begin_term *= nu_begin;
      nu_end_term *= nu_end;
      // ( mu_i^t - nu_i^t ) / (t [ mu_i - nu_i ] ) with t=j+1
      const double value =
          (nu_end_term - nu_begin_term) / (double(j + 1) * (nu_end - nu_begin));
      gsl_matrix_set(x, i, j, value);
    }
    gsl_vector_set(y, i, data_points_[i].second);
  }

  const int result = gsl_multifit_linear(x, y, c, cov, &chisq, work);
  if (result == 0) {
    terms.resize(n_terms);
    for (size_t j = 0; j != n_terms; ++j) {
      terms[j] = gsl_vector_get(c, j);
    }
  }

  gsl_matrix_free(x);
  gsl_vector_free(c);
  gsl_vector_free(y);

  gsl_matrix_free(cov);

  gsl_multifit_linear_free(work);

  if (result != 0) throw std::runtime_error("Linear fit failed");
}

double PolynomialChannelFitter::Evaluate(double x,
                                         const std::vector<double>& terms) {
  assert(!terms.empty());
  double value = terms.front();
  double f = 1.0;
  for (size_t i = 1; i != terms.size(); ++i) {
    f *= x;
    value += f * terms[i];
  }
  return value;
}