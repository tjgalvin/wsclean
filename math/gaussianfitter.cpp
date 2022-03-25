#include "gaussianfitter.h"

#include <cmath>
#include <iostream>

#include <aocommon/matrix2x2.h>
#include <aocommon/uvector.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>

using aocommon::Matrix2x2;

namespace {

const long double kSigmaToBeam = 2.0L * std::sqrt(2.0L * std::log(2.0L));

void ToAnglesAndFwhm(double sx, double sy, double beta, double& ellipseMaj,
                     double& ellipseMin, double& ellipsePA) {
  const double betaFact = 1.0 - beta * beta;
  double cov[4];
  cov[0] = sx * sx / betaFact;
  cov[1] = beta * sx * sy / betaFact;
  cov[2] = cov[1];
  cov[3] = sy * sy / betaFact;

  double e1, e2, vec1[2], vec2[2];
  Matrix2x2::EigenValuesAndVectors(cov, e1, e2, vec1, vec2);
  if (std::isfinite(e1)) {
    ellipseMaj = std::sqrt(std::fabs(e1)) * kSigmaToBeam;
    ellipseMin = std::sqrt(std::fabs(e2)) * kSigmaToBeam;
    if (ellipseMaj < ellipseMin) {
      std::swap(ellipseMaj, ellipseMin);
      vec1[0] = vec2[0];
      vec1[1] = vec2[1];
    }
    ellipsePA = -std::atan2(vec1[0], vec1[1]);
  } else {
    ellipseMaj = std::sqrt(std::fabs(sx)) * kSigmaToBeam;
    ellipseMin = std::sqrt(std::fabs(sx)) * kSigmaToBeam;
    ellipsePA = 0.0;
  }
}

/**
 * Calculates a two dimensional gaussian with the specified parameters.
 * @param x X coordinate.
 * @param y Y coordinate.
 * @param sx Sigma value for the x direction.
 * @param sy Sigma value for the y direction.
 * @param beta Beta value.
 */
double GaussCentered(double x, double y, double sx, double sy, double beta) {
  return std::exp(-x * x / (2.0 * sx * sx) + beta * x * y / (sx * sy) -
                  y * y / (2.0 * sy * sy));
}

/**
 * Calculates a circular two dimensional gaussian with the specified parameters.
 * @param x X coordinate.
 * @param y Y coordinate.
 * @param s Sigma value for both x and y directions.
 */
double GaussCircularCentered(double x, double y, double s) {
  return std::exp((-x * x - y * y) / (2.0 * s * s));
}

/**
 * Fitting function for SingleFit2DGaussianCentred(). Calculates the sum of the
 * squared errors(/residuals).
 */
int FittingCentered(const gsl_vector* xvec, void* data, gsl_vector* f) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double sx = gsl_vector_get(xvec, 0);
  const double sy = gsl_vector_get(xvec, 1);
  const double beta = gsl_vector_get(xvec, 2);
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int x_mid = width / 2;
  const int y_mid = height / 2;
  const double scale = 1.0 / fitter.ScaleFactor();

  size_t data_index = 0;
  for (size_t yi = 0; yi != height; ++yi) {
    double y = (yi - y_mid) * scale;
    for (size_t xi = 0; xi != width; ++xi) {
      double x = (xi - x_mid) * scale;
      double e = GaussCentered(x, y, sx, sy, beta) - fitter.Image()[data_index];
      gsl_vector_set(f, data_index, e);
      ++data_index;
    }
  }
  return GSL_SUCCESS;
}

int FittingCircularCentered(const gsl_vector* xvec, void* data, gsl_vector* f) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double s = gsl_vector_get(xvec, 0);
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int x_mid = width / 2;
  const int y_mid = height / 2;
  const double scale = 1.0 / fitter.ScaleFactor();

  size_t data_index = 0;
  for (size_t yi = 0; yi != height; ++yi) {
    double y = (yi - y_mid) * scale;
    for (size_t xi = 0; xi != width; ++xi) {
      double x = (xi - x_mid) * scale;
      double e = GaussCircularCentered(x, y, s) - fitter.Image()[data_index];
      gsl_vector_set(f, data_index, e);
      ++data_index;
    }
  }
  return GSL_SUCCESS;
}

/**
 * Derivative function belong with SingleFit2DGaussianCentred().
 */
int FittingDerivativeCentered(const gsl_vector* xvec, void* data,
                              gsl_matrix* J) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double sx = gsl_vector_get(xvec, 0);
  const double sy = gsl_vector_get(xvec, 1);
  const double beta = gsl_vector_get(xvec, 2);
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int x_mid = width / 2;
  const int y_mid = height / 2;
  const double scale = 1.0 / fitter.ScaleFactor();

  size_t data_index = 0;
  for (size_t yi = 0; yi != height; ++yi) {
    double y = (yi - y_mid) * scale;
    for (size_t xi = 0; xi != width; ++xi) {
      double x = (xi - x_mid) * scale;
      double exp_term = GaussCentered(x, y, sx, sy, beta);
      double dsx =
          (beta * x * y / (sx * sx * sy) + x * x / (sx * sx * sx)) * exp_term;
      double dsy =
          (beta * x * y / (sy * sy * sx) + y * y / (sy * sy * sy)) * exp_term;
      double dbeta = x * y / (sx * sy) * exp_term;
      gsl_matrix_set(J, data_index, 0, dsx);
      gsl_matrix_set(J, data_index, 1, dsy);
      gsl_matrix_set(J, data_index, 2, dbeta);
      ++data_index;
    }
  }
  return GSL_SUCCESS;
}

int FittingDerivativeCircularCentered(const gsl_vector* xvec, void* data,
                                      gsl_matrix* J) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double s = gsl_vector_get(xvec, 0);
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int x_mid = width / 2, y_mid = height / 2;
  const double scale = 1.0 / fitter.ScaleFactor();

  size_t data_index = 0;
  for (size_t yi = 0; yi != height; ++yi) {
    double y = (yi - y_mid) * scale;
    for (size_t xi = 0; xi != width; ++xi) {
      double x = (xi - x_mid) * scale;
      double exp_term = GaussCircularCentered(x, y, s);
      // derivative of exp((-x*x - y*y)/(2.0*s*s)) to s
      // = (-x*x - y*y)/2.0*-2/(s*s*s)
      // = (-x*x - y*y)/(-s*s*s)
      // = (x*x + y*y)/(s*s*s)
      double ds = ((x * x + y * y) / (s * s * s)) * exp_term;
      gsl_matrix_set(J, data_index, 0, ds);
      ++data_index;
    }
  }
  return GSL_SUCCESS;
}

/**
 * Squared error and derivative function together.
 */
int FittingBothCentered(const gsl_vector* x, void* data, gsl_vector* f,
                        gsl_matrix* J) {
  FittingCentered(x, data, f);
  FittingDerivativeCentered(x, data, J);
  return GSL_SUCCESS;
}

int FittingBothCircularCentered(const gsl_vector* x, void* data, gsl_vector* f,
                                gsl_matrix* J) {
  FittingCircularCentered(x, data, f);
  FittingDerivativeCircularCentered(x, data, J);
  return GSL_SUCCESS;
}

int FittingWithAmplitude(const gsl_vector* xvec, void* data, gsl_vector* f) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double v = gsl_vector_get(xvec, 0);
  const double xc = gsl_vector_get(xvec, 1);
  const double yc = gsl_vector_get(xvec, 2);
  const double sx = gsl_vector_get(xvec, 3);
  const double sy = gsl_vector_get(xvec, 4);
  const double beta = gsl_vector_get(xvec, 5);
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int x_mid = width / 2;
  const int y_mid = height / 2;
  const double scale = 1.0 / fitter.ScaleFactor();

  size_t data_index = 0;
  for (int yi = 0; yi != int(height); ++yi) {
    double yS = yc + (yi - y_mid) * scale;
    for (int xi = 0; xi != int(width); ++xi) {
      double xS = xc + (xi - x_mid) * scale;
      double e =
          GaussCentered(xS, yS, sx, sy, beta) * v - fitter.Image()[data_index];
      gsl_vector_set(f, data_index, e);
      ++data_index;
    }
  }
  return GSL_SUCCESS;
}

int FittingDerivativeWithAmplitude(const gsl_vector* xvec, void* data,
                                   gsl_matrix* J) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double scale = 1.0 / fitter.ScaleFactor();
  const double v = gsl_vector_get(xvec, 0);
  const double xc = gsl_vector_get(xvec, 1);
  const double yc = gsl_vector_get(xvec, 2);
  const double sx = gsl_vector_get(xvec, 3);
  const double sy = gsl_vector_get(xvec, 4);
  const double beta = gsl_vector_get(xvec, 5);
  if (fitter.PosConstrained() != 0.0 &&
      (std::fabs(xc - fitter.XInit()) > fitter.PosConstrained() * scale ||
       std::fabs(yc - fitter.YInit()) > fitter.PosConstrained() * scale)) {
    std::cout << "GSL_EDOM\n";
    return GSL_EDOM;
  }
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int x_mid = width / 2;
  const int y_mid = height / 2;

  size_t data_index = 0;
  for (int yi = 0; yi != int(height); ++yi) {
    double y = yc + (yi - y_mid) * scale;
    for (int xi = 0; xi != int(width); ++xi) {
      // TODO I need to go over the signs -- ds, dy, dsx, dsy in particular
      double x = xc + (xi - x_mid) * scale;
      double exp_term = GaussCentered(x, y, sx, sy, beta);
      double dv = exp_term;
      exp_term *= v;
      double dx = (-beta * y / (sx * sy) - x / (sx * sx)) * exp_term;
      double dy = (-beta * x / (sy * sx) - y / (sy * sy)) * exp_term;
      double dsx =
          (beta * x * y / (sx * sx * sy) + x * x / (sx * sx * sx)) * exp_term;
      double dsy =
          (beta * x * y / (sy * sy * sx) + y * y / (sy * sy * sy)) * exp_term;
      double dbeta = x * y / (sx * sy) * exp_term;
      gsl_matrix_set(J, data_index, 0, dv);
      gsl_matrix_set(J, data_index, 1, dx);
      gsl_matrix_set(J, data_index, 2, dy);
      gsl_matrix_set(J, data_index, 3, dsx);
      gsl_matrix_set(J, data_index, 4, dsy);
      gsl_matrix_set(J, data_index, 5, dbeta);
      ++data_index;
    }
  }
  return GSL_SUCCESS;
}

int FittingBothWithAmplitude(const gsl_vector* x, void* data, gsl_vector* f,
                             gsl_matrix* J) {
  FittingWithAmplitude(x, data, f);
  FittingDerivativeWithAmplitude(x, data, J);
  return GSL_SUCCESS;
}

int FittingWithAmplitudeAndFloor(const gsl_vector* xvec, void* data,
                                 gsl_vector* f) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double scale = 1.0 / fitter.ScaleFactor();
  const double v = gsl_vector_get(xvec, 0);
  const double xc = gsl_vector_get(xvec, 1);
  const double yc = gsl_vector_get(xvec, 2);
  const double sx = gsl_vector_get(xvec, 3);
  const double sy = gsl_vector_get(xvec, 4);
  const double beta = gsl_vector_get(xvec, 5);
  const double fl = gsl_vector_get(xvec, 6);
  if (fitter.PosConstrained() != 0.0 &&
      (std::fabs(xc - fitter.XInit()) > fitter.PosConstrained() * scale ||
       std::fabs(yc - fitter.YInit()) > fitter.PosConstrained() * scale))
    return GSL_EDOM;
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int x_mid = width / 2;
  const int y_mid = height / 2;

  size_t data_index = 0;
  for (int yi = 0; yi != int(height); ++yi) {
    double yS = yc + (yi - y_mid) * scale;
    for (int xi = 0; xi != int(width); ++xi) {
      double xS = xc + (xi - x_mid) * scale;
      double e = GaussCentered(xS, yS, sx, sy, beta) * v -
                 fitter.Image()[data_index] + fl;
      gsl_vector_set(f, data_index, e);
      ++data_index;
    }
  }
  return GSL_SUCCESS;
}

int FittingDerivativeWithAmplitudeAndFloor(const gsl_vector* xvec, void* data,
                                           gsl_matrix* J) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double v = gsl_vector_get(xvec, 0);
  const double xc = gsl_vector_get(xvec, 1);
  const double yc = gsl_vector_get(xvec, 2);
  const double sx = gsl_vector_get(xvec, 3);
  const double sy = gsl_vector_get(xvec, 4);
  const double beta = gsl_vector_get(xvec, 5);
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int x_mid = width / 2;
  const int y_mid = height / 2;
  const double scale = 1.0 / fitter.ScaleFactor();

  size_t data_index = 0;
  for (int yi = 0; yi != int(height); ++yi) {
    double y = yc + (yi - y_mid) * scale;
    for (int xi = 0; xi != int(width); ++xi) {
      double x = xc + (xi - x_mid) * scale;
      double exp_term = GaussCentered(x, y, sx, sy, beta);
      double dv = exp_term;
      exp_term *= v;
      double dx = (-beta * y / (sx * sy) - x / (sx * sx)) * exp_term;
      double dy = (-beta * x / (sy * sx) - y / (sy * sy)) * exp_term;
      double dsx =
          (beta * x * y / (sx * sx * sy) + x * x / (sx * sx * sx)) * exp_term;
      double dsy =
          (beta * x * y / (sy * sy * sx) + y * y / (sy * sy * sy)) * exp_term;
      double dbeta = x * y / (sx * sy) * exp_term;
      double dfl = 1.0;
      gsl_matrix_set(J, data_index, 0, dv);
      gsl_matrix_set(J, data_index, 1, dx);
      gsl_matrix_set(J, data_index, 2, dy);
      gsl_matrix_set(J, data_index, 3, dsx);
      gsl_matrix_set(J, data_index, 4, dsy);
      gsl_matrix_set(J, data_index, 5, dbeta);
      gsl_matrix_set(J, data_index, 6, dfl);
      ++data_index;
    }
  }
  return GSL_SUCCESS;
}

int FittingBothWithAmplitudeAndFloor(const gsl_vector* x, void* data,
                                     gsl_vector* f, gsl_matrix* J) {
  FittingWithAmplitudeAndFloor(x, data, f);
  FittingDerivativeWithAmplitudeAndFloor(x, data, J);
  return GSL_SUCCESS;
}

}  // namespace

void GaussianFitter::Fit2DGaussianCentred(const float* image, size_t width,
                                          size_t height, double beamEst,
                                          double& beamMaj, double& beamMin,
                                          double& beamPA, double boxScaleFactor,
                                          bool verbose) {
  size_t prefSize = std::max<size_t>(std::ceil(boxScaleFactor),
                                     std::ceil(beamEst * boxScaleFactor));
  if (prefSize % 2 != 0) ++prefSize;
  if (prefSize < width || prefSize < height) {
    size_t nIter = 0;
    bool boxWasLargeEnough;
    do {
      size_t boxWidth = std::min(prefSize, width);
      size_t boxHeight = std::min(prefSize, height);
      if (verbose) std::cout << "Fit initial value:" << beamEst << "\n";
      Fit2DGaussianCentredInBox(image, width, height, beamEst, beamMaj, beamMin,
                                beamPA, boxWidth, boxHeight, verbose);
      if (verbose)
        std::cout << "Fit result:" << beamMaj << " x " << beamMin << " px, "
                  << beamPA << " (box was " << boxWidth << " x " << boxHeight
                  << ")\n";

      boxWasLargeEnough =
          (beamMaj * boxScaleFactor * 0.8 < boxWidth || boxWidth >= width) &&
          (beamMaj * boxScaleFactor * 0.8 < boxHeight || boxHeight >= height);
      if (!boxWasLargeEnough) {
        prefSize = std::max<size_t>(std::ceil(boxScaleFactor),
                                    std::ceil(beamMaj * boxScaleFactor));
        if (prefSize % 2 != 0) ++prefSize;
        beamEst = std::max(beamMaj, beamEst);
      }
      ++nIter;
    } while (!boxWasLargeEnough && nIter < 5);
  } else {
    if (verbose) std::cout << "Image is as large as the fitting box.\n";
    SingleFit2DGaussianCentred(image, width, height, beamEst, beamMaj, beamMin,
                               beamPA, verbose);
  }
}

void GaussianFitter::Fit2DCircularGaussianCentred(const float* image,
                                                  size_t width, size_t height,
                                                  double& beamSize,
                                                  double boxScaleFactor) {
  double initialValue = beamSize;
  size_t prefSize = std::max<size_t>(std::ceil(boxScaleFactor),
                                     std::ceil(beamSize * boxScaleFactor));
  if (prefSize % 2 != 0) ++prefSize;
  if (prefSize < width || prefSize < height) {
    size_t boxWidth = std::min(prefSize, width);
    size_t boxHeight = std::min(prefSize, height);
    size_t nIter = 0;
    bool boxWasLargeEnough;
    do {
      Fit2DCircularGaussianCentredInBox(image, width, height, beamSize,
                                        boxWidth, boxHeight);

      boxWasLargeEnough =
          (beamSize * boxScaleFactor * 0.8 < boxWidth || width >= boxWidth) &&
          (beamSize * boxScaleFactor * 0.8 < boxHeight || height >= boxHeight);
      if (!boxWasLargeEnough) {
        prefSize = std::max<size_t>(std::ceil(boxScaleFactor),
                                    std::ceil(beamSize * boxScaleFactor));
        if (prefSize % 2 != 0) ++prefSize;
        beamSize = std::max(initialValue, beamSize);
      }
      ++nIter;
    } while (!boxWasLargeEnough && nIter < 5);
  } else {
    SingleFit2DCircularGaussianCentred(image, width, height, beamSize);
  }
}

void GaussianFitter::Fit2DGaussianFull(const float* image, size_t width,
                                       size_t height, double& val, double& posX,
                                       double& posY, double& beamMaj,
                                       double& beamMin, double& beamPA,
                                       double* floorLevel) {
  size_t prefSize = std::max<size_t>(10, std::ceil(beamMaj * 10.0));
  if (prefSize % 2 != 0) ++prefSize;
  if (prefSize < width || prefSize < height) {
    size_t xStart = std::max<int>(0, int(std::round(posX)) - int(prefSize) / 2);
    size_t xEnd = std::min(width, size_t(std::round(posX)) + prefSize / 2);
    size_t yStart = std::max<int>(0, int(std::round(posY)) - int(prefSize) / 2);
    size_t yEnd = std::min(height, size_t(std::round(posY)) + prefSize / 2);
    size_t nIter = 0;
    bool boxWasLargeEnough;
    do {
      Fit2DGaussianWithAmplitudeInBox(image, width, height, val, posX, posY,
                                      beamMaj, beamMin, beamPA, floorLevel,
                                      xStart, xEnd, yStart, yEnd);

      size_t boxWidth = xEnd - xStart;
      size_t boxHeight = yEnd - yStart;
      boxWasLargeEnough = (beamMaj * 4.0 < boxWidth || width >= boxWidth) &&
                          (beamMaj * 4.0 < boxHeight || height >= boxHeight);
      if (!boxWasLargeEnough) {
        prefSize = std::max<size_t>(10, std::ceil(beamMaj * 10.0));
        if (prefSize % 2 != 0) ++prefSize;
      }
      ++nIter;
    } while (!boxWasLargeEnough && nIter < 5);
  } else {
    Fit2DGaussianWithAmplitude(image, width, height, val, posX, posY, beamMaj,
                               beamMin, beamPA, floorLevel);
  }
}

void GaussianFitter::Fit2DGaussianCentredInBox(const float* image, size_t width,
                                               size_t height, double beamEst,
                                               double& beamMaj, double& beamMin,
                                               double& beamPA, size_t boxWidth,
                                               size_t boxHeight, bool verbose) {
  size_t startX = (width - boxWidth) / 2;
  size_t startY = (height - boxHeight) / 2;
  aocommon::UVector<float> smallImage(boxWidth * boxHeight);
  for (size_t y = startY; y != (height + boxHeight) / 2; ++y) {
    std::copy_n(&image[y * width + startX], boxWidth,
                &smallImage[(y - startY) * boxWidth]);
  }

  SingleFit2DGaussianCentred(&smallImage[0], boxWidth, boxHeight, beamEst,
                             beamMaj, beamMin, beamPA, verbose);
}

void GaussianFitter::Fit2DCircularGaussianCentredInBox(
    const float* image, size_t width, size_t height, double& beamSize,
    size_t boxWidth, size_t boxHeight) {
  size_t startX = (width - boxWidth) / 2;
  size_t startY = (height - boxHeight) / 2;
  aocommon::UVector<float> smallImage(boxWidth * boxHeight);
  for (size_t y = startY; y != (height + boxHeight) / 2; ++y) {
    std::copy_n(&image[y * width + startX], boxWidth,
                &smallImage[(y - startY) * boxWidth]);
  }

  SingleFit2DCircularGaussianCentred(&smallImage[0], boxWidth, boxHeight,
                                     beamSize);
}

void GaussianFitter::SingleFit2DGaussianCentred(const float* image,
                                                size_t width, size_t height,
                                                double beamEst, double& beamMaj,
                                                double& beamMin, double& beamPA,
                                                bool verbose) {
  width_ = width;
  height_ = height;
  image_ = image;
  scale_factor_ = (width + height) / 2;

  const gsl_multifit_fdfsolver_type* T = gsl_multifit_fdfsolver_lmsder;
  gsl_multifit_fdfsolver* solver =
      gsl_multifit_fdfsolver_alloc(T, width_ * height_, 3);

  gsl_multifit_function_fdf fdf;
  fdf.f = &FittingCentered;
  fdf.df = &FittingDerivativeCentered;
  fdf.fdf = &FittingBothCentered;
  fdf.n = width_ * height_;
  fdf.p = 3;
  fdf.params = this;

  // Using the FWHM formula for a Gaussian:
  double initialValsArray[3] = {
      beamEst / (scale_factor_ * double(kSigmaToBeam)),
      beamEst / (scale_factor_ * double(kSigmaToBeam)), 0.0};
  gsl_vector_view initialVals = gsl_vector_view_array(initialValsArray, 3);
  gsl_multifit_fdfsolver_set(solver, &fdf, &initialVals.vector);

  int status;
  size_t iter = 0;
  do {
    if (verbose) std::cout << "Iteration " << iter << ": ";
    iter++;
    status = gsl_multifit_fdfsolver_iterate(solver);

    if (status) break;

    status = gsl_multifit_test_delta(solver->dx, solver->x, 1e-7, 1e-7);

  } while (status == GSL_CONTINUE && iter < 500);

  double sx = gsl_vector_get(solver->x, 0), sy = gsl_vector_get(solver->x, 1),
         beta = gsl_vector_get(solver->x, 2);

  gsl_multifit_fdfsolver_free(solver);

  ToAnglesAndFwhm(sx, sy, beta, beamMaj, beamMin, beamPA);
  beamMaj *= scale_factor_;
  beamMin *= scale_factor_;
}

void GaussianFitter::SingleFit2DCircularGaussianCentred(const float* image,
                                                        size_t width,
                                                        size_t height,
                                                        double& beamSize) {
  width_ = width;
  height_ = height;
  image_ = image;
  scale_factor_ = (width + height) / 2;

  const gsl_multifit_fdfsolver_type* T = gsl_multifit_fdfsolver_lmsder;
  gsl_multifit_fdfsolver* solver =
      gsl_multifit_fdfsolver_alloc(T, width_ * height_, 1);

  gsl_multifit_function_fdf fdf;
  fdf.f = &FittingCircularCentered;
  fdf.df = &FittingDerivativeCircularCentered;
  fdf.fdf = &FittingBothCircularCentered;
  fdf.n = width_ * height_;
  fdf.p = 1;
  fdf.params = this;

  // Using the FWHM formula for a Gaussian:
  double initialValsArray[1] = {beamSize /
                                (scale_factor_ * double(kSigmaToBeam))};
  gsl_vector_view initialVals = gsl_vector_view_array(initialValsArray, 1);
  gsl_multifit_fdfsolver_set(solver, &fdf, &initialVals.vector);

  int status;
  size_t iter = 0;
  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(solver);

    if (status) break;

    status = gsl_multifit_test_delta(solver->dx, solver->x, 1e-7, 1e-7);

  } while (status == GSL_CONTINUE && iter < 500);

  const double s = gsl_vector_get(solver->x, 0);
  gsl_multifit_fdfsolver_free(solver);

  beamSize = s * kSigmaToBeam * scale_factor_;
}

void GaussianFitter::Fit2DGaussianWithAmplitudeInBox(
    const float* image, size_t width, size_t /*height*/, double& val,
    double& posX, double& posY, double& beamMaj, double& beamMin,
    double& beamPA, double* floorLevel, size_t xStart, size_t xEnd,
    size_t yStart, size_t yEnd) {
  size_t boxWidth = xEnd - xStart;
  size_t boxHeight = yEnd - yStart;
  aocommon::UVector<float> smallImage(boxWidth * boxHeight);
  for (size_t y = yStart; y != yEnd; ++y) {
    std::copy_n(&image[y * width + xStart], boxWidth,
                &smallImage[(y - yStart) * boxWidth]);
  }

  posX -= xStart;
  posY -= yStart;
  Fit2DGaussianWithAmplitude(&smallImage[0], boxWidth, boxHeight, val, posX,
                             posY, beamMaj, beamMin, beamPA, floorLevel);
  posX += xStart;
  posY += yStart;
}

/**
 * Fits the position, size and amplitude of a Gaussian. If floorLevel is not
 * a nullptr, the floor (background level, or zero level) is fitted too.
 */
void GaussianFitter::Fit2DGaussianWithAmplitude(const float* image,
                                                size_t width, size_t height,
                                                double& val, double& posX,
                                                double& posY, double& beamMaj,
                                                double& beamMin, double& beamPA,
                                                double* floorLevel) {
  width_ = width;
  height_ = height;
  image_ = image;
  scale_factor_ = (width + height) / 2;

  if (floorLevel == nullptr)
    Fit2DGaussianWithAmplitude(val, posX, posY, beamMaj, beamMin, beamPA);
  else
    Fit2DGaussianWithAmplitudeWithFloor(val, posX, posY, beamMaj, beamMin,
                                        beamPA, *floorLevel);
}

void GaussianFitter::Fit2DGaussianWithAmplitude(double& val, double& posX,
                                                double& posY, double& beamMaj,
                                                double& beamMin,
                                                double& beamPA) {
  const gsl_multifit_fdfsolver_type* T = gsl_multifit_fdfsolver_lmsder;
  gsl_multifit_fdfsolver* solver =
      gsl_multifit_fdfsolver_alloc(T, width_ * height_, 6);

  gsl_multifit_function_fdf fdf;
  fdf.f = &FittingWithAmplitude;
  fdf.df = &FittingDerivativeWithAmplitude;
  fdf.fdf = &FittingBothWithAmplitude;
  fdf.n = width_ * height_;
  fdf.p = 6;
  fdf.params = this;

  // Using the FWHM formula for a Gaussian:
  x_init_ = -(posX - width_ / 2) / scale_factor_;
  y_init_ = -(posY - height_ / 2) / scale_factor_;
  double initialValsArray[6] = {
      val,
      x_init_,
      y_init_,
      beamMaj / (scale_factor_ * double(kSigmaToBeam)),
      beamMaj / (scale_factor_ * double(kSigmaToBeam)),
      0.0};
  gsl_vector_view initialVals = gsl_vector_view_array(initialValsArray, 6);
  gsl_multifit_fdfsolver_set(solver, &fdf, &initialVals.vector);

  int status;
  size_t iter = 0;
  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(solver);

    if (status) break;

    status = gsl_multifit_test_delta(solver->dx, solver->x, 1e-7, 1e-7);

  } while (status == GSL_CONTINUE && iter < 500);

  val = gsl_vector_get(solver->x, 0);
  posX = -1.0 * gsl_vector_get(solver->x, 1) * scale_factor_ + width_ / 2;
  posY = -1.0 * gsl_vector_get(solver->x, 2) * scale_factor_ + height_ / 2;
  double sx = gsl_vector_get(solver->x, 3), sy = gsl_vector_get(solver->x, 4),
         beta = gsl_vector_get(solver->x, 5);

  gsl_multifit_fdfsolver_free(solver);

  ToAnglesAndFwhm(sx, sy, beta, beamMaj, beamMin, beamPA);
  beamMaj *= scale_factor_;
  beamMin *= scale_factor_;
}

void GaussianFitter::Fit2DGaussianWithAmplitudeWithFloor(
    double& val, double& posX, double& posY, double& beamMaj, double& beamMin,
    double& beamPA, double& floorLevel) {
  const gsl_multifit_fdfsolver_type* T = gsl_multifit_fdfsolver_lmsder;
  gsl_multifit_fdfsolver* solver =
      gsl_multifit_fdfsolver_alloc(T, width_ * height_, 7);

  gsl_multifit_function_fdf fdf;
  fdf.f = &FittingWithAmplitudeAndFloor;
  fdf.df = &FittingDerivativeWithAmplitudeAndFloor;
  fdf.fdf = &FittingBothWithAmplitudeAndFloor;
  fdf.n = width_ * height_;
  fdf.p = 7;
  fdf.params = this;

  // Using the FWHM formula for a Gaussian:
  x_init_ = -(posX - width_ / 2) / scale_factor_;
  y_init_ = -(posY - height_ / 2) / scale_factor_;
  double initialValsArray[7] = {
      val,
      x_init_,
      y_init_,
      beamMaj / (scale_factor_ * double(kSigmaToBeam)),
      beamMaj / (scale_factor_ * double(kSigmaToBeam)),
      0.0,
      0.0};
  gsl_vector_view initialVals = gsl_vector_view_array(initialValsArray, 7);
  gsl_multifit_fdfsolver_set(solver, &fdf, &initialVals.vector);

  int status;
  size_t iter = 0;
  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(solver);

    if (status) break;

    status = gsl_multifit_test_delta(solver->dx, solver->x, 1e-7, 1e-7);

  } while (status == GSL_CONTINUE && iter < 500);

  val = gsl_vector_get(solver->x, 0);
  posX = -1.0 * gsl_vector_get(solver->x, 1) * scale_factor_ + width_ / 2;
  posY = -1.0 * gsl_vector_get(solver->x, 2) * scale_factor_ + height_ / 2;
  double sx = gsl_vector_get(solver->x, 3);
  double sy = gsl_vector_get(solver->x, 4);
  double beta = gsl_vector_get(solver->x, 5);
  floorLevel = gsl_vector_get(solver->x, 6);

  gsl_multifit_fdfsolver_free(solver);

  ToAnglesAndFwhm(sx, sy, beta, beamMaj, beamMin, beamPA);
  beamMaj *= scale_factor_;
  beamMin *= scale_factor_;
}
