#ifndef GAUSSIAN_FITTER_H
#define GAUSSIAN_FITTER_H

#include <cstring>

class GaussianFitter {
 public:
  void Fit2DGaussianCentred(const float* image, size_t width, size_t height,
                            double beam_est, double& beam_major,
                            double& beam_minor, double& beam_pa,
                            double box_scale_factor = 10.0,
                            bool verbose = false);

  void Fit2DCircularGaussianCentred(const float* image, size_t width,
                                    size_t height, double& beam_size,
                                    double box_scale_factor = 10.0);

  void Fit2DGaussianFull(const float* image, size_t width, size_t height,
                         double& val, double& pos_x, double& pos_y,
                         double& beam_major, double& beam_minor,
                         double& beam_pa, double* floor_level = nullptr);

  const float* Image() const { return image_; }
  size_t Width() const { return width_; }
  size_t Height() const { return height_; }
  size_t ScaleFactor() const { return scale_factor_; }
  double XInit() const { return x_init_; };
  double YInit() const { return y_init_; };
  double PosConstrained() const { return pos_constrained_; };

 private:
  void Fit2DGaussianCentredInBox(const float* image, size_t width,
                                 size_t height, double beam_est,
                                 double& beam_major, double& beam_minor,
                                 double& beam_pa, size_t box_width,
                                 size_t box_height, bool verbose);

  void Fit2DCircularGaussianCentredInBox(const float* image, size_t width,
                                         size_t height, double& beam_size,
                                         size_t box_width, size_t box_height);

  /**
   * This function performs a single fit of a Gaussian. The position of the
   * Gaussian is constrained to be in the centre of the image. The Gaussian is
   * fitted such that the squared residuals (data - model) are minimal.
   *
   * This function is typically used to find the beam-shape of the point-spread
   * function. The beam estimate is used as initial value for the minor and
   * major shape.
   */
  void SingleFit2DGaussianCentred(const float* image, size_t width,
                                  size_t height, double beam_est,
                                  double& beam_major, double& beam_minor,
                                  double& beam_pa, bool verbose);

  void SingleFit2DCircularGaussianCentred(const float* image, size_t width,
                                          size_t height, double& beam_size);

  void Fit2DGaussianWithAmplitudeInBox(const float* image, size_t width,
                                       size_t height, double& val,
                                       double& pos_x, double& pos_y,
                                       double& beam_major, double& beam_minor,
                                       double& beam_pa, double* floor_level,
                                       size_t x_start, size_t x_end,
                                       size_t y_start, size_t y_end);

  /**
   * Fits the position, size and amplitude of a Gaussian. If floor_level is not
   * a nullptr, the floor (background level, or zero level) is fitted too.
   */
  void Fit2DGaussianWithAmplitude(const float* image, size_t width,
                                  size_t height, double& val, double& pos_x,
                                  double& pos_y, double& beam_major,
                                  double& beam_minor, double& beam_pa,
                                  double* floor_level);

  /**
   * Like SingleFit2DGaussianCentred(), but includes Gaussian centre X and Y
   * position and amplitude in the fitted parameters.
   *
   * This function can typically be used for source fitting.
   */
  void Fit2DGaussianWithAmplitude(double& val, double& pos_x, double& pos_y,
                                  double& beam_major, double& beam_minor,
                                  double& beam_pa);

  /**
   * Like Fit2DGaussianWithAmplitude(), but includes floor_level as fitted
   * parameter. Floor is the background/zero level on which the Gaussian
   * resides.
   */
  void Fit2DGaussianWithAmplitudeWithFloor(double& val, double& pos_x,
                                           double& pos_y, double& beam_major,
                                           double& beam_minor, double& beam_pa,
                                           double& floor_level);

  const float* image_ = nullptr;
  size_t width_ = 0;
  size_t height_ = 0;
  size_t scale_factor_ = 0;
  double x_init_ = 0.0;
  double y_init_ = 0.0;
  double pos_constrained_ = 0.0;
};

#endif
