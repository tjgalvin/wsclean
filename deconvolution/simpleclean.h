#ifndef SIMPLE_CLEAN_H
#define SIMPLE_CLEAN_H

#include "deconvolutionalgorithm.h"

#include <cstring>

#ifdef __SSE__
#define USE_INTRINSICS
#endif

class SimpleClean {
 public:
  SimpleClean() = delete;
  static void SubtractImage(float *image, const float *psf, size_t width,
                            size_t height, size_t x, size_t y, float factor);

  static void PartialSubtractImage(float *image, const float *psf, size_t width,
                                   size_t height, size_t x, size_t y,
                                   float factor, size_t startY, size_t endY);

  static void PartialSubtractImage(float *image, size_t imgWidth,
                                   size_t imgHeight, const float *psf,
                                   size_t psfWidth, size_t psfHeight, size_t x,
                                   size_t y, float factor, size_t startY,
                                   size_t endY);

#if defined __AVX__ && defined USE_INTRINSICS
  static void PartialSubtractImageAVX(double *image, size_t imgWidth,
                                      size_t imgHeight, const double *psf,
                                      size_t psfWidth, size_t psfHeight,
                                      size_t x, size_t y, double factor,
                                      size_t startY, size_t endY);
#endif
};

#endif
