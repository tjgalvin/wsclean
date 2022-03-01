#include "peakfinder.h"

#ifdef __SSE__
#define USE_INTRINSICS
#endif

#ifdef USE_INTRINSICS
#include <emmintrin.h>
#include <immintrin.h>
#endif

#include <limits>

std::optional<float> PeakFinder::Simple(const float *image, size_t width,
                                        size_t height, size_t &x, size_t &y,
                                        bool allowNegativeComponents,
                                        size_t startY, size_t endY,
                                        size_t horizontalBorder,
                                        size_t verticalBorder) {
  float peakMax = std::numeric_limits<float>::min();
  size_t peakIndex = width * height;

  size_t xiStart = horizontalBorder, xiEnd = width - horizontalBorder;
  size_t yiStart = std::max(startY, verticalBorder),
         yiEnd = std::min(endY, height - verticalBorder);
  if (xiEnd < xiStart) xiEnd = xiStart;
  if (yiEnd < yiStart) yiEnd = yiStart;

  for (size_t yi = yiStart; yi != yiEnd; ++yi) {
    size_t index = yi * width + xiStart;
    for (size_t xi = xiStart; xi != xiEnd; ++xi) {
      float value = image[index];
      if (allowNegativeComponents) value = std::fabs(value);
      if (value > peakMax) {
        peakIndex = index;
        peakMax = std::fabs(value);
      }
      ++value;
      ++index;
    }
  }
  if (peakIndex == width * height) {
    x = width;
    y = height;
    return std::optional<float>();
  } else {
    x = peakIndex % width;
    y = peakIndex / width;
    return image[x + y * width];
  }
}

std::optional<double> PeakFinder::Simple(const double *image, size_t width,
                                         size_t height, size_t &x, size_t &y,
                                         bool allowNegativeComponents,
                                         size_t startY, size_t endY,
                                         size_t horizontalBorder,
                                         size_t verticalBorder) {
  double peakMax = std::numeric_limits<double>::min();
  size_t peakIndex = width * height;

  size_t xiStart = horizontalBorder, xiEnd = width - horizontalBorder;
  size_t yiStart = std::max(startY, verticalBorder),
         yiEnd = std::min(endY, height - verticalBorder);
  if (xiEnd < xiStart) xiEnd = xiStart;
  if (yiEnd < yiStart) yiEnd = yiStart;

  for (size_t yi = yiStart; yi != yiEnd; ++yi) {
    size_t index = yi * width + xiStart;
    for (size_t xi = xiStart; xi != xiEnd; ++xi) {
      double value = image[index];
      if (allowNegativeComponents) value = std::fabs(value);
      if (value > peakMax) {
        peakIndex = index;
        peakMax = std::fabs(value);
      }
      ++value;
      ++index;
    }
  }
  if (peakIndex == width * height) {
    x = width;
    y = height;
    return {};
  } else {
    x = peakIndex % width;
    y = peakIndex / width;
    return image[x + y * width];
  }
}

std::optional<float> PeakFinder::FindWithMask(
    const float *image, size_t width, size_t height, size_t &x, size_t &y,
    bool allowNegativeComponents, size_t startY, size_t endY,
    const bool *cleanMask, size_t horizontalBorder, size_t verticalBorder) {
  float peakMax = std::numeric_limits<float>::min();
  x = width;
  y = height;

  size_t xiStart = horizontalBorder, xiEnd = width - horizontalBorder;
  size_t yiStart = std::max(startY, verticalBorder),
         yiEnd = std::min(endY, height - verticalBorder);
  if (xiEnd < xiStart) xiEnd = xiStart;
  if (yiEnd < yiStart) yiEnd = yiStart;

  for (size_t yi = yiStart; yi != yiEnd; ++yi) {
    const float *imgIter = &image[yi * width + xiStart];
    const bool *cleanMaskPtr = &cleanMask[yi * width + xiStart];
    for (size_t xi = xiStart; xi != xiEnd; ++xi) {
      float value = *imgIter;
      if (allowNegativeComponents) value = std::fabs(value);
      if (value > peakMax && *cleanMaskPtr) {
        x = xi;
        y = yi;
        peakMax = std::fabs(value);
      }
      ++imgIter;
      ++cleanMaskPtr;
    }
  }
  if (y == height)
    return std::optional<float>();
  else
    return image[x + y * width];
}

#if defined __AVX__ && defined USE_INTRINSICS && !defined FORCE_NON_AVX
template <bool AllowNegativeComponent>
std::optional<double> PeakFinder::AVX(const double *image, size_t width,
                                      size_t height, size_t &x, size_t &y,
                                      size_t startY, size_t endY,
                                      size_t horizontalBorder,
                                      size_t verticalBorder) {
  double peakMax = std::numeric_limits<double>::min();
  size_t peakIndex = 0;

  __m256d mPeakMax = _mm256_set1_pd(peakMax);

  size_t xiStart = horizontalBorder, xiEnd = width - horizontalBorder;
  size_t yiStart = std::max(startY, verticalBorder),
         yiEnd = std::min(endY, height - verticalBorder);
  if (xiEnd < xiStart) xiEnd = xiStart;
  if (yiEnd < yiStart) yiEnd = yiStart;

  for (size_t yi = yiStart; yi != yiEnd; ++yi) {
    size_t index = yi * width + xiStart;
    const double *const endPtr = image + yi * width + xiEnd - 4;
    const double *i = image + index;
    for (; i < endPtr; i += 4) {
      __m256d val = _mm256_loadu_pd(i);
      if (AllowNegativeComponent) {
        __m256d negVal = _mm256_sub_pd(_mm256_set1_pd(0.0), val);
        val = _mm256_max_pd(val, negVal);
      }
      int mask = _mm256_movemask_pd(_mm256_cmp_pd(val, mPeakMax, _CMP_GT_OQ));
      if (mask != 0) {
        for (size_t di = 0; di != 4; ++di) {
          double value = i[di];
          if (AllowNegativeComponent) value = std::fabs(value);
          if (value > peakMax) {
            peakIndex = index + di;
            peakMax = std::fabs(i[di]);
            mPeakMax = _mm256_set1_pd(peakMax);
          }
        }
      }
      index += 4;
    }
    for (; i != endPtr + 4; ++i) {
      double value = *i;
      if (AllowNegativeComponent) value = std::fabs(value);
      if (value > peakMax) {
        peakIndex = index;
        peakMax = std::fabs(*i);
      }
      ++index;
    }
  }
  x = peakIndex % width;
  y = peakIndex / width;
  return image[x + y * width];
}

template std::optional<double> PeakFinder::AVX<false>(
    const double *image, size_t width, size_t height, size_t &x, size_t &y,
    size_t startY, size_t endY, size_t horizontalBorder, size_t verticalBorder);
template std::optional<double> PeakFinder::AVX<true>(
    const double *image, size_t width, size_t height, size_t &x, size_t &y,
    size_t startY, size_t endY, size_t horizontalBorder, size_t verticalBorder);

template <bool AllowNegativeComponent>
std::optional<float> PeakFinder::AVX(const float *image, size_t width,
                                     size_t height, size_t &x, size_t &y,
                                     size_t startY, size_t endY,
                                     size_t horizontalBorder,
                                     size_t verticalBorder) {
  float peakMax = std::numeric_limits<float>::min();
  size_t peakIndex = 0;

  __m256 mPeakMax = _mm256_set1_ps(peakMax);

  size_t xiStart = horizontalBorder, xiEnd = width - horizontalBorder;
  size_t yiStart = std::max(startY, verticalBorder),
         yiEnd = std::min(endY, height - verticalBorder);
  if (xiEnd < xiStart) xiEnd = xiStart;
  if (yiEnd < yiStart) yiEnd = yiStart;

  for (size_t yi = yiStart; yi != yiEnd; ++yi) {
    size_t index = yi * width + xiStart;
    const float *const endPtr = image + yi * width + xiEnd - 8;
    const float *i = image + index;
    for (; i < endPtr; i += 8) {
      __m256 val = _mm256_loadu_ps(i);
      if (AllowNegativeComponent) {
        __m256 negVal = _mm256_sub_ps(_mm256_set1_ps(0.0), val);
        val = _mm256_max_ps(val, negVal);
      }
      int mask = _mm256_movemask_ps(_mm256_cmp_ps(val, mPeakMax, _CMP_GT_OQ));
      if (mask != 0) {
        for (size_t di = 0; di != 8; ++di) {
          double value = i[di];
          if (AllowNegativeComponent) value = std::fabs(value);
          if (value > peakMax) {
            peakIndex = index + di;
            peakMax = std::fabs(i[di]);
            mPeakMax = _mm256_set1_ps(peakMax);
          }
        }
      }
      index += 8;
    }
    for (; i != endPtr + 8; ++i) {
      double value = *i;
      if (AllowNegativeComponent) value = std::fabs(value);
      if (value > peakMax) {
        peakIndex = index;
        peakMax = std::fabs(*i);
      }
      ++index;
    }
  }
  x = peakIndex % width;
  y = peakIndex / width;
  return image[x + y * width];
}

template std::optional<float> PeakFinder::AVX<false>(
    const float *image, size_t width, size_t height, size_t &x, size_t &y,
    size_t startY, size_t endY, size_t horizontalBorder, size_t verticalBorder);
template std::optional<float> PeakFinder::AVX<true>(
    const float *image, size_t width, size_t height, size_t &x, size_t &y,
    size_t startY, size_t endY, size_t horizontalBorder, size_t verticalBorder);

#else
#warning "Not using AVX optimized version of FindPeak()!"
#endif  // __AVX__
