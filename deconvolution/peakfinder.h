#ifndef PEAK_FINDER_H
#define PEAK_FINDER_H

#include <cmath>
#include <cstring>
#include <optional>

#ifdef __SSE__
#define USE_INTRINSICS
#endif

class PeakFinder {
 public:
  PeakFinder() = delete;

  static std::optional<float> Simple(const float *image, size_t width,
                                     size_t height, size_t &x, size_t &y,
                                     bool allowNegativeComponents,
                                     size_t startY, size_t endY,
                                     size_t horizontalBorder,
                                     size_t verticalBorder);

  static std::optional<double> Simple(const double *image, size_t width,
                                      size_t height, size_t &x, size_t &y,
                                      bool allowNegativeComponents,
                                      size_t startY, size_t endY,
                                      size_t horizontalBorder,
                                      size_t verticalBorder);

#if defined __AVX__ && defined USE_INTRINSICS && !defined FORCE_NON_AVX
  template <bool AllowNegativeComponent>
  static std::optional<float> AVX(const float *image, size_t width,
                                  size_t height, size_t &x, size_t &y,
                                  size_t startY, size_t endY,
                                  size_t horizontalBorder,
                                  size_t verticalBorder);

  static std::optional<float> AVX(const float *image, size_t width,
                                  size_t height, size_t &x, size_t &y,
                                  bool allowNegativeComponents, size_t startY,
                                  size_t endY, size_t horizontalBorder,
                                  size_t verticalBorder) {
    if (allowNegativeComponents)
      return AVX<true>(image, width, height, x, y, startY, endY,
                       horizontalBorder, verticalBorder);
    else
      return AVX<false>(image, width, height, x, y, startY, endY,
                        horizontalBorder, verticalBorder);
  }

  template <bool AllowNegativeComponent>
  static std::optional<double> AVX(const double *image, size_t width,
                                   size_t height, size_t &x, size_t &y,
                                   size_t startY, size_t endY,
                                   size_t horizontalBorder,
                                   size_t verticalBorder);

  static std::optional<double> AVX(const double *image, size_t width,
                                   size_t height, size_t &x, size_t &y,
                                   bool allowNegativeComponents, size_t startY,
                                   size_t endY, size_t horizontalBorder,
                                   size_t verticalBorder) {
    if (allowNegativeComponents)
      return AVX<true>(image, width, height, x, y, startY, endY,
                       horizontalBorder, verticalBorder);
    else
      return AVX<false>(image, width, height, x, y, startY, endY,
                        horizontalBorder, verticalBorder);
  }
#endif

  /**
   * Find peaks with a relative border ratio.
   */
  template <typename NumT>
  static std::optional<NumT> Find(const NumT *image, size_t width,
                                  size_t height, size_t &x, size_t &y,
                                  bool allowNegativeComponents, size_t startY,
                                  size_t endY, float borderRatio) {
    return Find(image, width, height, x, y, allowNegativeComponents, startY,
                endY, round(width * borderRatio), round(height * borderRatio));
  }

  /**
   * Find peaks with a fixed border.
   */
  template <typename NumT>
  static std::optional<NumT> Find(const NumT *image, size_t width,
                                  size_t height, size_t &x, size_t &y,
                                  bool allowNegativeComponents, size_t startY,
                                  size_t endY, size_t horizontalBorder,
                                  size_t verticalBorder) {
#if defined __AVX__ && defined USE_INTRINSICS && !defined FORCE_NON_AVX
    return AVX(image, width, height, x, y, allowNegativeComponents, startY,
               endY, horizontalBorder, verticalBorder);
#else
    return Simple(image, width, height, x, y, allowNegativeComponents, startY,
                  endY, horizontalBorder, verticalBorder);
#endif
  }

  static std::optional<float> FindWithMask(const float *image, size_t width,
                                           size_t height, size_t &x, size_t &y,
                                           bool allowNegativeComponents,
                                           const bool *cleanMask);

  static std::optional<float> FindWithMask(const float *image, size_t width,
                                           size_t height, size_t &x, size_t &y,
                                           bool allowNegativeComponents,
                                           size_t startY, size_t endY,
                                           const bool *cleanMask,
                                           float borderRatio) {
    return FindWithMask(image, width, height, x, y, allowNegativeComponents,
                        startY, endY, cleanMask, round(width * borderRatio),
                        round(height * borderRatio));
  }

  static std::optional<float> FindWithMask(
      const float *image, size_t width, size_t height, size_t &x, size_t &y,
      bool allowNegativeComponents, size_t startY, size_t endY,
      const bool *cleanMask, size_t horizontalBorder, size_t verticalBorder);
};

#endif
