#ifndef PEAK_FINDER_H
#define PEAK_FINDER_H

#include <boost/optional/optional.hpp>

#include <cmath>
#include <cstring>

#ifdef __SSE__
#define USE_INTRINSICS
#endif

class PeakFinder {
 public:
  PeakFinder() = delete;

  static boost::optional<float> Simple(const float *image, size_t width,
                                       size_t height, size_t &x, size_t &y,
                                       bool allowNegativeComponents,
                                       size_t startY, size_t endY,
                                       size_t horizontalBorder,
                                       size_t verticalBorder);

  static boost::optional<double> Simple(const double *image, size_t width,
                                        size_t height, size_t &x, size_t &y,
                                        bool allowNegativeComponents,
                                        size_t startY, size_t endY,
                                        size_t horizontalBorder,
                                        size_t verticalBorder);

#if defined __AVX__ && defined USE_INTRINSICS && !defined FORCE_NON_AVX
  template <bool AllowNegativeComponent>
  static boost::optional<float> AVX(const float *image, size_t width,
                                    size_t height, size_t &x, size_t &y,
                                    size_t startY, size_t endY,
                                    size_t horizontalBorder,
                                    size_t verticalBorder);

  static boost::optional<float> AVX(const float *image, size_t width,
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
  static boost::optional<double> AVX(const double *image, size_t width,
                                     size_t height, size_t &x, size_t &y,
                                     size_t startY, size_t endY,
                                     size_t horizontalBorder,
                                     size_t verticalBorder);

  static boost::optional<double> AVX(const double *image, size_t width,
                                     size_t height, size_t &x, size_t &y,
                                     bool allowNegativeComponents,
                                     size_t startY, size_t endY,
                                     size_t horizontalBorder,
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
  static boost::optional<NumT> Find(const NumT *image, size_t width,
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
  static boost::optional<NumT> Find(const NumT *image, size_t width,
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

  static boost::optional<float> FindWithMask(const float *image, size_t width,
                                             size_t height, size_t &x,
                                             size_t &y,
                                             bool allowNegativeComponents,
                                             const bool *cleanMask);

  static boost::optional<float> FindWithMask(
      const float *image, size_t width, size_t height, size_t &x, size_t &y,
      bool allowNegativeComponents, size_t startY, size_t endY,
      const bool *cleanMask, float borderRatio) {
    return FindWithMask(image, width, height, x, y, allowNegativeComponents,
                        startY, endY, cleanMask, round(width * borderRatio),
                        round(height * borderRatio));
  }

  static boost::optional<float> FindWithMask(
      const float *image, size_t width, size_t height, size_t &x, size_t &y,
      bool allowNegativeComponents, size_t startY, size_t endY,
      const bool *cleanMask, size_t horizontalBorder, size_t verticalBorder);
};

#endif
