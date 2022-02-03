#ifndef MULTI_SCALE_TRANSFORMS_H
#define MULTI_SCALE_TRANSFORMS_H

#include <cmath>
#include <initializer_list>

#include <aocommon/image.h>
#include <aocommon/uvector.h>

#include <vector>

class MultiScaleTransforms {
 public:
  enum Shape { TaperedQuadraticShape, GaussianShape };

  MultiScaleTransforms(class FFTWManager& fftwManager, size_t width,
                       size_t height, Shape shape)
      : _fftwManager(fftwManager),
        _width(width),
        _height(height),
        _shape(shape),
        _threadCount(1) {}

  void PrepareTransform(float* kernel, float scale);
  void FinishTransform(float* image, const float* kernel);

  void Transform(aocommon::Image& image, aocommon::Image& scratch,
                 float scale) {
    std::vector<aocommon::Image> images(1, std::move(image));
    Transform(images, scratch, scale);
    image = std::move(images[0]);
  }

  void Transform(std::vector<aocommon::Image>& images, aocommon::Image& scratch,
                 float scale);

  size_t Width() const { return _width; }
  size_t Height() const { return _height; }

  static float KernelIntegratedValue(float scaleInPixels, size_t maxN,
                                     Shape shape) {
    size_t n;
    aocommon::Image kernel = MakeShapeFunction(scaleInPixels, n, maxN, shape);

    float value = 0.0;
    for (float& x : kernel) value += x;

    return value;
  }

  static float KernelPeakValue(double scaleInPixels, size_t maxN, Shape shape) {
    size_t n;
    aocommon::Image kernel = MakeShapeFunction(scaleInPixels, n, maxN, shape);
    return kernel[n / 2 + (n / 2) * n];
  }

  static void AddShapeComponent(float* image, size_t width, size_t height,
                                float scaleSizeInPixels, size_t x, size_t y,
                                float gain, Shape shape) {
    size_t n;
    aocommon::Image kernel =
        MakeShapeFunction(scaleSizeInPixels, n, std::min(width, height), shape);
    int left;
    if (x > n / 2)
      left = x - n / 2;
    else
      left = 0;
    int top;
    if (y > n / 2)
      top = y - n / 2;
    else
      top = 0;
    size_t right = std::min(x + (n + 1) / 2, width);
    size_t bottom = std::min(y + (n + 1) / 2, height);
    for (size_t yi = top; yi != bottom; ++yi) {
      float* imagePtr = &image[yi * width];
      const float* kernelPtr =
          &kernel.Data()[(yi + n / 2 - y) * n + left + n / 2 - x];
      for (size_t xi = left; xi != right; ++xi) {
        imagePtr[xi] += *kernelPtr * gain;
        ++kernelPtr;
      }
    }
  }

  static aocommon::Image MakeShapeFunction(float scaleSizeInPixels, size_t& n,
                                           size_t maxN, Shape shape) {
    switch (shape) {
      default:
      case TaperedQuadraticShape:
        return makeTaperedQuadraticShapeFunction(scaleSizeInPixels, n);
      case GaussianShape:
        return makeGaussianFunction(scaleSizeInPixels, n, maxN);
    }
  }

  aocommon::Image MakeShapeFunction(float scaleSizeInPixels, size_t& n) {
    return MakeShapeFunction(scaleSizeInPixels, n, std::min(_width, _height),
                             _shape);
  }

  static float GaussianSigma(float scaleSizeInPixels) {
    return scaleSizeInPixels * (3.0 / 16.0);
  }

  void SetThreadCount(size_t threadCount) { _threadCount = threadCount; }

 private:
  class FFTWManager& _fftwManager;
  size_t _width, _height;
  enum Shape _shape;
  size_t _threadCount;

  static size_t taperedQuadraticKernelSize(double scaleInPixels) {
    return size_t(ceil(scaleInPixels * 0.5) * 2.0) + 1;
  }

  static aocommon::Image makeTaperedQuadraticShapeFunction(
      double scaleSizeInPixels, size_t& n) {
    n = taperedQuadraticKernelSize(scaleSizeInPixels);
    aocommon::Image output(n, n);
    taperedQuadraticShapeFunction(n, output, scaleSizeInPixels);
    return output;
  }

  static aocommon::Image makeGaussianFunction(double scaleSizeInPixels,
                                              size_t& n, size_t maxN) {
    float sigma = GaussianSigma(scaleSizeInPixels);

    n = int(ceil(sigma * 12.0 / 2.0)) * 2 + 1;  // bounding box of 12 sigma
    if (n > maxN) {
      n = maxN;
      if ((n % 2) == 0 && n > 0) --n;
    }
    if (n < 1) n = 1;
    if (sigma == 0.0) {
      sigma = 1.0;
      n = 1;
    }
    aocommon::Image output(n, n);
    const float mu = int(n / 2);
    const float twoSigmaSquared = 2.0 * sigma * sigma;
    float sum = 0.0;
    float* outputPtr = output.Data();
    aocommon::UVector<float> gaus(n);
    for (int i = 0; i != int(n); ++i) {
      float vI = float(i) - mu;
      gaus[i] = std::exp(-vI * vI / twoSigmaSquared);
    }
    for (int y = 0; y != int(n); ++y) {
      for (int x = 0; x != int(n); ++x) {
        *outputPtr = gaus[x] * gaus[y];
        sum += *outputPtr;
        ++outputPtr;
      }
    }
    float normFactor = 1.0 / sum;
    output *= normFactor;
    return output;
  }

  static void taperedQuadraticShapeFunction(size_t n, aocommon::Image& output2d,
                                            double scaleSizeInPixels) {
    if (scaleSizeInPixels == 0.0)
      output2d[0] = 1.0;
    else {
      float sum = 0.0;
      float* outputPtr = output2d.Data();
      for (int y = 0; y != int(n); ++y) {
        float dy = y - 0.5 * (n - 1);
        float dydy = dy * dy;
        for (int x = 0; x != int(n); ++x) {
          float dx = x - 0.5 * (n - 1);
          float r = std::sqrt(dx * dx + dydy);
          *outputPtr =
              hannWindowFunction(r, n) * shapeFunction(r / scaleSizeInPixels);
          sum += *outputPtr;
          ++outputPtr;
        }
      }
      float normFactor = 1.0 / sum;
      output2d *= normFactor;
    }
  }

  static float hannWindowFunction(float x, size_t n) {
    return (x * 2 <= n + 1)
               ? (0.5 * (1.0 + std::cos(2.0 * M_PI * x / double(n + 1))))
               : 0.0;
  }

  static float shapeFunction(float x) {
    return (x < 1.0) ? (1.0 - x * x) : 0.0;
  }
};

#endif
