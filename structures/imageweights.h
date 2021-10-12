#ifndef IMAGE_WEIGHTS_H
#define IMAGE_WEIGHTS_H

#include <cstddef>
#include <complex>

#include <aocommon/io/serialstreamfwd.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include "msselection.h"
#include "weightmode.h"

class ImageWeights {
 public:
  /**
   * Constructor to be used when needing to unserialize an object
   */
  ImageWeights();

  ImageWeights(const WeightMode& weightMode, size_t imageWidth,
               size_t imageHeight, double pixelScaleX, double pixelScaleY,
               bool weightsAsTaper, double superWeight, size_t threadCount = 1);

  ImageWeights(const ImageWeights&) = delete;
  ImageWeights& operator=(const ImageWeights&) = delete;

  double GetWeight(double u, double v) const { return sampleGridValue(u, v); }

  void Grid(class MSProvider& ms, const MSSelection& selection);
  void Grid(double u, double v, double weight) {
    int x, y;
    uvToXY(u, v, x, y);

    if (isWithinLimits(x, y)) {
      size_t index = (size_t)x + (size_t)y * _imageWidth;
      _grid[index] += weight;
      _totalSum += weight;
    }
  }

  void FinishGridding();

  void SetMaxUVRange(double maxUVInLambda);
  void SetMinUVRange(double minUVInLambda);
  void SetGaussianTaper(double beamSize);
  void SetTukeyTaper(double transitionSizeInLambda, double maxUVInLambda);
  void SetTukeyInnerTaper(double transitionSizeInLambda, double minUVInLambda);
  void SetEdgeTaper(double sizeInLambda);
  void SetEdgeTukeyTaper(double transitionSizeInLambda,
                         double edgeSizeInLambda);
  void SetThreadCount(size_t threadCount) { _threadCount = threadCount; }

  void SetAllValues(double newValue) { _grid.assign(_grid.size(), newValue); }
  void GetGrid(double* image) const;
  void Save(const std::string& filename) const;

  void RankFilter(double rankLimit, size_t windowSize);

  size_t Width() const { return _imageWidth; }
  size_t Height() const { return _imageHeight; }
  WeightMode GetWeightMode() const { return _weightMode; }

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);

 private:
  void uvToXY(double u, double v, int& x, int& y) const {
    if (v < 0.0) {
      u = -u;
      v = -v;
    }
    x = int(floor(u * _imageWidth * _pixelScaleX + _imageWidth / 2));
    y = int(floor(v * _imageHeight * _pixelScaleY));
  }

  void xyToUV(int x, int y, double& u, double& v) const {
    if (y < 0.0) {
      x = -x;
      y = -y;
    }

    u = double(x - int(_imageWidth / 2)) / double(_imageWidth * _pixelScaleX);
    v = double(y) / double(_imageHeight * _pixelScaleY);
  }

  bool isWithinLimits(int x, int y) const {
    return x >= 0 && x < int(_imageWidth) && y < int(_imageHeight / 2);
  }

  double sampleGridValue(double u, double v) const {
    int x, y;
    uvToXY(u, v, x, y);
    if (isWithinLimits(x, y))
      return _grid[(size_t)x + (size_t)y * _imageWidth];
    else {
      return 0.0;
    }
  }

  double windowMean(size_t x, size_t y, size_t windowSize);

  /**
   * Returns Tukey tapering function. This function is
   * 0 when x=0 and 1 when x=n.
   */
  double tukeyFrom0ToN(double x, double n) {
    return 0.5 * (1.0 + cos((M_PI / n) * (x - n)));
  }

  template <typename T>
  static T frequencyToWavelength(const T frequency) {
    return speedOfLight() / frequency;
  }
  static long double speedOfLight() { return 299792458.0L; }
  WeightMode _weightMode;
  std::size_t _imageWidth, _imageHeight;
  double _pixelScaleX, _pixelScaleY;

  std::vector<double> _grid;
  double _totalSum;
  bool _isGriddingFinished, _weightsAsTaper;
  size_t _threadCount;
};

#endif
