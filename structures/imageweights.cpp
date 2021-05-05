#include "imageweights.h"

#include "multibanddata.h"

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

#include "../msproviders/msprovider.h"
#include "../msproviders/msreaders/msreader.h"

#include "../units/angle.h"

#include "../io/logger.h"

#include <aocommon/fits/fitswriter.h>
#include <aocommon/banddata.h>

#include <cmath>
#include <iostream>
#include <cstring>

ImageWeights::ImageWeights()
    : _weightMode(WeightMode::UniformWeighted),
      _imageWidth(0),
      _imageHeight(0),
      _pixelScaleX(0),
      _pixelScaleY(0),
      _totalSum(0.0),
      _isGriddingFinished(false),
      _weightsAsTaper(false) {}

ImageWeights::ImageWeights(const WeightMode& weightMode, size_t imageWidth,
                           size_t imageHeight, double pixelScaleX,
                           double pixelScaleY, bool weightsAsTaper,
                           double superWeight)
    : _weightMode(weightMode),
      _imageWidth(round(double(imageWidth) / superWeight)),
      _imageHeight(round(double(imageHeight) / superWeight)),
      _pixelScaleX(pixelScaleX),
      _pixelScaleY(pixelScaleY),
      _totalSum(0.0),
      _isGriddingFinished(false),
      _weightsAsTaper(weightsAsTaper) {
  if (_imageWidth % 2 != 0) ++_imageWidth;
  if (_imageHeight % 2 != 0) ++_imageHeight;
  _grid.assign(_imageWidth * _imageHeight / 2, 0.0);
}

void ImageWeights::Serialize(aocommon::SerialOStream& stream) const {
  _weightMode.Serialize(stream);
  stream.UInt64(_imageWidth)
      .UInt64(_imageHeight)
      .Double(_pixelScaleX)
      .Double(_pixelScaleY)
      .Vector(_grid)
      .Double(_totalSum)
      .Bool(_isGriddingFinished)
      .Bool(_weightsAsTaper);
}

void ImageWeights::Unserialize(aocommon::SerialIStream& stream) {
  _weightMode.Unserialize(stream);
  stream.UInt64(_imageWidth)
      .UInt64(_imageHeight)
      .Double(_pixelScaleX)
      .Double(_pixelScaleY)
      .Vector(_grid)
      .Double(_totalSum)
      .Bool(_isGriddingFinished)
      .Bool(_weightsAsTaper);
}

void ImageWeights::Grid(MSProvider& msProvider, const MSSelection& selection) {
  if (_isGriddingFinished)
    throw std::runtime_error("Grid() called after a call to FinishGridding()");
  size_t polarizationCount =
      (msProvider.Polarization() == aocommon::Polarization::Instrumental) ? 4
                                                                          : 1;
  if (_weightMode.RequiresGridding()) {
    SynchronizedMS ms(msProvider.MS());
    const MultiBandData bandData(ms->spectralWindow(), ms->dataDescription());
    MultiBandData selectedBand;
    if (selection.HasChannelRange())
      selectedBand = MultiBandData(bandData, selection.ChannelRangeStart(),
                                   selection.ChannelRangeEnd());
    else
      selectedBand = bandData;
    aocommon::UVector<float> weightBuffer(selectedBand.MaxChannels() *
                                          polarizationCount);

    std::unique_ptr<MSReader> msReader = msProvider.MakeReader();
    while (msReader->CurrentRowAvailable()) {
      double uInM, vInM, wInM;
      size_t dataDescId;
      msReader->ReadMeta(uInM, vInM, wInM, dataDescId);
      msReader->ReadWeights(weightBuffer.data());
      if (_weightsAsTaper) {
        for (float& w : weightBuffer) {
          if (w != 0.0) w = 1.0;
        }
      }
      const BandData& curBand = selectedBand[dataDescId];
      if (vInM < 0.0) {
        uInM = -uInM;
        vInM = -vInM;
      }

      const float* weightIter = weightBuffer.data();
      for (size_t ch = 0; ch != curBand.ChannelCount(); ++ch) {
        double u = uInM / curBand.ChannelWavelength(ch),
               v = vInM / curBand.ChannelWavelength(ch);
        for (size_t p = 0; p != polarizationCount; ++p) {
          Grid(u, v, *weightIter);
          ++weightIter;
        }
      }

      msReader->NextInputRow();
    }
  }
}

void ImageWeights::FinishGridding() {
  if (_isGriddingFinished)
    throw std::runtime_error("FinishGridding() called twice");
  _isGriddingFinished = true;

  switch (_weightMode.Mode()) {
    case WeightMode::BriggsWeighted: {
      double avgW = 0.0;
      for (double val : _grid) avgW += val * val;
      avgW /= _totalSum;
      double numeratorSqrt =
          5.0 * std::exp((-M_LN10) * _weightMode.BriggsRobustness());
      double sSq = numeratorSqrt * numeratorSqrt / avgW;
      for (double& val : _grid) {
        val = 1.0 / (1.0 + val * sSq);
      }
    } break;
    case WeightMode::UniformWeighted: {
      for (double& val : _grid) {
        if (val != 0.0)
          val = 1.0 / val;
        else
          val = 0.0;
      }
    } break;
    case WeightMode::NaturalWeighted: {
      for (double& val : _grid) {
        if (val != 0.0) val = 1.0;
      }
    } break;
  }
}

void ImageWeights::SetMinUVRange(double minUVInLambda) {
  auto iter = _grid.begin();
  const double minSq = minUVInLambda * minUVInLambda;
  int halfWidth = _imageWidth / 2;
  for (size_t y = 0; y != _imageHeight / 2; ++y) {
    for (size_t x = 0; x != _imageWidth; ++x) {
      int xi = int(x) - halfWidth;
      double u = double(xi) / (_imageWidth * _pixelScaleX);
      double v = double(y) / (_imageHeight * _pixelScaleY);
      if (u * u + v * v < minSq) *iter = 0.0;
      ++iter;
    }
  }
}

void ImageWeights::SetMaxUVRange(double maxUVInLambda) {
  auto iter = _grid.begin();
  const double maxSq = maxUVInLambda * maxUVInLambda;
  int halfWidth = _imageWidth / 2;
  for (size_t y = 0; y != _imageHeight / 2; ++y) {
    for (size_t x = 0; x != _imageWidth; ++x) {
      int xi = int(x) - halfWidth;
      double u = double(xi) / (_imageWidth * _pixelScaleX);
      double v = double(y) / (_imageHeight * _pixelScaleY);
      if (u * u + v * v > maxSq) *iter = 0.0;
      ++iter;
    }
  }
}

void ImageWeights::SetTukeyTaper(double transitionSizeInLambda,
                                 double maxUVInLambda) {
  auto iter = _grid.begin();
  const double maxUVSq = maxUVInLambda * maxUVInLambda;
  const double transitionDistSq = (maxUVInLambda - transitionSizeInLambda) *
                                  (maxUVInLambda - transitionSizeInLambda);
  for (size_t y = 0; y != _imageHeight / 2; ++y) {
    for (size_t x = 0; x != _imageWidth; ++x) {
      double u, v;
      xyToUV(x, y, u, v);
      double distSq = u * u + v * v;
      if (distSq > maxUVSq)
        *iter = 0.0;
      else if (distSq > transitionDistSq) {
        *iter *=
            tukeyFrom0ToN(maxUVInLambda - sqrt(distSq), transitionSizeInLambda);
      }
      ++iter;
    }
  }
}

void ImageWeights::SetTukeyInnerTaper(double transitionSizeInLambda,
                                      double minUVInLambda) {
  auto iter = _grid.begin();
  const double minUVSq = minUVInLambda * minUVInLambda;
  const double totalSizeSq = (minUVInLambda + transitionSizeInLambda) *
                             (minUVInLambda + transitionSizeInLambda);
  for (size_t y = 0; y != _imageHeight / 2; ++y) {
    for (size_t x = 0; x != _imageWidth; ++x) {
      double u, v;
      xyToUV(x, y, u, v);
      double distSq = u * u + v * v;
      if (distSq < minUVSq)
        *iter = 0.0;
      else if (distSq < totalSizeSq) {
        *iter *=
            tukeyFrom0ToN(sqrt(distSq) - minUVInLambda, transitionSizeInLambda);
      }
      ++iter;
    }
  }
}

void ImageWeights::SetEdgeTaper(double sizeInLambda) {
  auto iter = _grid.begin();
  double maxU, maxV;
  xyToUV(_imageWidth, _imageHeight / 2, maxU, maxV);
  for (size_t y = 0; y != _imageHeight / 2; ++y) {
    for (size_t x = 0; x != _imageWidth; ++x) {
      double u, v;
      xyToUV(x, y, u, v);
      if (maxU - std::fabs(u) < sizeInLambda ||
          maxV - std::fabs(v) < sizeInLambda)
        *iter = 0.0;
      ++iter;
    }
  }
}

void ImageWeights::SetEdgeTukeyTaper(double transitionSizeInLambda,
                                     double edgeSizeInLambda) {
  auto iter = _grid.begin();
  double maxU, maxV;
  xyToUV(_imageWidth, _imageHeight / 2, maxU, maxV);
  double totalSize = transitionSizeInLambda + edgeSizeInLambda;
  for (size_t y = 0; y != _imageHeight / 2; ++y) {
    for (size_t x = 0; x != _imageWidth; ++x) {
      double u, v;
      xyToUV(x, y, u, v);
      double uDist = maxU - std::fabs(u);
      double vDist = maxV - std::fabs(v);
      if (uDist < edgeSizeInLambda || vDist < edgeSizeInLambda)
        *iter = 0.0;
      else if (uDist < totalSize || vDist < totalSize) {
        double ru = uDist - edgeSizeInLambda;
        double rv = vDist - edgeSizeInLambda;
        if (ru > transitionSizeInLambda) ru = transitionSizeInLambda;
        if (rv > transitionSizeInLambda) rv = transitionSizeInLambda;
        *iter *= tukeyFrom0ToN(ru, transitionSizeInLambda) *
                 tukeyFrom0ToN(rv, transitionSizeInLambda);
      }
      ++iter;
    }
  }
}

void ImageWeights::GetGrid(double* image) const {
  const double* srcPtr = _grid.data();
  for (size_t y = 0; y != _imageHeight / 2; ++y) {
    size_t yUpper = _imageHeight / 2 - 1 - y;
    size_t yLower = _imageHeight / 2 + y;
    double* upperRow = &image[yUpper * _imageWidth];
    double* lowerRow = &image[yLower * _imageWidth];
    for (size_t x = 0; x != _imageWidth; ++x) {
      upperRow[_imageWidth - x - 1] = *srcPtr;
      lowerRow[x] = *srcPtr;
      ++srcPtr;
    }
  }
}

void ImageWeights::Save(const string& filename) const {
  aocommon::UVector<double> image(_imageWidth * _imageHeight);
  GetGrid(&image[0]);
  aocommon::FitsWriter writer;
  writer.SetImageDimensions(_imageWidth, _imageHeight);
  writer.Write(filename, image.data());
}

void ImageWeights::RankFilter(double rankLimit, size_t windowSize) {
  std::vector<double> newGrid(_grid);
  for (size_t y = 0; y != _imageHeight / 2; ++y) {
    for (size_t x = 0; x != _imageWidth; ++x) {
      double w = _grid[y * _imageWidth + x];
      if (w != 0.0) {
        double mean = windowMean(x, y, windowSize);
        if (w > mean * rankLimit)
          newGrid[y * _imageWidth + x] = mean * rankLimit;
      }
    }
  }
  _grid = newGrid;
}

void ImageWeights::SetGaussianTaper(double beamSize) {
  Logger::Debug << "Applying " << Angle::ToNiceString(beamSize)
                << " Gaussian taper...\n";
  double halfPowerUV = 1.0 / (beamSize * 2.0 * M_PI);
  const long double sigmaToHP = 2.0L * sqrtl(2.0L * logl(2.0L));
  double minusTwoSigmaSq = halfPowerUV * sigmaToHP;
  Logger::Debug << "UV taper: " << minusTwoSigmaSq << '\n';
  minusTwoSigmaSq *= -2.0 * minusTwoSigmaSq;
  for (size_t y = 0; y != _imageHeight / 2; ++y) {
    for (size_t x = 0; x != _imageWidth; ++x) {
      double val = _grid[y * _imageWidth + x];
      if (val != 0.0) {
        double u, v;
        xyToUV(x, y, u, v);
        double gaus = exp((u * u + v * v) / minusTwoSigmaSq);
        _grid[y * _imageWidth + x] = val * gaus;
      }
    }
  }
}

double ImageWeights::windowMean(size_t x, size_t y, size_t windowSize) {
  size_t d = windowSize / 2;
  size_t x1, y1, x2, y2;
  if (x <= d)
    x1 = 0;
  else
    x1 = x - d;

  if (y <= d)
    y1 = 0;
  else
    y1 = y - d;

  if (x + d >= _imageWidth)
    x2 = _imageWidth;
  else
    x2 = x + d;

  if (y + d >= _imageHeight / 2)
    y2 = _imageHeight / 2;
  else
    y2 = y + d;

  size_t windowCount = 0;
  double windowSum = 0.0;
  for (size_t yi = y1; yi < y2; ++yi) {
    for (size_t xi = x1; xi < x2; ++xi) {
      double w = _grid[yi * _imageWidth + xi];
      if (w != 0.0) {
        ++windowCount;
        windowSum += w;
      }
    }
  }
  return windowSum / double(windowCount);
}
