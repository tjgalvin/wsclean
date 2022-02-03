#include "iuwtdeconvolutionalgorithm.h"

#include "imageanalysis.h"

#include "../system/fftwmanager.h"

#include "../math/fftconvolver.h"
#include "../math/gaussianfitter.h"

#include "../deconvolution/imageset.h"

#include <aocommon/image.h>
#include <aocommon/system.h>

#include <boost/numeric/conversion/bounds.hpp>

#include <algorithm>
#include <iostream>

using aocommon::Image;

IUWTDeconvolutionAlgorithm::IUWTDeconvolutionAlgorithm(
    FFTWManager& fftwManager, size_t width, size_t height, float gain,
    float mGain, float cleanBorder, bool allowNegativeComponents,
    const bool* mask, float absoluteThreshold, float thresholdSigmaLevel,
    float tolerance, bool useSNRTest)
    : _fftwManager(fftwManager),
      _width(width),
      _height(height),
      _gain(gain),
      _mGain(mGain),
      _cleanBorder(cleanBorder),
      _mask(mask),
      _absoluteThreshold(absoluteThreshold),
      _thresholdSigmaLevel(thresholdSigmaLevel),
      _tolerance(tolerance),
      _allowNegativeComponents(allowNegativeComponents),
      _useSNRTest(useSNRTest) {}

void IUWTDeconvolutionAlgorithm::measureRMSPerScale(
    const float* image, const float* convolvedImage, float* scratch,
    size_t endScale, std::vector<ScaleResponse>& psfResponse) {
  IUWTDecomposition imageIUWT(endScale, _width, _height);
  imageIUWT.Decompose(*_staticFor, image, scratch, false);

  _psfMaj = 2.0;
  _psfMin = 2.0;
  _psfPA = 0.0;
  double fl = 0.0;
  GaussianFitter fitter;
  fitter.Fit2DGaussianCentred(image, _width, _height, 2.0, _psfMaj, _psfMin,
                              _psfPA);
  _psfVolume = (M_PI / 4.0) * _psfMaj * _psfMin / M_LOG2E;

  double v = 1.0, x = _width / 2, y = _height / 2;
  double bMaj = _psfMaj, bMin = _psfMin, bPA = _psfPA;
  fitter.Fit2DGaussianFull(image, _width, _height, v, x, y, bMaj, bMin, bPA,
                           &fl);

  psfResponse.resize(endScale);
  for (size_t scale = 0; scale != endScale; ++scale) {
    psfResponse[scale].rms = rms(imageIUWT[scale].Coefficients());
    psfResponse[scale].peakResponse =
        centralPeak(imageIUWT[scale].Coefficients());
    bMaj = 2.0;
    bMin = 2.0;
    bPA = 0.0;
    v = 1.0;
    x = _width / 2;
    y = _height / 2;
    fitter.Fit2DGaussianFull(imageIUWT[scale].Coefficients().Data(), _width,
                             _height, v, x, y, bMaj, bMin, bPA, &fl);
    psfResponse[scale].bMaj = bMaj;
    psfResponse[scale].bMin = bMin;
    psfResponse[scale].bPA = bPA;
  }

  imageIUWT.Decompose(*_staticFor, imageIUWT[1].Coefficients().Data(), scratch,
                      false);
  for (size_t scale = 0; scale != endScale; ++scale) {
    psfResponse[scale].peakResponseToNextScale =
        centralPeak(imageIUWT[scale].Coefficients());
  }

  imageIUWT.Decompose(*_staticFor, convolvedImage, scratch, false);

  for (size_t scale = 0; scale != endScale; ++scale) {
    psfResponse[scale].convolvedPeakResponse =
        centralPeak(imageIUWT[scale].Coefficients());
  }

  aocommon::UVector<float> thresholds(imageIUWT.NScales());
  for (size_t i = 0; i != imageIUWT.NScales(); ++i) {
    thresholds[i] = psfResponse[0].convolvedPeakResponse * _tolerance;
  }
  IUWTMask mask(imageIUWT.NScales(), _width, _height);
  ImageAnalysis::Component component(_width / 2, _height / 2, 0);
  size_t areaSize;
  ImageAnalysis::Floodfill(imageIUWT, mask, thresholds, 0,
                           std::min<size_t>(endScale, 2), component, 0.0,
                           areaSize);
  aocommon::UVector<bool> markedMask0(mask[0].size(), false);
  ImageAnalysis::Component2D c2D(_width / 2, _height / 2);
  float threshold = psfResponse[0].convolvedPeakResponse * _tolerance;
  ImageAnalysis::FloodFill2D(imageIUWT[0].Coefficients().Data(),
                             markedMask0.data(), threshold, c2D, _width,
                             _height, psfResponse[0].convolvedArea);
}

float IUWTDeconvolutionAlgorithm::mad(const float* dest) {
  Image v(_width, _height);
  for (size_t i = 0; i != _width * _height; ++i) v[i] = std::fabs(dest[i]);
  size_t mid = (_width * _height) / 2;
  std::nth_element(v.begin(), v.begin() + mid, v.end());
  return v[mid] / 0.674559;
}

float IUWTDeconvolutionAlgorithm::getMaxAbsWithoutMask(const Image& data,
                                                       size_t& x, size_t& y,
                                                       size_t width) {
  size_t height = data.Size() / width;
  size_t xBorder = _cleanBorder * width;
  size_t yBorder = _cleanBorder * height;
  size_t minX = xBorder, maxX = width - xBorder;
  size_t minY = yBorder, maxY = height - yBorder;
  x = width;
  y = height;

  float maxVal = boost::numeric::bounds<float>::lowest();
  for (size_t yi = minY; yi != maxY; ++yi) {
    const float* dataPtr = data.Data() + yi * width;
    for (size_t xi = minX; xi != maxX; ++xi) {
      float val =
          _allowNegativeComponents ? std::fabs(dataPtr[xi]) : dataPtr[xi];
      if (val > maxVal) {
        maxVal = val;
        x = xi;
        y = yi;
      }
    }
  }
  return maxVal;
}

float IUWTDeconvolutionAlgorithm::getMaxAbsWithMask(const Image& data,
                                                    size_t& x, size_t& y,
                                                    size_t width) {
  size_t height = data.Size() / width;
  size_t xBorder = _cleanBorder * width;
  size_t yBorder = _cleanBorder * height;
  size_t minX = xBorder, maxX = width - xBorder;
  size_t minY = yBorder, maxY = height - yBorder;
  x = width;
  y = height;

  float maxVal = boost::numeric::bounds<float>::lowest();
  for (size_t yi = minY; yi != maxY; ++yi) {
    const float* dataPtr = data.Data() + yi * width;
    const bool* maskPtr = _mask + yi * width;
    for (size_t xi = minX; xi != maxX; ++xi) {
      if (maskPtr[xi]) {
        float val =
            _allowNegativeComponents ? std::fabs(dataPtr[xi]) : dataPtr[xi];
        if (val > maxVal) {
          maxVal = val;
          x = xi;
          y = yi;
        }
      }
    }
  }
  return maxVal;
}

float IUWTDeconvolutionAlgorithm::dotProduct(const Image& lhs,
                                             const Image& rhs) {
  float sum = 0.0;
  for (size_t i = 0; i != lhs.Size(); ++i) sum += lhs[i] * rhs[i];
  return sum;
}

void IUWTDeconvolutionAlgorithm::factorAdd(float* dest, const float* rhs,
                                           float factor, size_t width,
                                           size_t height) {
  for (size_t i = 0; i != width * height; ++i) dest[i] += rhs[i] * factor;
}

void IUWTDeconvolutionAlgorithm::factorAdd(Image& dest, const Image& rhs,
                                           float factor) {
  for (size_t i = 0; i != dest.Size(); ++i) dest[i] += rhs[i] * factor;
}

void IUWTDeconvolutionAlgorithm::Subtract(float* dest, const Image& rhs) {
  for (size_t i = 0; i != rhs.Size(); ++i) dest[i] -= rhs[i];
}

void IUWTDeconvolutionAlgorithm::boundingBox(size_t& x1, size_t& y1, size_t& x2,
                                             size_t& y2, const Image& image,
                                             size_t width, size_t height) {
  float mP = *std::max_element(image.begin(), image.end());
  float mN = *std::min_element(image.begin(), image.end());
  float m = std::max(mP, -mN);
  x1 = width;
  x2 = 0;
  y1 = height;
  y2 = 0;
  for (size_t y = 0; y != height; ++y) {
    const float* ptr = image.Data() + y * width;
    for (size_t x = 0; x != x1; ++x) {
      if (std::fabs(ptr[x]) > m * 0.01) {
        x1 = x;
        break;
      }
    }
    for (size_t x = width - 1; x != x2; --x) {
      if (std::fabs(ptr[x]) > m * 0.01) {
        x2 = x;
        break;
      }
    }
  }
  x2++;
  for (size_t y = 0; y != height; ++y) {
    const float* ptr = image.Data() + y * width;
    for (size_t x = 0; x != width; ++x) {
      if (std::fabs(ptr[x]) > m * 0.01) {
        if (y1 > y) y1 = y;
        if (y2 < y) y2 = y + 1;
      }
    }
  }
}

void IUWTDeconvolutionAlgorithm::adjustBox(size_t& x1, size_t& y1, size_t& x2,
                                           size_t& y2, size_t width,
                                           size_t height, int endScale) {
  const int minBoxSize = std::max<int>(
      128, IUWTDecomposition::MinImageDimension(endScale) * 3 / 2);

  int boxWidth = x2 - x1;
  int boxHeight = y2 - y1;
  int newX1 = x1 - 0.5 * boxWidth, newX2 = x2 + 0.5 * boxWidth,
      newY1 = y1 - 0.5 * boxHeight, newY2 = y2 + 0.5 * boxHeight;

  if (newX2 - newX1 < minBoxSize) {
    int mid = 0.5 * (int(x1) + int(x2));
    newX1 = mid - minBoxSize / 2;
    newX2 = mid + minBoxSize / 2;
  }
  if (newY2 - newY1 < minBoxSize) {
    int mid = 0.5 * (int(y1) + int(y2));
    newY1 = mid - minBoxSize / 2;
    newY2 = mid + minBoxSize / 2;
  }
  if (newX1 >= 0)
    x1 = newX1;
  else
    x1 = 0;
  if (newX2 < int(width))
    x2 = newX2;
  else
    x2 = width;
  if (newY1 >= 0)
    y1 = newY1;
  else
    y1 = 0;
  if (newY2 < int(height))
    y2 = newY2;
  else
    y2 = height;
  while ((x2 - x1) % 8 != 0) x2--;
  while ((y2 - y1) % 8 != 0) y2--;
}

void IUWTDeconvolutionAlgorithm::trim(Image& dest, const float* source,
                                      size_t oldWidth, size_t /*oldHeight*/,
                                      size_t x1, size_t y1, size_t x2,
                                      size_t y2) {
  // We do this so that dest and source can be the same image.
  Image scratch((x2 - x1), (y2 - y1));
  for (size_t y = y1; y != y2; ++y) {
    const float* oldPtr = &source[y * oldWidth];
    float* newPtr = &scratch[(y - y1) * (x2 - x1)];
    for (size_t x = x1; x != x2; ++x) {
      newPtr[x - x1] = oldPtr[x];
    }
  }
  dest = scratch;
}

void IUWTDeconvolutionAlgorithm::untrim(Image& image, size_t width,
                                        size_t height, size_t x1, size_t y1,
                                        size_t x2, size_t y2) {
  Image scratch(width, height, 0.0);
  std::copy_n(image.Data(), image.Width() * image.Height(), scratch.Data());
  image = scratch;
  size_t y = y2;
  while (y != y1) {
    --y;
    float* newPtr = &image[y * width];
    float* oldPtr = &image[(y - y1) * (x2 - x1)];
    size_t x = x2;
    while (x != x1) {
      --x;
      newPtr[x] = oldPtr[x - x1];
    }
  }
  for (size_t y = 0; y != y1; ++y) {
    float* ptr = &image[y * width];
    for (size_t x = 0; x != width; ++x) ptr[x] = 0;
  }
  for (size_t y = y1; y != y2; ++y) {
    float* ptr = &image[y * width];
    for (size_t x = 0; x != x1; ++x) ptr[x] = 0.0;
    for (size_t x = x2; x != width; ++x) ptr[x] = 0.0;
  }
  for (size_t y = y2; y != height; ++y) {
    float* ptr = &image[y * width];
    for (size_t x = 0; x != width; ++x) ptr[x] = 0;
  }
}

float IUWTDeconvolutionAlgorithm::sum(const Image& img) const {
  float s = 0.0;
  for (size_t i = 0; i != img.Size(); ++i) s += img[i];
  return s;
}

float IUWTDeconvolutionAlgorithm::snr(const IUWTDecomposition& noisyImg,
                                      const IUWTDecomposition& model) const {
  float mSum = 0.0, nSum = 0.0;
  for (size_t scale = 0; scale != noisyImg.NScales(); ++scale) {
    const Image& n = noisyImg[scale].Coefficients();
    const Image& m = model[scale].Coefficients();
    for (size_t i = 0; i != n.Size(); ++i) {
      mSum += m[i] * m[i];
      float d = m[i] - n[i];
      nSum += d * d;
    }
  }
  return mSum / nSum;
}

float IUWTDeconvolutionAlgorithm::rmsDiff(const Image& a, const Image& b) {
  float sum = 0.0;
  for (size_t i = 0; i != a.Size(); ++i) {
    float d = a[i] - b[i];
    sum += d * d;
  }
  return std::sqrt(sum / a.Size());
}

float IUWTDeconvolutionAlgorithm::rms(const Image& image) {
  float sum = 0.0;
  for (size_t i = 0; i != image.Size(); ++i) {
    float v = image[i];
    sum += v * v;
  }
  return std::sqrt(sum / image.Size());
}

bool IUWTDeconvolutionAlgorithm::runConjugateGradient(
    IUWTDecomposition& iuwt, const IUWTMask& mask, Image& maskedDirty,
    Image& structureModel, Image& scratch, const Image& psfKernel, size_t width,
    size_t height) {
  Image gradient = maskedDirty;
  float modelSNR = 0.0;

  IUWTDecomposition initialDirtyIUWT(iuwt);

  for (size_t minorIter = 0; minorIter != 20; ++minorIter) {
    // scratch = gradient (x) psf
    scratch = gradient;
    FFTConvolver::ConvolveSameSize(_fftwManager, scratch.Data(),
                                   psfKernel.Data(), width, height,
                                   _staticFor->NThreads());

    // calc: IUWT gradient (x) psf
    iuwt.Decompose(*_staticFor, scratch.Data(), scratch.Data(), false);

    // calc: mask IUWT gradient (x) psf
    iuwt.ApplyMask(mask);

    // scratch = IUWT^-1 mask IUWT gradient (x) psf
    iuwt.Recompose(scratch, false);

    // stepsize = <residual, residual> / <gradient, scratch>
    float gradientDotScratch = dotProduct(gradient, scratch);
    if (gradientDotScratch == 0.0) return false;
    float stepSize = dotProduct(maskedDirty, maskedDirty) / gradientDotScratch;

    // model_i+1 = model_i + stepsize * gradient
    factorAdd(structureModel.Data(), gradient.Data(), stepSize, width, height);

    // For Dabbech's approach (see below) :
    //  Image scratch2 = maskedDirty;

    float gradStepDen = dotProduct(maskedDirty, maskedDirty);
    if (gradStepDen == 0.0) return false;
    // residual_i+1 = residual_i - stepsize * scratch
    factorAdd(maskedDirty.Data(), scratch.Data(), -stepSize, width, height);

    // PyMORESANE uses this:
    // gradstep = <residual_i+1, residual_i+1> / <residual_i, residual_i>
    // float gradStep = dotProduct(maskedDirty, maskedDirty) / gradStepDen;
    // But in MORESANE's paper A. Dabbech says this:
    // gradstep = <residual_i+1 - residual_i, residual_i+1> / <residual_i,
    // residual_i> scratch = maskedDirty; subtract(scratch, scratch2); float
    // gradStep = dotProduct(scratch, maskedDirty) / gradStepDen;
    float gradStep = dotProduct(maskedDirty, maskedDirty) / gradStepDen;

    // gradient_i+1 = residual_i+1 + gradstep * gradient_i
    scratch = gradient;
    gradient = maskedDirty;
    factorAdd(gradient.Data(), scratch.Data(), gradStep, width, height);

    // scratch = mask IUWT PSF (x) model
    scratch = structureModel;
    FFTConvolver::ConvolveSameSize(_fftwManager, scratch.Data(),
                                   psfKernel.Data(), width, height,
                                   _staticFor->NThreads());
    iuwt.Decompose(*_staticFor, scratch.Data(), scratch.Data(), false);
    iuwt.ApplyMask(mask);

    float previousSNR = modelSNR;
    modelSNR = snr(iuwt, initialDirtyIUWT);
    if (modelSNR > 100 && minorIter > 2) {
      std::cout << "Converged after " << minorIter << " iterations.\n";
      return true;
    } else if (modelSNR < previousSNR && minorIter > 5) {
      if (modelSNR > 3) {
        std::cout << "SNR decreased after " << minorIter
                  << " iterations (SNR=" << modelSNR << ").\n";
        return true;
      }
    }
  }
  if (modelSNR <= 3.0) {
    std::cout << "Failed to converge (SNR=" << modelSNR << ").\n";
    structureModel = Image(width, height, 0.0);
    return false;
  }
  return true;
}

struct PointSource {
  float x, y, flux;
  bool operator<(const PointSource& rhs) const { return flux < rhs.flux; }
};

void IUWTDeconvolutionAlgorithm::constrainedPSFConvolve(float* image,
                                                        const float* psf,
                                                        size_t width,
                                                        size_t height) {
  Image smallerPsf(width, height, 0.0);
  Image kernel(width, height);
  size_t s = std::round(std::sqrt(_psfResponse[0].convolvedArea * 25.0));
  size_t smallWidth = std::min(s, width);
  size_t smallHeight = std::min(s, height);
  std::cout << "Constrained PSF=" << smallWidth << " x " << smallHeight << '\n';
  size_t xMin = width / 2 - smallWidth / 2, xMax = width / 2 + smallWidth / 2;
  size_t yMin = height / 2 - smallHeight / 2,
         yMax = height / 2 + smallHeight / 2;
  for (size_t y = yMin; y != yMax; ++y) {
    for (size_t x = xMin; x != xMax; ++x) {
      smallerPsf[y * width + x] = psf[y * width + x];
    }
  }
  FFTConvolver::PrepareKernel(kernel.Data(), smallerPsf.Data(), width, height,
                              _staticFor->NThreads());
  FFTConvolver::ConvolveSameSize(_fftwManager, image, kernel.Data(), width,
                                 height, _staticFor->NThreads());
}

bool IUWTDeconvolutionAlgorithm::findAndDeconvolveStructure(
    IUWTDecomposition& iuwt, Image& dirty, const Image& psf,
    const Image& psfKernel, Image& scratch, ImageSet& structureModel,
    size_t curEndScale, size_t curMinScale,
    std::vector<IUWTDeconvolutionAlgorithm::ValComponent>& maxComponents) {
  iuwt.Decompose(*_staticFor, dirty.Data(), scratch.Data(), false);
  aocommon::UVector<float> thresholds(curEndScale);
  _rmses.resize(curEndScale);
  for (size_t scale = 0; scale != curEndScale; ++scale) {
    float r = mad(iuwt[scale].Coefficients().Data());
    _rmses[scale] = r;
    thresholds[scale] = r * (_thresholdSigmaLevel * 4.0 / 5.0);
  }

  scratch = dirty;
  maxComponents.resize(curEndScale);
  for (size_t scale = 0; scale != curEndScale; ++scale) {
    size_t x, y;
    float maxAbsCoef = getMaxAbs(iuwt[scale].Coefficients(), x, y, _width);
    maxComponents[scale].x = x;
    maxComponents[scale].y = y;
    maxComponents[scale].scale = scale;
    maxComponents[scale].val = maxAbsCoef;
  }

  float maxVal = -1.0;
  size_t maxX = 0, maxY = 0;
  int maxValScale = -1;
  for (size_t scale = 0; scale != curEndScale; ++scale) {
    // Considerations for this section:
    // - Scale 0 should be chosen if the input corresponds to the PSF.
    //   Therefore, a peak on scale 1 should be at least:
    //   (PSF peak on scale 1) * (peak on scale 0) / (PSF (x) scale 1 peak
    //   response) Such that anything smaller than scale 1 will be considered
    //   scale 0.

    const ValComponent& val = maxComponents[scale];
    float absCoef = val.val / _psfResponse[scale].rms;
    // std::cout << scale << ">=" << curMinScale << " && " << absCoef << " > "
    // << maxVal << " && " << val.val << " > " << _rmses[scale]*_thresholdLevel
    // << "\n";
    if (scale >= curMinScale && absCoef > maxVal &&
        val.val > _rmses[scale] * _thresholdSigmaLevel &&
        val.val > _rmses[scale] / _rmses[0] * _absoluteThreshold) {
      maxX = val.x;
      maxY = val.y;
      maxValScale = scale;
      if (scale == 0) {
        float lowestRMS = std::min(_psfResponse[0].rms, _psfResponse[1].rms);
        maxVal = val.val / lowestRMS * _psfResponse[1].peakResponse /
                 _psfResponse[0].peakResponseToNextScale;
      } else
        maxVal = absCoef;
    }
  }
  if (maxValScale == -1) {
    std::cout << "No significant pixel found.\n";
    return false;
  }

  maxVal = iuwt[maxValScale][maxX + maxY * _width];
  std::cout << "Most significant pixel: " << maxX << ',' << maxY << "="
            << maxVal << " (" << maxVal / _rmses[maxValScale] << "σ) on scale "
            << maxValScale << '\n';

  if (std::fabs(maxVal) < thresholds[maxValScale]) {
    std::cout << "Most significant pixel is in the noise, stopping.\n";
    return false;
  }

  float scaleMaxAbsVal = std::fabs(maxVal);
  for (size_t scale = 0; scale != curEndScale; ++scale) {
    if (thresholds[scale] < _tolerance * scaleMaxAbsVal) {
      thresholds[scale] = _tolerance * scaleMaxAbsVal;
    }
    if (maxVal < 0.0) thresholds[scale] = -thresholds[scale];
  }

  ImageAnalysis::Component maxComp(maxX, maxY, maxValScale);
  return fillAndDeconvolveStructure(iuwt, dirty, structureModel, scratch, psf,
                                    psfKernel, curEndScale, curMinScale, _width,
                                    _height, thresholds, maxComp, true, _mask);
}

bool IUWTDeconvolutionAlgorithm::fillAndDeconvolveStructure(
    IUWTDecomposition& iuwt, Image& dirty, ImageSet& structureModelFull,
    Image& scratch, const Image& psf, const Image& psfKernel,
    size_t curEndScale, size_t curMinScale, size_t width, size_t height,
    const aocommon::UVector<float>& thresholds,
    const ImageAnalysis::Component& maxComp, bool allowTrimming,
    const bool* priorMask) {
  IUWTMask mask(curEndScale, width, height);
  size_t areaSize;
  ImageAnalysis::SelectStructures(iuwt, mask, thresholds, curMinScale,
                                  curEndScale, _cleanBorder, priorMask,
                                  areaSize);
  std::cout << "Flood-filled area contains " << areaSize
            << " significant components.\n";

  iuwt.ApplyMask(mask);
  iuwt.Recompose(scratch, false);

  // Find bounding box
  size_t x1, y1, x2, y2;
  boundingBox(x1, y1, x2, y2, scratch, width, height);
  adjustBox(x1, y1, x2, y2, width, height, maxComp.scale + 1);
  if (allowTrimming && ((x2 - x1) < width || (y2 - y1) < height)) {
    _curBoxXStart = x1;
    _curBoxXEnd = x2;
    _curBoxYStart = y1;
    _curBoxYEnd = y2;
    std::cout << "Bounding box: (" << x1 << ',' << y1 << ")-(" << x2 << ','
              << y2 << ")\n";
    size_t newWidth = x2 - x1, newHeight = y2 - y1;
    trim(dirty, dirty, width, height, x1, y1, x2, y2);
    Image smallPSF;

    trimPsf(smallPSF, psf.Data(), width, height, newWidth, newHeight);

    Image smallPSFKernel(smallPSF.Width(), smallPSF.Height());
    FFTConvolver::PrepareKernel(smallPSFKernel.Data(), smallPSF.Data(),
                                newWidth, newHeight, _staticFor->NThreads());

    scratch = Image(dirty.Width(), dirty.Height());

    int maxScale =
        std::max(IUWTDecomposition::EndScale(std::min(x2 - x1, y2 - y1)),
                 maxComp.scale + 1);
    if (maxScale < int(curEndScale)) {
      std::cout << "Bounding box too small for largest scale of " << curEndScale
                << " -- ignoring scales>=" << maxScale
                << " in this iteration.\n";
      curEndScale = maxScale;
    }
    std::unique_ptr<IUWTDecomposition> trimmedIUWT(
        iuwt.CreateTrimmed(curEndScale, x1, y1, x2, y2));

    std::unique_ptr<ImageSet> trimmedStructureModel(
        structureModelFull.Trim(x1, y1, x2, y2, width));

    aocommon::UVector<bool> trimmedPriorMask;
    bool* trimmedPriorMaskPtr;
    if (priorMask == nullptr)
      trimmedPriorMaskPtr = nullptr;
    else {
      trimmedPriorMask.resize(newWidth * newHeight);
      trimmedPriorMaskPtr = trimmedPriorMask.data();
      Image::TrimBox(trimmedPriorMaskPtr, x1, y1, newWidth, newHeight,
                     priorMask, width, height);
    }

    ImageAnalysis::Component newMaxComp(maxComp.x - x1, maxComp.y - y1,
                                        maxComp.scale);
    bool result = fillAndDeconvolveStructure(
        *trimmedIUWT, dirty, *trimmedStructureModel, scratch, smallPSF,
        smallPSFKernel, curEndScale, curMinScale, x2 - x1, y2 - y1, thresholds,
        newMaxComp, false, trimmedPriorMaskPtr);
    for (size_t i = 0; i != structureModelFull.size(); ++i) {
      std::copy_n((*trimmedStructureModel)[i], (y2 - y1) * (x2 - x1),
                  scratch.Data());
      untrim(scratch, width, height, x1, y1, x2, y2);
      std::copy_n(scratch.Data(), width * height, structureModelFull[i]);
    }

    dirty = Image(scratch.Width(), scratch.Height());
    _curBoxXStart = 0;
    _curBoxXEnd = width;
    _curBoxYStart = 0;
    _curBoxYEnd = height;
    return result;
  } else {
    if (curEndScale <= 3) {
      // TODO: remove?
      // bool pointSourcesWereFound = extractPointSources(iuwt, mask,
      // dirty.Data(), structureModel.Data()); if(pointSourcesWereFound)
      // return
      // true;
    }

    // get undeconvolved dirty
    iuwt.Decompose(*_staticFor, dirty.Data(), scratch.Data(), false);

    iuwt.ApplyMask(mask);
    iuwt.Recompose(scratch, false);

    Image maskedDirty = scratch;

    Image structureModel(width, height, 0.0);
    bool success = runConjugateGradient(iuwt, mask, maskedDirty, structureModel,
                                        scratch, psfKernel, width, height);
    if (!success) return false;

    float rmsBefore = rms(dirty);
    scratch = structureModel;
    FFTConvolver::ConvolveSameSize(_fftwManager, scratch.Data(),
                                   psfKernel.Data(), width, height,
                                   _staticFor->NThreads());
    maskedDirty = dirty;  // we use maskedDirty as temporary
    factorAdd(maskedDirty.Data(), scratch.Data(), -_gain, width, height);
    float rmsAfter = rms(maskedDirty);
    if (rmsAfter > rmsBefore) {
      std::cout << "RMS got worse: " << rmsBefore << " -> " << rmsAfter << '\n';
      return false;
    }

    // TODO when only one image is available, this is not necessary
    performSubImageFitAll(iuwt, mask, structureModel, scratch, maskedDirty,
                          maxComp, structureModelFull, psf.Data(), dirty);

    return true;
  }
}

void IUWTDeconvolutionAlgorithm::performSubImageFitAll(
    IUWTDecomposition& iuwt, const IUWTMask& mask, const Image& structureModel,
    Image& scratchA, Image& scratchB, const ImageAnalysis::Component& maxComp,
    ImageSet& fittedModel, const float* psf, const Image& dirty) {
  size_t width = iuwt.Width(), height = iuwt.Height();

  if (_dirtySet->size() == 1) {
    // With only one image, we don't have to refit
    Image img(width, height);
    std::copy_n(structureModel.Data(), width * height, img.begin());
    fittedModel.SetImage(0, std::move(img));
  } else {
    std::cout << "Fitting structure in images: " << std::flush;
    aocommon::UVector<float> correctionFactors;
    scratchA = dirty;
    performSubImageFitSingle(iuwt, mask, structureModel, scratchB, maxComp, psf,
                             scratchA, nullptr, correctionFactors);

    fittedModel = 0.0;

    for (size_t imgIndex = 0; imgIndex != _dirtySet->size(); ++imgIndex) {
      std::cout << '.' << std::flush;
      const float* subPsf = _psfs[_dirtySet->PSFIndex(imgIndex)];

      trim(scratchA, (*_dirtySet)[imgIndex], _width, _height, _curBoxXStart,
           _curBoxYStart, _curBoxXEnd, _curBoxYEnd);

      Image smallSubPsf;
      const float* subPsfData;
      if (_width != width || _height != height) {
        trimPsf(smallSubPsf, subPsf, _width, _height, width, height);
        subPsfData = smallSubPsf.Data();
      } else {
        subPsfData = subPsf;
      }

      performSubImageFitSingle(iuwt, mask, structureModel, scratchB, maxComp,
                               subPsfData, scratchA, fittedModel[imgIndex],
                               correctionFactors);
    }
    std::cout << '\n';
  }
}

void IUWTDeconvolutionAlgorithm::performSubImageFitSingle(
    IUWTDecomposition& iuwt, const IUWTMask& mask, const Image& structureModel,
    Image& scratchB, const ImageAnalysis::Component& maxComp, const float* psf,
    Image& subDirty, float* fittedSubModel,
    aocommon::UVector<float>& correctionFactors) {
  size_t width = iuwt.Width(), height = iuwt.Height();

  Image psfKernel(width, height);
  FFTConvolver::PrepareKernel(psfKernel.Data(), psf, width, height,
                              _staticFor->NThreads());

  Image& maskedDirty = scratchB;

  iuwt.Decompose(*_staticFor, subDirty.Data(), subDirty.Data(), false);
  iuwt.ApplyMask(mask);
  iuwt.Recompose(maskedDirty, false);
  aocommon::UVector<bool> mask2d(structureModel.Size(), false);
  float peakLevel = std::fabs(structureModel[maxComp.y * width + maxComp.x]);
  size_t componentIndex = 0;
  for (size_t y = 0; y != height; ++y) {
    bool* maskRow = &mask2d[y * width];
    const float* modelRow = &structureModel[y * width];
    for (size_t x = 0; x != width; ++x) {
      if (!maskRow[x] && std::fabs(modelRow[x]) > peakLevel * 1e-4) {
        std::vector<ImageAnalysis::Component2D> area;
        ImageAnalysis::Component2D comp(x, y);
        ImageAnalysis::FloodFill2D(structureModel.Data(), mask2d.data(),
                                   peakLevel * 1e-4, comp, width, height, area);
        // Find bounding box and copy active pixels to subDirty
        subDirty = Image(width, height, 0.0);
        size_t boxX1 = width, boxX2 = 0, boxY1 = height, boxY2 = 0;
        for (std::vector<ImageAnalysis::Component2D>::const_iterator a =
                 area.begin();
             a != area.end(); ++a) {
          size_t index = a->x + a->y * width;
          boxX1 = std::min(a->x, boxX1);
          boxX2 = std::max(a->x, boxX2);
          boxY1 = std::min(a->y, boxY1);
          boxY2 = std::max(a->y, boxY2);
          subDirty[index] = structureModel[index];
        }
        adjustBox(boxX1, boxY1, boxX2, boxY2, width, height, iuwt.NScales());

        float factor = performSubImageComponentFitBoxed(
            iuwt, mask, area, subDirty, maskedDirty, psf, psfKernel, boxX1,
            boxY1, boxX2, boxY2);

        // if no fittedSubModel was given, we just need to store the factors.
        // Otherwise, scale the deconvolved model and add it to the contaminated
        // model.
        if (fittedSubModel != nullptr) {
          const float integratedFactor = correctionFactors[componentIndex];
          if (std::isfinite(factor) && std::isfinite(integratedFactor) &&
              integratedFactor != 0.0) {
            for (std::vector<ImageAnalysis::Component2D>::const_iterator a =
                     area.begin();
                 a != area.end(); ++a) {
              size_t index = a->x + a->y * width;
              fittedSubModel[index] +=
                  structureModel[index] * factor / integratedFactor;
            }
          }
          ++componentIndex;
        } else {
          correctionFactors.push_back(factor);
        }
      }
    }
  }
}

float IUWTDeconvolutionAlgorithm::performSubImageComponentFitBoxed(
    IUWTDecomposition& iuwt, const IUWTMask& mask,
    const std::vector<ImageAnalysis::Component2D>& area, Image& model,
    Image& maskedDirty, const float* psf, const Image& psfKernel, size_t x1,
    size_t y1, size_t x2, size_t y2) {
  const size_t width = iuwt.Width(), height = iuwt.Height();
  if (x1 > 0 || y1 > 0 || x2 < width || y2 < height) {
    size_t newWidth = x2 - x1, newHeight = y2 - y1;
    IUWTDecomposition smallIUWTW(iuwt.NScales(), newWidth, newHeight);
    std::unique_ptr<IUWTMask> smallMask(mask.CreateTrimmed(x1, y1, x2, y2));
    Image smallModel;
    trim(smallModel, model, width, height, x1, y1, x2, y2);

    Image smallPsf;
    trimPsf(smallPsf, psf, width, height, newWidth, newHeight);
    Image smallPsfKernel(smallPsf.Width(), smallPsf.Height());
    FFTConvolver::PrepareKernel(smallPsfKernel.Data(), smallPsf.Data(),
                                newWidth, newHeight, _staticFor->NThreads());

    Image smallMaskedDirty;
    trim(smallMaskedDirty, maskedDirty, width, height, x1, y1, x2, y2);

    float factor =
        performSubImageComponentFit(smallIUWTW, *smallMask, area, smallModel,
                                    smallMaskedDirty, smallPsfKernel, x1, y1);
    return factor;
  } else {
    return performSubImageComponentFit(iuwt, mask, area, model, maskedDirty,
                                       psfKernel, 0, 0);
  }
}

float IUWTDeconvolutionAlgorithm::performSubImageComponentFit(
    IUWTDecomposition& iuwt, const IUWTMask& mask,
    const std::vector<ImageAnalysis::Component2D>& area, Image& model,
    Image& maskedDirty, const Image& psfKernel, size_t xOffset,
    size_t yOffset) {
  const size_t width = iuwt.Width(), height = iuwt.Height();
  // Calculate IUWT^-1 mask IUWT model (x) PSF
  FFTConvolver::ConvolveSameSize(_fftwManager, model.Data(), psfKernel.Data(),
                                 width, height, _staticFor->NThreads());
  iuwt.Decompose(*_staticFor, model.Data(), model.Data(), false);
  iuwt.ApplyMask(mask);
  iuwt.Recompose(model, false);

  float modelSum = 0.0, dirtySum = 0.0;
  for (std::vector<ImageAnalysis::Component2D>::const_iterator a = area.begin();
       a != area.end(); ++a) {
    size_t index = (a->x - xOffset) + (a->y - yOffset) * width;
    modelSum += model[index];
    dirtySum += maskedDirty[index];
  }
  // std::cout << "factor=" << dirtySum << " / " << modelSum << " = " <<
  // dirtySum/modelSum << '\n';
  if (modelSum == 0.0 || !std::isfinite(dirtySum) || !std::isfinite(modelSum))
    return 0.0;
  else
    return dirtySum / modelSum;
}

float IUWTDeconvolutionAlgorithm::PerformMajorIteration(
    size_t& iterCounter, size_t nIter, ImageSet& modelSet, ImageSet& dirtySet,
    const aocommon::UVector<const float*>& psfs, bool& reachedMajorThreshold) {
  aocommon::StaticFor<size_t> static_for(aocommon::system::ProcessorCount());
  _staticFor = &static_for;

  reachedMajorThreshold = false;
  if (iterCounter == nIter) return 0.0;

  _modelSet = &modelSet;
  _dirtySet = &dirtySet;
  _psfs = psfs;

  _curBoxXStart = 0;
  _curBoxXEnd = _width;
  _curBoxYStart = 0;
  _curBoxYEnd = _height;

  Image dirty(_width, _height);
  dirtySet.GetLinearIntegrated(dirty);
  Image psf(_width, _height);
  dirtySet.GetIntegratedPSF(psf, psfs);

  int maxScale = IUWTDecomposition::EndScale(std::min(_width, _height));
  int curEndScale = 2;

  // Prepare the PSF for convolutions later on
  Image psfKernel(_width, _height);
  FFTConvolver::PrepareKernel(psfKernel.Data(), psf.Data(), _width, _height,
                              static_for.NThreads());

  std::cout << "Measuring PSF...\n";
  {
    Image convolvedPSF(psf);
    Image scratch(_width, _height);

    FFTConvolver::ConvolveSameSize(_fftwManager, convolvedPSF.Data(),
                                   psfKernel.Data(), _width, _height,
                                   static_for.NThreads());
    measureRMSPerScale(psf.Data(), convolvedPSF.Data(), scratch.Data(),
                       maxScale, _psfResponse);
  }

  ImageSet structureModel(&modelSet.Table(), modelSet.Settings(), _width,
                          _height);

  std::unique_ptr<IUWTDecomposition> iuwt(
      new IUWTDecomposition(curEndScale, _width, _height));

  Image dirtyBeforeIteration;

  float maxValue = 0.0;
  size_t curMinScale = 0;
  reachedMajorThreshold = false;
  bool doContinue = true;
  std::vector<ValComponent> initialComponents;
  do {
    std::cout << "*** Deconvolution iteration " << iterCounter << " ***\n";
    dirtyBeforeIteration = dirty;
    FFTConvolver::PrepareKernel(psfKernel.Data(), psf.Data(), _width, _height,
                                static_for.NThreads());
    std::vector<ValComponent> maxComponents;
    Image scratch(_width, _height);
    bool succeeded = findAndDeconvolveStructure(
        *iuwt, dirty, psf, psfKernel, scratch, structureModel, curEndScale,
        curMinScale, maxComponents);

    if (succeeded) {
      structureModel *= _gain;
      modelSet += structureModel;

      // Calculate: dirty = dirty - structureModel (x) psf
      for (size_t i = 0; i != dirtySet.size(); ++i) {
        scratch.Assign(structureModel[i], structureModel[i] + _width * _height);
        size_t psfIndex = dirtySet.PSFIndex(i);
        FFTConvolver::PrepareKernel(psfKernel.Data(), psfs[psfIndex], _width,
                                    _height, static_for.NThreads());
        FFTConvolver::ConvolveSameSize(_fftwManager, scratch.Data(),
                                       psfKernel.Data(), _width, _height,
                                       static_for.NThreads());
        Subtract(dirtySet[i], scratch);
      }
      dirtySet.GetLinearIntegrated(dirty);

      while (maxComponents.size() > initialComponents.size()) {
        initialComponents.push_back(maxComponents[initialComponents.size()]);
      }
      maxValue = 0.0;
      for (size_t c = 0; c != initialComponents.size(); ++c) {
        std::cout << initialComponents[c].val << " now " << maxComponents[c].val
                  << '\n';
        maxValue = std::max(maxValue, maxComponents[c].val);
        if (std::fabs(maxComponents[c].val) <
            std::fabs(initialComponents[c].val) * (1.0 - _mGain)) {
          std::cout << "Scale " << c << " reached mGain (starting level: "
                    << initialComponents[c].val
                    << ", now: " << maxComponents[c].val << ").\n";
          reachedMajorThreshold = true;
        }
      }
      if (reachedMajorThreshold) break;
    } else {
      if (int(curMinScale) + 1 < curEndScale) {
        ++curMinScale;
        std::cout << "=> Min scale now " << curMinScale << '\n';
      } else {
        curMinScale = 0;
        if (curEndScale != maxScale) {
          ++curEndScale;
          std::cout << "=> Scale now " << curEndScale << ".\n";
          iuwt.reset(new IUWTDecomposition(curEndScale, _width, _height));
        } else {
          std::cout << "Max scale reached: finished all scales, quiting.\n";
          doContinue = false;
        }
      }
      dirty = dirtyBeforeIteration;
    }

    ++iterCounter;
  } while (iterCounter != nIter && doContinue);
  return maxValue;
}
