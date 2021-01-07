#ifndef IUWT_DECONVOLUTION_ALGORITHM_H
#define IUWT_DECONVOLUTION_ALGORITHM_H

#include <aocommon/uvector.h>

#include "iuwtdecomposition.h"
#include "imageanalysis.h"

#include "../structures/image.h"

#include <vector>

class IUWTDeconvolutionAlgorithm {
 public:
  IUWTDeconvolutionAlgorithm(class FFTWManager& fftwManager, size_t width,
                             size_t height, float gain, float mGain,
                             float cleanBorder, bool allowNegativeComponents,
                             const bool* mask, float absoluteThreshold,
                             float thresholdSigmaLevel = 4.0,
                             float tolerance = 0.75, bool useSNRTest = true);

  float PerformMajorIteration(size_t& iterCounter, size_t nIter,
                              class ImageSet& modelSet,
                              class ImageSet& dirtySet,
                              const aocommon::UVector<const float*>& psfs,
                              bool& reachedMajorThreshold);

  void Subtract(float* dest, const ImageF& rhs);
  void Subtract(ImageF& dest, const ImageF& rhs) { Subtract(dest.data(), rhs); }

 private:
  struct ValComponent {
    ValComponent() {}
    ValComponent(size_t _x, size_t _y, int _scale, float _val = 0.0)
        : x(_x), y(_y), scale(_scale), val(_val) {}

    std::string ToString() const {
      std::ostringstream str;
      str << x << ',' << y << ", scale " << scale;
      return str.str();
    }

    size_t x, y;
    int scale;
    float val;
  };

  struct ScaleResponse {
    float rms, peakResponse, peakResponseToNextScale, convolvedPeakResponse;
    double bMaj, bMin, bPA;
    size_t convolvedArea;
  };

  float getMaxAbsWithoutMask(const ImageF& data, size_t& x, size_t& y,
                             size_t width);
  float getMaxAbsWithMask(const ImageF& data, size_t& x, size_t& y,
                          size_t width);
  float getMaxAbs(const ImageF& data, size_t& x, size_t& y, size_t width) {
    if (_mask == 0)
      return getMaxAbsWithoutMask(data, x, y, width);
    else
      return getMaxAbsWithMask(data, x, y, width);
  }

  void measureRMSPerScale(const float* image, const float* convolvedImage,
                          float* scratch, size_t endScale,
                          std::vector<ScaleResponse>& psfResponse);

  float mad(const float* dest);

  float dotProduct(const ImageF& lhs, const ImageF& rhs);

  void factorAdd(float* dest, const float* rhs, float factor, size_t width,
                 size_t height);

  void factorAdd(ImageF& dest, const ImageF& rhs, float factor);

  void boundingBox(size_t& x1, size_t& y1, size_t& x2, size_t& y2,
                   const ImageF& image, size_t width, size_t height);

  void adjustBox(size_t& x1, size_t& y1, size_t& x2, size_t& y2, size_t width,
                 size_t height, int endScale);

  void trim(ImageF& dest, const float* source, size_t oldWidth,
            size_t oldHeight, size_t x1, size_t y1, size_t x2, size_t y2);

  void trim(ImageF& dest, const ImageF& source, size_t oldWidth,
            size_t oldHeight, size_t x1, size_t y1, size_t x2, size_t y2) {
    trim(dest, source.data(), oldWidth, oldHeight, x1, y1, x2, y2);
  }

  void trimPsf(ImageF& dest, const float* source, size_t oldWidth,
               size_t oldHeight, size_t newWidth, size_t newHeight) {
    trim(dest, source, oldWidth, oldHeight, (oldWidth - newWidth) / 2,
         (oldHeight - newHeight) / 2, (oldWidth + newWidth) / 2,
         (oldHeight + newHeight) / 2);
  }

  void untrim(ImageF& image, size_t width, size_t height, size_t x1, size_t y1,
              size_t x2, size_t y2);

  float sum(const ImageF& img) const;

  float snr(const IUWTDecomposition& noisyImg,
            const IUWTDecomposition& model) const;

  float rmsDiff(const ImageF& a, const ImageF& b);

  float rms(const ImageF& image);

  bool runConjugateGradient(IUWTDecomposition& iuwt, const IUWTMask& mask,
                            ImageF& maskedDirty, ImageF& structureModel,
                            ImageF& scratch, const ImageF& psfKernel,
                            size_t width, size_t height);

  bool fillAndDeconvolveStructure(IUWTDecomposition& iuwt, ImageF& dirty,
                                  class ImageSet& structureModelFull,
                                  ImageF& scratch, const ImageF& psf,
                                  const ImageF& psfKernel, size_t curEndScale,
                                  size_t curMinScale, size_t width,
                                  size_t height,
                                  const aocommon::UVector<float>& thresholds,
                                  const ImageAnalysis::Component& maxComp,
                                  bool allowTrimming, const bool* priorMask);

  bool findAndDeconvolveStructure(IUWTDecomposition& iuwt, ImageF& dirty,
                                  const ImageF& psf, const ImageF& psfKernel,
                                  ImageF& scratch,
                                  class ImageSet& structureModelFull,
                                  size_t curEndScale, size_t curMinScale,
                                  std::vector<ValComponent>& maxComponents);

  void performSubImageFitAll(IUWTDecomposition& iuwt, const IUWTMask& mask,
                             const ImageF& structureModel, ImageF& scratchA,
                             ImageF& scratchB,
                             const ImageAnalysis::Component& maxComp,
                             ImageSet& fittedModel, const float* psf,
                             const ImageF& dirty);

  void performSubImageFitSingle(IUWTDecomposition& iuwt, const IUWTMask& mask,
                                const ImageF& structureModel, ImageF& scratchB,
                                const ImageAnalysis::Component& maxComp,
                                const float* psf, ImageF& subDirty,
                                float* fittedSubModel,
                                aocommon::UVector<float>& correctionFactors);

  float performSubImageComponentFitBoxed(
      IUWTDecomposition& iuwt, const IUWTMask& mask,
      const std::vector<ImageAnalysis::Component2D>& area, ImageF& scratch,
      ImageF& maskedDirty, const float* psf, const ImageF& psfKernel, size_t x1,
      size_t y1, size_t x2, size_t y2);

  float performSubImageComponentFit(
      IUWTDecomposition& iuwt, const IUWTMask& mask,
      const std::vector<ImageAnalysis::Component2D>& area, ImageF& scratch,
      ImageF& maskedDirty, const ImageF& psfKernel, size_t xOffset,
      size_t yOffset);

  float centralPeak(const ImageF& data) {
    return data[_width / 2 + (_height / 2) * _width];
  }
  void constrainedPSFConvolve(float* image, const float* psf, size_t width,
                              size_t height);

  class FFTWManager& _fftwManager;
  size_t _width, _height;
  size_t _curBoxXStart, _curBoxXEnd;
  size_t _curBoxYStart, _curBoxYEnd;
  float _gain, _mGain, _cleanBorder;
  const bool* _mask;
  float _absoluteThreshold, _thresholdSigmaLevel, _tolerance;
  double _psfMaj, _psfMin, _psfPA, _psfVolume;
  aocommon::UVector<float> _rmses;
  FitsWriter _writer;
  std::vector<ScaleResponse> _psfResponse;
  bool _allowNegativeComponents, _useSNRTest;
  class ImageSet* _modelSet;
  class ImageSet* _dirtySet;
  aocommon::UVector<const float*> _psfs;
  class ThreadPool* _threadPool;
};

#endif
