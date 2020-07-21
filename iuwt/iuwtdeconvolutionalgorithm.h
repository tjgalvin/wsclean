#ifndef IUWT_DECONVOLUTION_ALGORITHM_H
#define IUWT_DECONVOLUTION_ALGORITHM_H

#include <aocommon/uvector.h>

#include "iuwtdecomposition.h"
#include "imageanalysis.h"

#include <vector>

class IUWTDeconvolutionAlgorithm {
 public:
  IUWTDeconvolutionAlgorithm(class FFTWManager& fftwManager, size_t width,
                             size_t height, double gain, double mGain,
                             double cleanBorder, bool allowNegativeComponents,
                             const bool* mask, double absoluteThreshold,
                             double thresholdSigmaLevel = 4.0,
                             double tolerance = 0.75, bool useSNRTest = true);

  double PerformMajorIteration(size_t& iterCounter, size_t nIter,
                               class ImageSet& modelSet,
                               class ImageSet& dirtySet,
                               const aocommon::UVector<const double*>& psfs,
                               bool& reachedMajorThreshold);

  void Subtract(double* dest, const aocommon::UVector<double>& rhs);
  void Subtract(aocommon::UVector<double>& dest,
                const aocommon::UVector<double>& rhs) {
    Subtract(dest.data(), rhs);
  }

 private:
  struct ValComponent {
    ValComponent() {}
    ValComponent(size_t _x, size_t _y, int _scale, double _val = 0.0)
        : x(_x), y(_y), scale(_scale), val(_val) {}

    std::string ToString() const {
      std::ostringstream str;
      str << x << ',' << y << ", scale " << scale;
      return str.str();
    }

    size_t x, y;
    int scale;
    double val;
  };

  struct ScaleResponse {
    double rms, peakResponse, peakResponseToNextScale, convolvedPeakResponse;
    double bMaj, bMin, bPA;
    size_t convolvedArea;
  };

  double getMaxAbsWithoutMask(const aocommon::UVector<double>& data, size_t& x,
                              size_t& y, size_t width);
  double getMaxAbsWithMask(const aocommon::UVector<double>& data, size_t& x,
                           size_t& y, size_t width);
  double getMaxAbs(const aocommon::UVector<double>& data, size_t& x, size_t& y,
                   size_t width) {
    if (_mask == 0)
      return getMaxAbsWithoutMask(data, x, y, width);
    else
      return getMaxAbsWithMask(data, x, y, width);
  }

  void measureRMSPerScale(const double* image, const double* convolvedImage,
                          double* scratch, size_t endScale,
                          std::vector<ScaleResponse>& psfResponse);

  double mad(const double* dest);

  double dotProduct(const aocommon::UVector<double>& lhs,
                    const aocommon::UVector<double>& rhs);

  void factorAdd(double* dest, const double* rhs, double factor, size_t width,
                 size_t height);

  void factorAdd(aocommon::UVector<double>& dest,
                 const aocommon::UVector<double>& rhs, double factor);

  void boundingBox(size_t& x1, size_t& y1, size_t& x2, size_t& y2,
                   const aocommon::UVector<double>& image, size_t width,
                   size_t height);

  void adjustBox(size_t& x1, size_t& y1, size_t& x2, size_t& y2, size_t width,
                 size_t height, int endScale);

  void trim(aocommon::UVector<double>& dest, const double* source,
            size_t oldWidth, size_t oldHeight, size_t x1, size_t y1, size_t x2,
            size_t y2);

  void trim(aocommon::UVector<double>& dest,
            const aocommon::UVector<double>& source, size_t oldWidth,
            size_t oldHeight, size_t x1, size_t y1, size_t x2, size_t y2) {
    trim(dest, source.data(), oldWidth, oldHeight, x1, y1, x2, y2);
  }

  void trimPsf(aocommon::UVector<double>& dest, const double* source,
               size_t oldWidth, size_t oldHeight, size_t newWidth,
               size_t newHeight) {
    trim(dest, source, oldWidth, oldHeight, (oldWidth - newWidth) / 2,
         (oldHeight - newHeight) / 2, (oldWidth + newWidth) / 2,
         (oldHeight + newHeight) / 2);
  }

  void untrim(aocommon::UVector<double>& image, size_t width, size_t height,
              size_t x1, size_t y1, size_t x2, size_t y2);

  double sum(const aocommon::UVector<double>& img) const;

  double snr(const IUWTDecomposition& noisyImg,
             const IUWTDecomposition& model) const;

  double rmsDiff(const aocommon::UVector<double>& a,
                 const aocommon::UVector<double>& b);

  double rms(const aocommon::UVector<double>& image);

  bool runConjugateGradient(IUWTDecomposition& iuwt, const IUWTMask& mask,
                            aocommon::UVector<double>& maskedDirty,
                            aocommon::UVector<double>& structureModel,
                            aocommon::UVector<double>& scratch,
                            const aocommon::UVector<double>& psfKernel,
                            size_t width, size_t height);

  bool fillAndDeconvolveStructure(IUWTDecomposition& iuwt,
                                  aocommon::UVector<double>& dirty,
                                  class ImageSet& structureModelFull,
                                  aocommon::UVector<double>& scratch,
                                  const aocommon::UVector<double>& psf,
                                  const aocommon::UVector<double>& psfKernel,
                                  size_t curEndScale, size_t curMinScale,
                                  size_t width, size_t height,
                                  const aocommon::UVector<double>& thresholds,
                                  const ImageAnalysis::Component& maxComp,
                                  bool allowTrimming, const bool* priorMask);

  bool findAndDeconvolveStructure(IUWTDecomposition& iuwt,
                                  aocommon::UVector<double>& dirty,
                                  const aocommon::UVector<double>& psf,
                                  const aocommon::UVector<double>& psfKernel,
                                  aocommon::UVector<double>& scratch,
                                  class ImageSet& structureModelFull,
                                  size_t curEndScale, size_t curMinScale,
                                  std::vector<ValComponent>& maxComponents);

  void performSubImageFitAll(IUWTDecomposition& iuwt, const IUWTMask& mask,
                             const aocommon::UVector<double>& structureModel,
                             aocommon::UVector<double>& scratchA,
                             aocommon::UVector<double>& scratchB,
                             const ImageAnalysis::Component& maxComp,
                             ImageSet& fittedModel, const double* psf,
                             const aocommon::UVector<double>& dirty);

  void performSubImageFitSingle(IUWTDecomposition& iuwt, const IUWTMask& mask,
                                const aocommon::UVector<double>& structureModel,
                                aocommon::UVector<double>& scratchB,
                                const ImageAnalysis::Component& maxComp,
                                const double* psf,
                                aocommon::UVector<double>& subDirty,
                                double* fittedSubModel,
                                aocommon::UVector<double>& correctionFactors);

  double performSubImageComponentFitBoxed(
      IUWTDecomposition& iuwt, const IUWTMask& mask,
      const std::vector<ImageAnalysis::Component2D>& area,
      aocommon::UVector<double>& scratch,
      aocommon::UVector<double>& maskedDirty, const double* psf,
      const aocommon::UVector<double>& psfKernel, size_t x1, size_t y1,
      size_t x2, size_t y2);

  double performSubImageComponentFit(
      IUWTDecomposition& iuwt, const IUWTMask& mask,
      const std::vector<ImageAnalysis::Component2D>& area,
      aocommon::UVector<double>& scratch,
      aocommon::UVector<double>& maskedDirty,
      const aocommon::UVector<double>& psfKernel, size_t xOffset,
      size_t yOffset);

  double centralPeak(const aocommon::UVector<double>& data) {
    return data[_width / 2 + (_height / 2) * _width];
  }
  void constrainedPSFConvolve(double* image, const double* psf, size_t width,
                              size_t height);

  bool extractPointSources(const IUWTDecomposition& iuwt, const IUWTMask& mask,
                           const double* dirty, double* model);

  class FFTWManager& _fftwManager;
  size_t _width, _height;
  size_t _curBoxXStart, _curBoxXEnd;
  size_t _curBoxYStart, _curBoxYEnd;
  double _gain, _mGain, _cleanBorder;
  const bool* _mask;
  double _absoluteThreshold, _thresholdSigmaLevel, _tolerance;
  double _psfMaj, _psfMin, _psfPA, _psfVolume;
  aocommon::UVector<double> _rmses;
  FitsWriter _writer;
  std::vector<ScaleResponse> _psfResponse;
  bool _allowNegativeComponents, _useSNRTest;
  class ImageSet* _modelSet;
  class ImageSet* _dirtySet;
  aocommon::UVector<const double*> _psfs;
  class ThreadPool* _threadPool;
};

#endif
