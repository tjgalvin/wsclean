#ifndef IUWT_DECONVOLUTION_ALGORITHM_H
#define IUWT_DECONVOLUTION_ALGORITHM_H

#include <aocommon/staticfor.h>
#include <aocommon/uvector.h>
#include <aocommon/fits/fitswriter.h>

#include "iuwtdecomposition.h"
#include "imageanalysis.h"

#include <vector>

namespace aocommon {
template <typename NumT>
class ImageBase;
using Image = ImageBase<float>;
}  // namespace aocommon

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
                              const std::vector<aocommon::Image>& psfs,
                              bool& reachedMajorThreshold);

  void Subtract(float* dest, const aocommon::Image& rhs);
  void Subtract(aocommon::Image& dest, const aocommon::Image& rhs) {
    Subtract(dest.Data(), rhs);
  }

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

  float getMaxAbsWithoutMask(const aocommon::Image& data, size_t& x, size_t& y,
                             size_t width);
  float getMaxAbsWithMask(const aocommon::Image& data, size_t& x, size_t& y,
                          size_t width);
  float getMaxAbs(const aocommon::Image& data, size_t& x, size_t& y,
                  size_t width) {
    if (_mask == nullptr)
      return getMaxAbsWithoutMask(data, x, y, width);
    else
      return getMaxAbsWithMask(data, x, y, width);
  }

  void measureRMSPerScale(const float* image, const float* convolvedImage,
                          float* scratch, size_t endScale,
                          std::vector<ScaleResponse>& psfResponse);

  float mad(const float* dest);

  float dotProduct(const aocommon::Image& lhs, const aocommon::Image& rhs);

  void factorAdd(float* dest, const float* rhs, float factor, size_t width,
                 size_t height);

  void factorAdd(aocommon::Image& dest, const aocommon::Image& rhs,
                 float factor);

  void boundingBox(size_t& x1, size_t& y1, size_t& x2, size_t& y2,
                   const aocommon::Image& image, size_t width, size_t height);

  void adjustBox(size_t& x1, size_t& y1, size_t& x2, size_t& y2, size_t width,
                 size_t height, int endScale);

  void trim(aocommon::Image& dest, const float* source, size_t oldWidth,
            size_t x1, size_t y1, size_t x2, size_t y2);

  void trim(aocommon::Image& dest, const aocommon::Image& source, size_t x1,
            size_t y1, size_t x2, size_t y2) {
    trim(dest, source.Data(), source.Width(), x1, y1, x2, y2);
  }

  void trimPsf(aocommon::Image& dest, const aocommon::Image& source,
               size_t newWidth, size_t newHeight) {
    const size_t oldWidth = source.Width();
    const size_t oldHeight = source.Height();
    trim(dest, source, (oldWidth - newWidth) / 2, (oldHeight - newHeight) / 2,
         (oldWidth + newWidth) / 2, (oldHeight + newHeight) / 2);
  }

  void untrim(aocommon::Image& image, size_t width, size_t height, size_t x1,
              size_t y1, size_t x2, size_t y2);

  float sum(const aocommon::Image& img) const;

  float snr(const IUWTDecomposition& noisyImg,
            const IUWTDecomposition& model) const;

  float rmsDiff(const aocommon::Image& a, const aocommon::Image& b);

  float rms(const aocommon::Image& image);

  bool runConjugateGradient(IUWTDecomposition& iuwt, const IUWTMask& mask,
                            aocommon::Image& maskedDirty,
                            aocommon::Image& structureModel,
                            aocommon::Image& scratch,
                            const aocommon::Image& psfKernel, size_t width,
                            size_t height);

  bool fillAndDeconvolveStructure(
      IUWTDecomposition& iuwt, aocommon::Image& dirty,
      class ImageSet& structureModelFull, aocommon::Image& scratch,
      const aocommon::Image& psf, const aocommon::Image& psfKernel,
      const std::vector<aocommon::Image>& psfs, size_t curEndScale,
      size_t curMinScale, size_t width, size_t height,
      const aocommon::UVector<float>& thresholds,
      const ImageAnalysis::Component& maxComp, bool allowTrimming,
      const bool* priorMask);

  bool findAndDeconvolveStructure(
      IUWTDecomposition& iuwt, aocommon::Image& dirty,
      const aocommon::Image& psf, const aocommon::Image& psfKernel,
      const std::vector<aocommon::Image>& psfs, aocommon::Image& scratch,
      class ImageSet& structureModelFull, size_t curEndScale,
      size_t curMinScale, std::vector<ValComponent>& maxComponents);

  void performSubImageFitAll(IUWTDecomposition& iuwt, const IUWTMask& mask,
                             const aocommon::Image& structureModel,
                             aocommon::Image& scratchA,
                             aocommon::Image& scratchB,
                             const ImageAnalysis::Component& maxComp,
                             ImageSet& fittedModel, const aocommon::Image& psf,
                             const std::vector<aocommon::Image>& psfs,
                             const aocommon::Image& dirty);

  void performSubImageFitSingle(IUWTDecomposition& iuwt, const IUWTMask& mask,
                                const aocommon::Image& structureModel,
                                aocommon::Image& scratchB,
                                const ImageAnalysis::Component& maxComp,
                                const aocommon::Image& psf,
                                aocommon::Image& subDirty,
                                float* fittedSubModel,
                                aocommon::UVector<float>& correctionFactors);

  float performSubImageComponentFitBoxed(
      IUWTDecomposition& iuwt, const IUWTMask& mask,
      const std::vector<ImageAnalysis::Component2D>& area,
      aocommon::Image& scratch, aocommon::Image& maskedDirty,
      const aocommon::Image& psf, const aocommon::Image& psfKernel, size_t x1,
      size_t y1, size_t x2, size_t y2);

  float performSubImageComponentFit(
      IUWTDecomposition& iuwt, const IUWTMask& mask,
      const std::vector<ImageAnalysis::Component2D>& area,
      aocommon::Image& scratch, aocommon::Image& maskedDirty,
      const aocommon::Image& psfKernel, size_t xOffset, size_t yOffset);

  float centralPeak(const aocommon::Image& data) {
    return data[_width / 2 + (_height / 2) * _width];
  }

  class FFTWManager& _fftwManager;
  size_t _width, _height;
  size_t _curBoxXStart, _curBoxXEnd;
  size_t _curBoxYStart, _curBoxYEnd;
  float _gain, _mGain, _cleanBorder;
  const bool* _mask;
  float _absoluteThreshold, _thresholdSigmaLevel, _tolerance;
  double _psfMaj, _psfMin, _psfPA, _psfVolume;
  aocommon::UVector<float> _rmses;
  aocommon::FitsWriter _writer;
  std::vector<ScaleResponse> _psfResponse;
  bool _allowNegativeComponents, _useSNRTest;
  class ImageSet* _modelSet;
  class ImageSet* _dirtySet;
  aocommon::StaticFor<size_t>* _staticFor;
};

#endif
