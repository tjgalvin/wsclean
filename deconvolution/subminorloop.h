#ifndef SUB_MINOR_LOOP_H
#define SUB_MINOR_LOOP_H

#include <cstring>
#include <optional>
#include <vector>

#include "../deconvolution/imageset.h"

#include <aocommon/image.h>
#include <aocommon/logger.h>

/**
 * In multi-scale, a subminor optimized loop looks like this:
 *
 * IterateAndMakeModel():
 * - Make a set S with positions of all the components larger than 'threshold',
 * which are also in the mask
 * - Find the largest component in S
 * Loop {
 * - Measure the largest component per frequency (from S)
 * - Store the model component in S
 * - Subtract this component multiplied with the twice convolved PSF and gain
 * from all components in S (per individual image)
 * - Find the new largest component in S
 * }
 *
 * CorrectResidualDirty():
 * For each individual image {
 * - Put the model components from S onto a full image (using
 * GetFullIndividualModel())
 * - Convolve the model with the SingleConvolvedPSF
 * - Subtract the convolved model from the residual
 * }
 *
 * Finalization:
 * - Put the model components from S onto a full image (using
 * GetFullIndividualModel())
 * - Convolve the model image with the scale kernel
 * - Add the model components to the full model
 *
 * A subminor loop has some correspondance with the so-called Clark
 * optimization. However, this implementation has some differences, e.g. by
 * collecting a list of threshold components prior of entering the subminor
 * loop.
 */

class SubMinorModel {
 public:
  SubMinorModel(size_t width, size_t /*height*/) : _width(width) {}

  void AddPosition(size_t x, size_t y) {
    _positions.push_back(std::make_pair(x, y));
  }

  /**
   * Return number of selected pixels.
   */
  size_t size() const { return _positions.size(); }

  void MakeSets(const ImageSet& templateSet);
  void MakeRMSFactorImage(aocommon::Image& rmsFactorImage);

  ImageSet& Residual() { return *_residual; }
  const ImageSet& Residual() const { return *_residual; }

  ImageSet& Model() { return *_model; }
  const ImageSet& Model() const { return *_model; }

  size_t X(size_t index) const { return _positions[index].first; }
  size_t Y(size_t index) const { return _positions[index].second; }
  size_t FullIndex(size_t index) const { return X(index) + Y(index) * _width; }
  template <bool AllowNegatives>
  size_t GetMaxComponent(aocommon::Image& scratch, float& maxValue) const;
  size_t GetMaxComponent(aocommon::Image& scratch, float& maxValue,
                         bool allowNegatives) const {
    if (allowNegatives)
      return GetMaxComponent<true>(scratch, maxValue);
    else
      return GetMaxComponent<false>(scratch, maxValue);
  }

 private:
  std::vector<std::pair<size_t, size_t>> _positions;
  std::unique_ptr<ImageSet> _residual, _model;
  aocommon::Image _rmsFactorImage;
  size_t _width;
};

class SubMinorLoop {
 public:
  SubMinorLoop(size_t width, size_t height, size_t convolutionWidth,
               size_t convolutionHeight, aocommon::LogReceiver& logReceiver)
      : _width(width),
        _height(height),
        _paddedWidth(convolutionWidth),
        _paddedHeight(convolutionHeight),
        _threshold(0.0),
        _consideredPixelThreshold(0.0),
        _gain(0.0),
        _horizontalBorder(0),
        _verticalBorder(0),
        _currentIteration(0),
        _maxIterations(0),
        _allowNegativeComponents(true),
        _stopOnNegativeComponent(false),
        _mask(nullptr),
        _fitter(nullptr),
        _subMinorModel(width, height),
        _fluxCleaned(0.0),
        _logReceiver(logReceiver),
        _threadCount(1) {}

  /**
   * @param threshold The threshold to which this subminor run should clean
   * @param consideredPixelThreshold The threshold that is used to determine
   * whether a pixel is considered. Typically, this is similar to threshold, but
   * it can be set lower if it is important that all peak values are below the
   * threshold, as otherwise some pixels might not be considered but get
   * increased by the cleaning, thereby stay above the threshold. This is
   * important for making multi-scale clean efficient near a stopping threshold.
   */
  void SetThreshold(float threshold, float consideredPixelThreshold) {
    _threshold = threshold;
    _consideredPixelThreshold = consideredPixelThreshold;
  }

  void SetIterationInfo(size_t currentIteration, size_t maxIterations) {
    _currentIteration = currentIteration;
    _maxIterations = maxIterations;
  }

  void SetGain(float gain) { _gain = gain; }

  void SetAllowNegativeComponents(bool allowNegativeComponents) {
    _allowNegativeComponents = allowNegativeComponents;
  }

  void SetStopOnNegativeComponent(bool stopOnNegativeComponent) {
    _stopOnNegativeComponent = stopOnNegativeComponent;
  }

  void SetSpectralFitter(const schaapcommon::fitters::SpectralFitter* fitter) {
    _fitter = fitter;
  }

  void SetCleanBorders(size_t horizontalBorder, size_t verticalBorder) {
    _horizontalBorder = horizontalBorder;
    _verticalBorder = verticalBorder;
  }

  void SetMask(const bool* mask) { _mask = mask; }

  void SetRMSFactorImage(const aocommon::Image& image) {
    _rmsFactorImage = image;
  }

  void SetThreadCount(size_t threadCount) { _threadCount = threadCount; }

  size_t CurrentIteration() const { return _currentIteration; }

  float FluxCleaned() const { return _fluxCleaned; }

  std::optional<float> Run(
      ImageSet& convolvedResidual,
      const std::vector<aocommon::Image>& twiceConvolvedPsfs);

  /**
   * The produced model is convolved with the given psf, and the result is
   * subtracted from the given residual image. To be called after Run(). After
   * this method, the residual will hold the result of the subminor loop run.
   * scratchA and scratchB need to be able to store the full padded image
   * (_untrimmedWidth x _untrimmedHeight). scratchC only needs to store the
   * trimmed size (_width x _height).
   */
  void CorrectResidualDirty(float* scratchA, float* scratchB, float* scratchC,
                            size_t imageIndex, float* residual,
                            const float* singleConvolvedPsf) const;

  void GetFullIndividualModel(size_t imageIndex,
                              float* individualModelImg) const;

  void UpdateAutoMask(bool* mask) const;

  void UpdateComponentList(class ComponentList& list, size_t scaleIndex) const;

 private:
  void findPeakPositions(ImageSet& convolvedResidual);

  size_t _width, _height, _paddedWidth, _paddedHeight;
  float _threshold, _consideredPixelThreshold, _gain;
  size_t _horizontalBorder, _verticalBorder;
  size_t _currentIteration, _maxIterations;
  bool _allowNegativeComponents, _stopOnNegativeComponent;
  const bool* _mask;
  const schaapcommon::fitters::SpectralFitter* _fitter;
  SubMinorModel _subMinorModel;
  float _fluxCleaned;
  aocommon::Image _rmsFactorImage;
  aocommon::LogReceiver& _logReceiver;
  size_t _threadCount;
};

#endif
