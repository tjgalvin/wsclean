#include "genericclean.h"

#include "subminorloop.h"

#include <aocommon/lane.h>

#include "../multiscale/threadeddeconvolutiontools.h"

#include "../units/fluxdensity.h"

#include "../wsclean/logger.h"

#include <boost/thread/thread.hpp>

GenericClean::GenericClean(class FFTWManager& fftwManager,
                           bool useSubMinorOptimization)
    : _convolutionPadding(1.1),
      _useSubMinorOptimization(useSubMinorOptimization),
      _fftwManager(fftwManager) {}

double GenericClean::ExecuteMajorIteration(
    ImageSet& dirtySet, ImageSet& modelSet,
    const aocommon::UVector<const double*>& psfs, size_t width, size_t height,
    bool& reachedMajorThreshold) {
  const size_t iterationCounterAtStart = _iterationNumber;
  if (_stopOnNegativeComponent) _allowNegativeComponents = true;
  _width = width;
  _height = height;
  _convolutionWidth = ceil(_convolutionPadding * _width);
  _convolutionHeight = ceil(_convolutionPadding * _height);
  if (_convolutionWidth % 2 != 0) ++_convolutionWidth;
  if (_convolutionHeight % 2 != 0) ++_convolutionHeight;

  Image integrated(width, height),
      scratchA(_convolutionWidth, _convolutionHeight),
      scratchB(_convolutionWidth, _convolutionHeight);
  dirtySet.GetLinearIntegrated(integrated.data());
  size_t componentX = 0, componentY = 0;
  boost::optional<double> maxValue =
      findPeak(integrated.data(), scratchA.data(), componentX, componentY);
  if (!maxValue) {
    _logReceiver->Info << "No peak found.\n";
    reachedMajorThreshold = false;
    return 0.0;
  }
  _logReceiver->Info << "Initial peak: "
                     << peakDescription(integrated.data(), componentX,
                                        componentY)
                     << '\n';
  double firstThreshold = this->_threshold;
  double majorIterThreshold = std::max(
      MajorIterThreshold(), std::fabs(*maxValue) * (1.0 - this->_mGain));
  if (majorIterThreshold > firstThreshold) {
    firstThreshold = majorIterThreshold;
    _logReceiver->Info << "Next major iteration at: "
                       << FluxDensity::ToNiceString(majorIterThreshold) << '\n';
  } else if (this->_mGain != 1.0) {
    _logReceiver->Info
        << "Major iteration threshold reached global threshold of "
        << FluxDensity::ToNiceString(this->_threshold)
        << ": final major iteration.\n";
  }

  if (_useSubMinorOptimization) {
    size_t startIteration = _iterationNumber;
    SubMinorLoop subMinorLoop(_width, _height, _convolutionWidth,
                              _convolutionHeight, *_logReceiver);
    subMinorLoop.SetIterationInfo(_iterationNumber, MaxNIter());
    subMinorLoop.SetThreshold(firstThreshold, firstThreshold * 0.99);
    subMinorLoop.SetGain(Gain());
    subMinorLoop.SetAllowNegativeComponents(AllowNegativeComponents());
    subMinorLoop.SetStopOnNegativeComponent(StopOnNegativeComponents());
    subMinorLoop.SetSpectralFitter(&Fitter());
    if (!_rmsFactorImage.empty())
      subMinorLoop.SetRMSFactorImage(_rmsFactorImage);
    if (_cleanMask) subMinorLoop.SetMask(_cleanMask);
    const size_t horBorderSize = round(_width * CleanBorderRatio()),
                 vertBorderSize = round(_height * CleanBorderRatio());
    subMinorLoop.SetCleanBorders(horBorderSize, vertBorderSize);

    maxValue = subMinorLoop.Run(dirtySet, psfs);

    _iterationNumber = subMinorLoop.CurrentIteration();

    _logReceiver->Info
        << "Performed " << _iterationNumber << " iterations in total, "
        << (_iterationNumber - startIteration)
        << " in this major iteration with sub-minor optimization.\n";

    for (size_t imageIndex = 0; imageIndex != dirtySet.size(); ++imageIndex) {
      // TODO this can be multi-threaded if each thread has its own temporaries
      const double* psf = psfs[dirtySet.PSFIndex(imageIndex)];
      subMinorLoop.CorrectResidualDirty(_fftwManager, scratchA.data(),
                                        scratchB.data(), integrated.data(),
                                        imageIndex, dirtySet[imageIndex], psf);

      subMinorLoop.GetFullIndividualModel(imageIndex, scratchA.data());
      double* model = modelSet[imageIndex];
      for (size_t i = 0; i != _width * _height; ++i)
        model[i] += scratchA.data()[i];
    }
  } else {
    ThreadedDeconvolutionTools tools(_threadCount);
    size_t peakIndex = componentX + componentY * _width;

    aocommon::UVector<double> peakValues(dirtySet.size());

    while (maxValue && fabs(*maxValue) > firstThreshold &&
           this->_iterationNumber < this->_maxIter &&
           !(maxValue < 0.0 && this->_stopOnNegativeComponent)) {
      if (this->_iterationNumber <= 10 ||
          (this->_iterationNumber <= 100 && this->_iterationNumber % 10 == 0) ||
          (this->_iterationNumber <= 1000 &&
           this->_iterationNumber % 100 == 0) ||
          this->_iterationNumber % 1000 == 0)
        _logReceiver->Info << "Iteration " << this->_iterationNumber << ": "
                           << peakDescription(integrated.data(), componentX,
                                              componentY)
                           << '\n';

      for (size_t i = 0; i != dirtySet.size(); ++i)
        peakValues[i] = dirtySet[i][peakIndex];

      PerformSpectralFit(peakValues.data());

      for (size_t i = 0; i != dirtySet.size(); ++i) {
        peakValues[i] *= this->_gain;
        modelSet[i][peakIndex] += peakValues[i];

        size_t psfIndex = dirtySet.PSFIndex(i);

        tools.SubtractImage(dirtySet[i], psfs[psfIndex], width, height,
                            componentX, componentY, peakValues[i]);
      }

      dirtySet.GetSquareIntegrated(integrated.data(), scratchA.data());
      maxValue =
          findPeak(integrated.data(), scratchA.data(), componentX, componentY);

      peakIndex = componentX + componentY * _width;

      ++this->_iterationNumber;
    }
  }
  if (maxValue) {
    _logReceiver->Info << "Stopped on peak "
                       << FluxDensity::ToNiceString(*maxValue) << ", because ";
    bool maxIterReached = _iterationNumber >= MaxNIter(),
         finalThresholdReached =
             std::fabs(*maxValue) <= _threshold || maxValue == 0.0,
         negativeReached = maxValue < 0.0 && this->_stopOnNegativeComponent,
         mgainReached = std::fabs(*maxValue) <= majorIterThreshold,
         didWork = (_iterationNumber - iterationCounterAtStart) != 0;

    if (maxIterReached)
      _logReceiver->Info << "maximum number of iterations was reached.\n";
    else if (finalThresholdReached)
      _logReceiver->Info << "the threshold was reached.\n";
    else if (negativeReached)
      _logReceiver->Info << "a negative component was found.\n";
    else if (!didWork)
      _logReceiver->Info << "no iterations could be performed.\n";
    else
      _logReceiver->Info << "the minor-loop threshold was reached. Continuing "
                            "cleaning after inversion/prediction round.\n";
    reachedMajorThreshold =
        mgainReached && didWork && !negativeReached && !finalThresholdReached;
    return *maxValue;
  } else {
    _logReceiver->Info << "Deconvolution aborted.\n";
    reachedMajorThreshold = false;
    return 0.0;
  }
}

std::string GenericClean::peakDescription(const double* image, size_t& x,
                                          size_t& y) {
  std::ostringstream str;
  size_t index = x + y * _width;
  double peak = image[index];
  str << FluxDensity::ToNiceString(peak) << " at " << x << "," << y;
  return str.str();
}

boost::optional<double> GenericClean::findPeak(const double* image,
                                               double* scratch, size_t& x,
                                               size_t& y) {
  if (_rmsFactorImage.empty()) {
    if (_cleanMask == 0)
      return SimpleClean::FindPeak(image, _width, _height, x, y,
                                   _allowNegativeComponents, 0, _height,
                                   _cleanBorderRatio);
    else
      return SimpleClean::FindPeakWithMask(image, _width, _height, x, y,
                                           _allowNegativeComponents, 0, _height,
                                           _cleanMask, _cleanBorderRatio);
  } else {
    for (size_t i = 0; i != _width * _height; ++i) {
      scratch[i] = image[i] * _rmsFactorImage[i];
    }
    if (_cleanMask == 0)
      return SimpleClean::FindPeak(scratch, _width, _height, x, y,
                                   _allowNegativeComponents, 0, _height,
                                   _cleanBorderRatio);
    else
      return SimpleClean::FindPeakWithMask(scratch, _width, _height, x, y,
                                           _allowNegativeComponents, 0, _height,
                                           _cleanMask, _cleanBorderRatio);
  }
}
