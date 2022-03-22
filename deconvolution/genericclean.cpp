#include "genericclean.h"

#include "subminorloop.h"
#include "peakfinder.h"

#include "../multiscale/threadeddeconvolutiontools.h"

#include <aocommon/image.h>
#include <aocommon/lane.h>
#include <aocommon/units/fluxdensity.h>

using aocommon::units::FluxDensity;
namespace {
std::string peakDescription(const aocommon::Image& image, size_t x, size_t y) {
  std::ostringstream str;
  const size_t index = x + y * image.Width();
  const float peak = image[index];
  str << FluxDensity::ToNiceString(peak) << " at " << x << "," << y;
  return str.str();
}
}  // namespace

GenericClean::GenericClean(bool useSubMinorOptimization)
    : _convolutionPadding(1.1),
      _useSubMinorOptimization(useSubMinorOptimization) {}

float GenericClean::ExecuteMajorIteration(
    ImageSet& dirtySet, ImageSet& modelSet,
    const std::vector<aocommon::Image>& psfs, bool& reachedMajorThreshold) {
  const size_t width = dirtySet.Width();
  const size_t height = dirtySet.Height();
  const size_t iterationCounterAtStart = _iterationNumber;
  if (_stopOnNegativeComponent) _allowNegativeComponents = true;
  _convolutionWidth = ceil(_convolutionPadding * width);
  _convolutionHeight = ceil(_convolutionPadding * height);
  if (_convolutionWidth % 2 != 0) ++_convolutionWidth;
  if (_convolutionHeight % 2 != 0) ++_convolutionHeight;

  aocommon::Image integrated(width, height);
  aocommon::Image scratchA(_convolutionWidth, _convolutionHeight);
  aocommon::Image scratchB(_convolutionWidth, _convolutionHeight);
  dirtySet.GetLinearIntegrated(integrated);
  size_t componentX = 0;
  size_t componentY = 0;
  std::optional<float> maxValue =
      findPeak(integrated, scratchA.Data(), componentX, componentY);
  if (!maxValue) {
    _logReceiver->Info << "No peak found.\n";
    reachedMajorThreshold = false;
    return 0.0;
  }
  _logReceiver->Info << "Initial peak: "
                     << peakDescription(integrated, componentX, componentY)
                     << '\n';
  float firstThreshold = this->_threshold;
  float majorIterThreshold = std::max<float>(
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
    SubMinorLoop subMinorLoop(width, height, _convolutionWidth,
                              _convolutionHeight, *_logReceiver);
    subMinorLoop.SetIterationInfo(_iterationNumber, MaxNIter());
    subMinorLoop.SetThreshold(firstThreshold, firstThreshold * 0.99);
    subMinorLoop.SetGain(Gain());
    subMinorLoop.SetAllowNegativeComponents(AllowNegativeComponents());
    subMinorLoop.SetStopOnNegativeComponent(StopOnNegativeComponents());
    subMinorLoop.SetSpectralFitter(&Fitter());
    if (!_rmsFactorImage.Empty())
      subMinorLoop.SetRMSFactorImage(_rmsFactorImage);
    if (_cleanMask) subMinorLoop.SetMask(_cleanMask);
    const size_t horBorderSize = std::round(width * CleanBorderRatio());
    const size_t vertBorderSize = std::round(height * CleanBorderRatio());
    subMinorLoop.SetCleanBorders(horBorderSize, vertBorderSize);
    subMinorLoop.SetThreadCount(_threadCount);

    maxValue = subMinorLoop.Run(dirtySet, psfs);

    _iterationNumber = subMinorLoop.CurrentIteration();

    _logReceiver->Info
        << "Performed " << _iterationNumber << " iterations in total, "
        << (_iterationNumber - startIteration)
        << " in this major iteration with sub-minor optimization.\n";

    for (size_t imageIndex = 0; imageIndex != dirtySet.size(); ++imageIndex) {
      // TODO this can be multi-threaded if each thread has its own temporaries
      const aocommon::Image& psf = psfs[dirtySet.PSFIndex(imageIndex)];
      subMinorLoop.CorrectResidualDirty(scratchA.Data(), scratchB.Data(),
                                        integrated.Data(), imageIndex,
                                        dirtySet.Data(imageIndex), psf.Data());

      subMinorLoop.GetFullIndividualModel(imageIndex, scratchA.Data());
      float* model = modelSet.Data(imageIndex);
      for (size_t i = 0; i != width * height; ++i)
        model[i] += scratchA.Data()[i];
    }
  } else {
    ThreadedDeconvolutionTools tools(_threadCount);
    size_t peakIndex = componentX + componentY * width;

    aocommon::UVector<float> peakValues(dirtySet.size());

    while (maxValue && fabs(*maxValue) > firstThreshold &&
           this->_iterationNumber < this->_maxIter &&
           !(maxValue < 0.0f && this->_stopOnNegativeComponent)) {
      if (this->_iterationNumber <= 10 ||
          (this->_iterationNumber <= 100 && this->_iterationNumber % 10 == 0) ||
          (this->_iterationNumber <= 1000 &&
           this->_iterationNumber % 100 == 0) ||
          this->_iterationNumber % 1000 == 0)
        _logReceiver->Info << "Iteration " << this->_iterationNumber << ": "
                           << peakDescription(integrated, componentX,
                                              componentY)
                           << '\n';

      for (size_t i = 0; i != dirtySet.size(); ++i)
        peakValues[i] = dirtySet[i][peakIndex];

      PerformSpectralFit(peakValues.data(), componentX, componentY);

      for (size_t i = 0; i != dirtySet.size(); ++i) {
        peakValues[i] *= this->_gain;
        modelSet.Data(i)[peakIndex] += peakValues[i];

        size_t psfIndex = dirtySet.PSFIndex(i);

        tools.SubtractImage(dirtySet.Data(i), psfs[psfIndex], componentX,
                            componentY, peakValues[i]);
      }

      dirtySet.GetSquareIntegrated(integrated, scratchA);
      maxValue = findPeak(integrated, scratchA.Data(), componentX, componentY);

      peakIndex = componentX + componentY * width;

      ++this->_iterationNumber;
    }
  }
  if (maxValue) {
    _logReceiver->Info << "Stopped on peak "
                       << FluxDensity::ToNiceString(*maxValue) << ", because ";
    bool maxIterReached = _iterationNumber >= MaxNIter(),
         finalThresholdReached =
             std::fabs(*maxValue) <= _threshold || maxValue == 0.0f,
         negativeReached = maxValue < 0.0f && this->_stopOnNegativeComponent,
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

std::optional<float> GenericClean::findPeak(const aocommon::Image& image,
                                            float* scratch_buffer, size_t& x,
                                            size_t& y) {
  const float* actual_image = image.Data();
  if (!_rmsFactorImage.Empty()) {
    std::copy_n(image.Data(), image.Size(), scratch_buffer);
    for (size_t i = 0; i != image.Size(); ++i) {
      scratch_buffer[i] *= _rmsFactorImage[i];
    }
    actual_image = scratch_buffer;
  }

  if (_cleanMask == nullptr)
    return PeakFinder::Find(actual_image, image.Width(), image.Height(), x, y,
                            _allowNegativeComponents, 0, image.Height(),
                            _cleanBorderRatio);
  else
    return PeakFinder::FindWithMask(actual_image, image.Width(), image.Height(),
                                    x, y, _allowNegativeComponents, 0,
                                    image.Height(), _cleanMask,
                                    _cleanBorderRatio);
}
