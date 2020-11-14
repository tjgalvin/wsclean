#include "deconvolutionalgorithm.h"

#include "../system/system.h"

#include "../model/modelsource.h"
#include "../model/model.h"
#include "../model/powerlawsed.h"

#include <aocommon/imagecoordinates.h>

using namespace aocommon;

DeconvolutionAlgorithm::DeconvolutionAlgorithm()
    : _threshold(0.0),
      _majorIterThreshold(0.0),
      _gain(0.1),
      _mGain(1.0),
      _cleanBorderRatio(0.05),
      _maxIter(500),
      _iterationNumber(0),
      _threadCount(System::ProcessorCount()),
      _allowNegativeComponents(true),
      _stopOnNegativeComponent(false),
      _cleanMask(0),
      _logReceiver(nullptr),
      _spectralFitter(NoSpectralFitting, 0) {}

void DeconvolutionAlgorithm::GetModelFromImage(
    Model& model, const float* image, size_t width, size_t height,
    double phaseCentreRA, double phaseCentreDec, double pixelSizeX,
    double pixelSizeY, double phaseCentreDL, double phaseCentreDM,
    double spectralIndex, double refFreq,
    aocommon::PolarizationEnum polarization) {
  for (size_t y = 0; y != height; ++y) {
    for (size_t x = 0; x != width; ++x) {
      double value = image[y * width + x];
      if (value != 0.0 && std::isfinite(value)) {
        long double l, m;
        ImageCoordinates::XYToLM<long double>(x, y, pixelSizeX, pixelSizeY,
                                              width, height, l, m);
        l += phaseCentreDL;
        m += phaseCentreDM;
        ModelComponent component;
        long double ra, dec;
        ImageCoordinates::LMToRaDec<long double>(l, m, phaseCentreRA,
                                                 phaseCentreDec, ra, dec);
        std::stringstream nameStr;
        nameStr << "component" << model.SourceCount();
        component.SetSED(
            MeasuredSED(value, refFreq, spectralIndex, polarization));
        component.SetPosRA(ra);
        component.SetPosDec(dec);

        ModelSource source;
        source.SetName(nameStr.str());
        source.AddComponent(component);
        model.AddSource(source);
      }
    }
  }
}

void DeconvolutionAlgorithm::ResizeImage(float* dest, size_t newWidth,
                                         size_t newHeight, const float* source,
                                         size_t width, size_t height) {
  size_t srcStartX = (width - newWidth) / 2,
         srcStartY = (height - newHeight) / 2;
  for (size_t y = 0; y != newHeight; ++y) {
    float* destPtr = dest + y * newWidth;
    const float* srcPtr = source + (y + srcStartY) * width + srcStartX;
    memcpy(destPtr, srcPtr, newWidth * sizeof(double));
  }
}

void DeconvolutionAlgorithm::RemoveNaNsInPSF(float* psf, size_t width,
                                             size_t height) {
  float* endPtr = psf + width * height;
  while (psf != endPtr) {
    if (!std::isfinite(*psf)) *psf = 0.0;
    ++psf;
  }
}

void DeconvolutionAlgorithm::PerformSpectralFit(float* values) {
  _spectralFitter.FitAndEvaluate(values);
}

float Evaluate(float x, const aocommon::UVector<float>& terms,
               double referenceFrequencyHz = 1.0);
