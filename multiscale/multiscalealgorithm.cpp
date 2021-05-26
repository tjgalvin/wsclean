#include "multiscalealgorithm.h"

#include "multiscaletransforms.h"

#include "../deconvolution/componentlist.h"
#include "../deconvolution/peakfinder.h"
#include "../deconvolution/subminorloop.h"

#include "../system/fftwmanager.h"

#include "../io/logger.h"

#include "../units/fluxdensity.h"

MultiScaleAlgorithm::MultiScaleAlgorithm(FFTWManager& fftwManager,
                                         double beamSize, double pixelScaleX,
                                         double pixelScaleY)
    : _fftwManager(fftwManager),
      _width(0),
      _height(0),
      _convolutionPadding(1.1),
      _beamSizeInPixels(beamSize / std::max(pixelScaleX, pixelScaleY)),
      _multiscaleScaleBias(0.6),
      _multiscaleGain(0.2),
      _scaleShape(MultiScaleTransforms::TaperedQuadraticShape),
      _maxScales(0),
      _trackPerScaleMasks(false),
      _usePerScaleMasks(false),
      _fastSubMinorLoop(true),
      _trackComponents(false) {
  if (_beamSizeInPixels <= 0.0) _beamSizeInPixels = 1;
}

MultiScaleAlgorithm::~MultiScaleAlgorithm() {
  Logger::Info << "Multi-scale cleaning summary:\n";
  size_t sumComponents = 0;
  float sumFlux = 0.0;
  for (size_t scaleIndex = 0; scaleIndex != _scaleInfos.size(); ++scaleIndex) {
    const ScaleInfo& scaleEntry = _scaleInfos[scaleIndex];
    Logger::Info << "- Scale " << round(scaleEntry.scale)
                 << " px, nr of components cleaned: "
                 << scaleEntry.nComponentsCleaned << " ("
                 << FluxDensity::ToNiceString(scaleEntry.totalFluxCleaned)
                 << ")\n";
    sumComponents += scaleEntry.nComponentsCleaned;
    sumFlux += scaleEntry.totalFluxCleaned;
  }
  Logger::Info << "Total: " << sumComponents << " components ("
               << FluxDensity::ToNiceString(sumFlux) << ")\n";
}

float MultiScaleAlgorithm::ExecuteMajorIteration(
    ImageSet& dirtySet, ImageSet& modelSet,
    const aocommon::UVector<const float*>& psfs, size_t width, size_t height,
    bool& reachedMajorThreshold) {
  // Rough overview of the procedure:
  // Convolve integrated image (all scales)
  // Find integrated peak & scale
  // Minor loop:
  // - Convolve individual images at fixed scale
  // - Subminor loop:
  //   - Measure individual peaks per individually convolved image
  //   - Subtract convolved PSF from individual images
  //   - Subtract twice convolved PSF from individually convolved images
  //   - Find integrated peak at fixed scale
  // - Convolve integrated image (all scales)
  // - Find integrated peak & scale
  //
  // (This excludes creating the convolved PSFs and twice-convolved PSFs
  //  at the appropriate moments).

  if (_stopOnNegativeComponent) _allowNegativeComponents = true;
  _width = width;
  _height = height;
  // The threads always need to be stopped at the end of this function, so we
  // use a scoped unique ptr.
  std::unique_ptr<ThreadedDeconvolutionTools> tools(
      new ThreadedDeconvolutionTools(_threadCount));

  initializeScaleInfo();

  if (_trackPerScaleMasks) {
    // Note that in a second round the nr of scales can be different (due to
    // different width/height, e.g. caused by a different subdivision in
    // parallel cleaning).
    for (const aocommon::UVector<bool>& mask : _scaleMasks) {
      if (mask.size() != _width * _height)
        throw std::runtime_error(
            "Invalid automask size in multiscale algorithm");
    }
    while (_scaleMasks.size() < _scaleInfos.size()) {
      _scaleMasks.emplace_back(_width * _height, false);
    }
  }
  if (_trackComponents) {
    if (_componentList == nullptr)
      _componentList.reset(new ComponentList(
          _width, _height, _scaleInfos.size(), dirtySet.size()));
    else if (_componentList->Width() != _width ||
             _componentList->Height() != _height) {
      throw std::runtime_error("Error in component list dimensions!");
    }
  }
  if (!_rmsFactorImage.empty() && (_rmsFactorImage.Width() != _width ||
                                   _rmsFactorImage.Height() != _height))
    throw std::runtime_error("Error in RMS factor image dimensions!");

  bool hasHitThresholdInSubLoop = false;
  size_t thresholdCountdown = std::max(size_t(8), _scaleInfos.size() * 3 / 2);

  Image scratch, scratchB, integratedScratch;
  // scratch and scratchB are used by the subminorloop, which convolves the
  // images and requires therefore more space. This space depends on the scale,
  // so here the required size for the largest scale is calculated.
  size_t scratchWidth, scratchHeight;
  getConvolutionDimensions(_scaleInfos.size() - 1, scratchWidth, scratchHeight);
  scratch = Image(scratchWidth, scratchHeight);
  scratchB = Image(scratchWidth, scratchHeight);
  integratedScratch = Image(_width, _height);
  std::unique_ptr<std::unique_ptr<Image[]>[]> convolvedPSFs(
      new std::unique_ptr<Image[]>[dirtySet.PSFCount()]);
  dirtySet.GetIntegratedPSF(integratedScratch, psfs);
  convolvePSFs(convolvedPSFs[0], integratedScratch.data(), scratch, true);

  // If there's only one, the integrated equals the first, so we can skip this
  if (dirtySet.PSFCount() > 1) {
    for (size_t i = 0; i != dirtySet.PSFCount(); ++i) {
      convolvePSFs(convolvedPSFs[i], psfs[i], scratch, false);
    }
  }

  MultiScaleTransforms msTransforms(_fftwManager, _width, _height, _scaleShape);

  size_t scaleWithPeak;
  findActiveScaleConvolvedMaxima(dirtySet, integratedScratch, scratch.data(),
                                 true, tools.get());
  if (!selectMaximumScale(scaleWithPeak)) {
    _logReceiver->Warn << "No peak found during multi-scale cleaning! Aborting "
                          "deconvolution.\n";
    reachedMajorThreshold = false;
    return 0.0;
  }

  bool isFinalThreshold = false;
  float mGainThreshold =
      std::fabs(_scaleInfos[scaleWithPeak].maxUnnormalizedImageValue *
                _scaleInfos[scaleWithPeak].biasFactor) *
      (1.0 - _mGain);
  mGainThreshold = std::max(mGainThreshold, MajorIterThreshold());
  float firstThreshold = mGainThreshold;
  if (_threshold > firstThreshold) {
    firstThreshold = _threshold;
    isFinalThreshold = true;
  }

  _logReceiver->Info
      << "Starting multi-scale cleaning. Start peak="
      << FluxDensity::ToNiceString(
             _scaleInfos[scaleWithPeak].maxUnnormalizedImageValue *
             _scaleInfos[scaleWithPeak].biasFactor)
      << ", major iteration threshold="
      << FluxDensity::ToNiceString(firstThreshold);
  if (isFinalThreshold) _logReceiver->Info << " (final)";
  _logReceiver->Info << '\n';

  std::unique_ptr<Image[]> twiceConvolvedPSFs(new Image[dirtySet.PSFCount()]);
  for (size_t i = 0; i != dirtySet.PSFCount(); ++i) {
    twiceConvolvedPSFs[i] = Image(_width, _height);
  }

  ImageSet individualConvolvedImages(&dirtySet.Table(), dirtySet.Settings(),
                                     _width, _height);

  //
  // The minor iteration loop
  //
  while (_iterationNumber < MaxNIter() &&
         std::fabs(_scaleInfos[scaleWithPeak].maxUnnormalizedImageValue *
                   _scaleInfos[scaleWithPeak].biasFactor) > firstThreshold &&
         (!StopOnNegativeComponents() ||
          _scaleInfos[scaleWithPeak].maxUnnormalizedImageValue >= 0.0) &&
         thresholdCountdown > 0) {
    // Create double-convolved PSFs & individually convolved images for this
    // scale
    std::vector<Image> transformList;
    for (size_t i = 0; i != dirtySet.PSFCount(); ++i) {
      float* psf = getConvolvedPSF(i, scaleWithPeak, convolvedPSFs);
      std::copy_n(psf, _width * _height, twiceConvolvedPSFs[i].data());
      transformList.emplace_back(std::move(twiceConvolvedPSFs[i]));
    }
    for (size_t i = 0; i != dirtySet.size(); ++i) {
      std::copy_n(dirtySet[i], _width * _height, individualConvolvedImages[i]);
      transformList.emplace_back(individualConvolvedImages.Release(i));
    }
    if (_scaleInfos[scaleWithPeak].scale != 0.0) {
      tools->MultiScaleTransform(&msTransforms, transformList, scratch,
                                 _scaleInfos[scaleWithPeak].scale);
    }
    for (size_t i = 0; i != dirtySet.PSFCount(); ++i)
      twiceConvolvedPSFs[i] = std::move(transformList[i]);
    for (size_t i = 0; i != dirtySet.size(); ++i)
      individualConvolvedImages.SetImage(
          i, std::move(transformList[i + dirtySet.PSFCount()]));

    //
    // The sub-minor iteration loop for this scale
    //
    float subIterationGainThreshold =
        std::fabs(_scaleInfos[scaleWithPeak].maxUnnormalizedImageValue *
                  _scaleInfos[scaleWithPeak].biasFactor) *
        (1.0 - _multiscaleGain);
    float firstSubIterationThreshold = subIterationGainThreshold;
    if (firstThreshold > firstSubIterationThreshold) {
      firstSubIterationThreshold = firstThreshold;
      if (!hasHitThresholdInSubLoop) {
        _logReceiver->Info << "Subminor loop is near minor loop threshold. "
                              "Initiating countdown.\n";
        hasHitThresholdInSubLoop = true;
      }
      thresholdCountdown--;
      _logReceiver->Info << '(' << thresholdCountdown << ") ";
    }
    // TODO we could chose to run the non-fast loop until we hit e.g. 10
    // iterations in a scale, because the fast loop takes more constant time and
    // is only efficient when doing many iterations.
    if (_fastSubMinorLoop) {
      FFTWManager::ThreadingScope fftwMT(_fftwManager);
      size_t subMinorStartIteration = _iterationNumber;
      size_t convolutionWidth, convolutionHeight;
      getConvolutionDimensions(scaleWithPeak, convolutionWidth,
                               convolutionHeight);
      SubMinorLoop subLoop(_width, _height, convolutionWidth, convolutionHeight,
                           *_logReceiver);
      subLoop.SetIterationInfo(_iterationNumber, MaxNIter());
      subLoop.SetThreshold(
          firstSubIterationThreshold / _scaleInfos[scaleWithPeak].biasFactor,
          subIterationGainThreshold / _scaleInfos[scaleWithPeak].biasFactor);
      subLoop.SetGain(_scaleInfos[scaleWithPeak].gain);
      subLoop.SetAllowNegativeComponents(AllowNegativeComponents());
      subLoop.SetStopOnNegativeComponent(StopOnNegativeComponents());
      const size_t scaleBorder =
                       size_t(ceil(_scaleInfos[scaleWithPeak].scale * 0.5)),
                   horBorderSize = std::max<size_t>(
                       round(width * _cleanBorderRatio), scaleBorder),
                   vertBorderSize = std::max<size_t>(
                       round(height * _cleanBorderRatio), scaleBorder);
      subLoop.SetCleanBorders(horBorderSize, vertBorderSize);
      if (!_rmsFactorImage.empty()) subLoop.SetRMSFactorImage(_rmsFactorImage);
      if (_usePerScaleMasks)
        subLoop.SetMask(_scaleMasks[scaleWithPeak].data());
      else if (_cleanMask)
        subLoop.SetMask(_cleanMask);
      subLoop.SetSpectralFitter(&Fitter());

      aocommon::UVector<const float*> subPSFs(dirtySet.PSFCount());
      for (size_t psfIndex = 0; psfIndex != subPSFs.size(); ++psfIndex)
        subPSFs[psfIndex] = twiceConvolvedPSFs[psfIndex].data();

      subLoop.Run(individualConvolvedImages, subPSFs);

      _iterationNumber = subLoop.CurrentIteration();
      _scaleInfos[scaleWithPeak].nComponentsCleaned +=
          (_iterationNumber - subMinorStartIteration);
      _scaleInfos[scaleWithPeak].totalFluxCleaned += subLoop.FluxCleaned();

      for (size_t imageIndex = 0; imageIndex != dirtySet.size(); ++imageIndex) {
        // TODO this can be multi-threaded if each thread has its own
        // temporaries
        float* psf = getConvolvedPSF(dirtySet.PSFIndex(imageIndex),
                                     scaleWithPeak, convolvedPSFs);
        subLoop.CorrectResidualDirty(_fftwManager, scratch.data(),
                                     scratchB.data(), integratedScratch.data(),
                                     imageIndex, dirtySet[imageIndex], psf);

        subLoop.GetFullIndividualModel(imageIndex, scratch.data());
        if (imageIndex == 0) {
          if (_trackPerScaleMasks)
            subLoop.UpdateAutoMask(_scaleMasks[scaleWithPeak].data());
          if (_trackComponents)
            subLoop.UpdateComponentList(*_componentList, scaleWithPeak);
        }
        if (_scaleInfos[scaleWithPeak].scale != 0.0) {
          std::vector<Image> transformList{std::move(scratch)};
          tools->MultiScaleTransform(&msTransforms, transformList,
                                     integratedScratch,
                                     _scaleInfos[scaleWithPeak].scale);
          scratch = std::move(transformList[0]);
        }
        float* model = modelSet[imageIndex];
        for (size_t i = 0; i != _width * _height; ++i)
          model[i] += scratch.data()[i];
      }

    } else {  // don't use the Clark optimization
      const ScaleInfo& maxScaleInfo = _scaleInfos[scaleWithPeak];
      while (_iterationNumber < MaxNIter() &&
             std::fabs(maxScaleInfo.maxUnnormalizedImageValue *
                       maxScaleInfo.biasFactor) > firstSubIterationThreshold &&
             (!StopOnNegativeComponents() ||
              _scaleInfos[scaleWithPeak].maxUnnormalizedImageValue >= 0.0)) {
        aocommon::UVector<float> componentValues;
        measureComponentValues(componentValues, scaleWithPeak,
                               individualConvolvedImages);
        const size_t x = maxScaleInfo.maxImageValueX;
        const size_t y = maxScaleInfo.maxImageValueY;
        PerformSpectralFit(componentValues.data(), x, y);

        for (size_t imgIndex = 0; imgIndex != dirtySet.size(); ++imgIndex) {
          // Subtract component from individual, non-deconvolved images
          componentValues[imgIndex] =
              componentValues[imgIndex] * maxScaleInfo.gain;

          float* psf = getConvolvedPSF(dirtySet.PSFIndex(imgIndex),
                                       scaleWithPeak, convolvedPSFs);
          tools->SubtractImage(dirtySet[imgIndex], psf, _width, _height, x, y,
                               componentValues[imgIndex]);

          // Subtract twice convolved PSFs from convolved images
          tools->SubtractImage(
              individualConvolvedImages[imgIndex],
              twiceConvolvedPSFs[dirtySet.PSFIndex(imgIndex)].data(), _width,
              _height, x, y, componentValues[imgIndex]);
          // TODO this is incorrect, but why is the residual without
          // Cotton-Schwab still OK ? Should test
          // tools->SubtractImage(individualConvolvedImages[imgIndex], psf,
          // _width, _height, x, y, componentValues[imgIndex]);

          // Adjust model
          addComponentToModel(modelSet[imgIndex], scaleWithPeak,
                              componentValues[imgIndex]);
        }
        if (_trackComponents) {
          _componentList->Add(x, y, scaleWithPeak, componentValues.data());
        }

        // Find maximum for this scale
        individualConvolvedImages.GetLinearIntegrated(integratedScratch);
        findPeakDirect(integratedScratch.data(), scratch.data(), scaleWithPeak);
        _logReceiver->Debug
            << "Scale now "
            << std::fabs(_scaleInfos[scaleWithPeak].maxUnnormalizedImageValue *
                         _scaleInfos[scaleWithPeak].biasFactor)
            << '\n';

        ++_iterationNumber;
      }
    }

    activateScales(scaleWithPeak);

    findActiveScaleConvolvedMaxima(dirtySet, integratedScratch, scratch.data(),
                                   false, tools.get());

    if (!selectMaximumScale(scaleWithPeak)) {
      _logReceiver->Warn << "No peak found in main loop of multi-scale "
                            "cleaning! Aborting deconvolution.\n";
      reachedMajorThreshold = false;
      return 0.0;
    }

    _logReceiver->Info
        << "Iteration " << _iterationNumber << ", scale "
        << round(_scaleInfos[scaleWithPeak].scale) << " px : "
        << FluxDensity::ToNiceString(
               _scaleInfos[scaleWithPeak].maxUnnormalizedImageValue *
               _scaleInfos[scaleWithPeak].biasFactor)
        << " at " << _scaleInfos[scaleWithPeak].maxImageValueX << ','
        << _scaleInfos[scaleWithPeak].maxImageValueY << '\n';
  }

  bool maxIterReached = _iterationNumber >= MaxNIter(),
       negativeReached =
           StopOnNegativeComponents() &&
           _scaleInfos[scaleWithPeak].maxUnnormalizedImageValue < 0.0;
  // finalThresholdReached =
  // std::fabs(_scaleInfos[scaleWithPeak].maxUnnormalizedImageValue *
  // _scaleInfos[scaleWithPeak].biasFactor) <= _threshold;

  if (maxIterReached)
    _logReceiver->Info << "Cleaning finished because maximum number of "
                          "iterations was reached.\n";
  else if (negativeReached)
    _logReceiver->Info
        << "Cleaning finished because a negative component was found.\n";
  else if (isFinalThreshold)
    _logReceiver->Info
        << "Cleaning finished because the final threshold was reached.\n";
  else
    _logReceiver->Info << "Minor loop finished, continuing cleaning after "
                          "inversion/prediction round.\n";

  reachedMajorThreshold =
      !maxIterReached && !isFinalThreshold && !negativeReached;
  return _scaleInfos[scaleWithPeak].maxUnnormalizedImageValue *
         _scaleInfos[scaleWithPeak].biasFactor;
}

void MultiScaleAlgorithm::initializeScaleInfo() {
  if (_manualScaleList.empty()) {
    if (_scaleInfos.empty()) {
      size_t scaleIndex = 0;
      double scale = _beamSizeInPixels * 2.0;
      do {
        _scaleInfos.push_back(ScaleInfo());
        ScaleInfo& newEntry = _scaleInfos.back();
        if (scaleIndex == 0)
          newEntry.scale = 0.0;
        else
          newEntry.scale = scale;
        newEntry.kernelPeak = MultiScaleTransforms::KernelPeakValue(
            scale, std::min(_width, _height), _scaleShape);

        scale *= 2.0;
        ++scaleIndex;
      } while (scale < std::min(_width, _height) * 0.5 &&
               (_maxScales == 0 || scaleIndex < _maxScales));
    } else {
      while (!_scaleInfos.empty() &&
             _scaleInfos.back().scale >= std::min(_width, _height) * 0.5) {
        _logReceiver->Info
            << "Scale size " << _scaleInfos.back().scale
            << " does not fit in cleaning region: removing scale.\n";
        _scaleInfos.erase(_scaleInfos.begin() + _scaleInfos.size() - 1);
      }
    }
  }
  if (!_manualScaleList.empty() && _scaleInfos.empty()) {
    std::sort(_manualScaleList.begin(), _manualScaleList.end());
    for (size_t scaleIndex = 0; scaleIndex != _manualScaleList.size();
         ++scaleIndex) {
      _scaleInfos.push_back(ScaleInfo());
      ScaleInfo& newEntry = _scaleInfos.back();
      newEntry.scale = _manualScaleList[scaleIndex];
      newEntry.kernelPeak = MultiScaleTransforms::KernelPeakValue(
          newEntry.scale, std::min(_width, _height), _scaleShape);
    }
  }
}

void MultiScaleAlgorithm::convolvePSFs(std::unique_ptr<Image[]>& convolvedPSFs,
                                       const float* psf, Image& scratch,
                                       bool isIntegrated) {
  MultiScaleTransforms msTransforms(_fftwManager, _width, _height, _scaleShape);
  convolvedPSFs.reset(new Image[_scaleInfos.size()]);
  if (isIntegrated) _logReceiver->Info << "Scale info:\n";
  const double firstAutoScaleSize = _beamSizeInPixels * 2.0;
  for (size_t scaleIndex = 0; scaleIndex != _scaleInfos.size(); ++scaleIndex) {
    ScaleInfo& scaleEntry = _scaleInfos[scaleIndex];

    convolvedPSFs[scaleIndex] = Image(_width, _height);
    std::copy_n(psf, _width * _height, convolvedPSFs[scaleIndex].data());

    if (isIntegrated) {
      if (scaleEntry.scale != 0.0)
        msTransforms.Transform(convolvedPSFs[scaleIndex], scratch,
                               scaleEntry.scale);

      scaleEntry.psfPeak =
          convolvedPSFs[scaleIndex][_width / 2 + (_height / 2) * _width];
      // We normalize this factor to 1 for scale 0, so:
      // factor = (psf / kernel) / (psf0 / kernel0) = psf * kernel0 / (kernel *
      // psf0)
      // scaleEntry.biasFactor = std::max(1.0,
      //	scaleEntry.psfPeak * scaleInfos[0].kernelPeak /
      //	(scaleEntry.kernelPeak * scaleInfos[0].psfPeak));
      double expTerm;
      if (scaleEntry.scale == 0.0 || _scaleInfos.size() < 2)
        expTerm = 0.0;
      else
        expTerm = std::log2(scaleEntry.scale / firstAutoScaleSize);
      scaleEntry.biasFactor =
          std::pow(_multiscaleScaleBias, -double(expTerm)) * 1.0;

      // I tried this, but wasn't perfect:
      // _gain * _scaleInfos[0].kernelPeak / scaleEntry.kernelPeak;
      scaleEntry.gain = _gain / scaleEntry.psfPeak;

      scaleEntry.isActive = true;

      if (scaleEntry.scale == 0.0)
        std::copy_n(psf, _width * _height, convolvedPSFs[scaleIndex].data());

      _logReceiver->Info << "- Scale " << round(scaleEntry.scale)
                         << ", bias factor="
                         << round(scaleEntry.biasFactor * 10.0) / 10.0
                         << ", psfpeak=" << scaleEntry.psfPeak
                         << ", gain=" << scaleEntry.gain
                         << ", kernel peak=" << scaleEntry.kernelPeak << '\n';
    } else {
      if (scaleEntry.scale != 0.0)
        msTransforms.Transform(convolvedPSFs[scaleIndex], scratch,
                               scaleEntry.scale);
    }
  }
}

void MultiScaleAlgorithm::findActiveScaleConvolvedMaxima(
    const ImageSet& imageSet, Image& integratedScratch, float* scratch,
    bool reportRMS, ThreadedDeconvolutionTools* tools) {
  MultiScaleTransforms msTransforms(_fftwManager, _width, _height, _scaleShape);
  // ImageBufferAllocator::Ptr convolvedImage;
  //_allocator.Allocate(_width*_height, convolvedImage);
  imageSet.GetLinearIntegrated(integratedScratch);
  aocommon::UVector<float> transformScales;
  aocommon::UVector<size_t> transformIndices;
  std::vector<aocommon::UVector<bool>> transformScaleMasks;
  for (size_t scaleIndex = 0; scaleIndex != _scaleInfos.size(); ++scaleIndex) {
    ScaleInfo& scaleEntry = _scaleInfos[scaleIndex];
    if (scaleEntry.isActive) {
      if (scaleEntry.scale == 0) {
        // Don't convolve scale 0: this is the delta function scale
        findPeakDirect(integratedScratch.data(), scratch, scaleIndex);
        if (reportRMS)
          scaleEntry.rms = ThreadedDeconvolutionTools::RMS(integratedScratch,
                                                           _width * _height);
      } else {
        transformScales.push_back(scaleEntry.scale);
        transformIndices.push_back(scaleIndex);
        if (_usePerScaleMasks)
          transformScaleMasks.push_back(_scaleMasks[scaleIndex]);
      }
    }
  }
  std::vector<ThreadedDeconvolutionTools::PeakData> results;

  tools->FindMultiScalePeak(&msTransforms, integratedScratch, transformScales,
                            results, _allowNegativeComponents, _cleanMask,
                            transformScaleMasks, _cleanBorderRatio,
                            _rmsFactorImage, reportRMS);

  for (size_t i = 0; i != results.size(); ++i) {
    ScaleInfo& scaleEntry = _scaleInfos[transformIndices[i]];
    scaleEntry.maxNormalizedImageValue =
        get_optional_value_or(results[i].normalizedValue, 0.0);
    scaleEntry.maxUnnormalizedImageValue =
        get_optional_value_or(results[i].unnormalizedValue, 0.0);
    scaleEntry.maxImageValueX = results[i].x;
    scaleEntry.maxImageValueY = results[i].y;
    if (reportRMS) scaleEntry.rms = results[i].rms;
  }
  if (reportRMS) {
    _logReceiver->Info << "RMS per scale: {";
    for (size_t scaleIndex = 0; scaleIndex != _scaleInfos.size();
         ++scaleIndex) {
      ScaleInfo& scaleEntry = _scaleInfos[scaleIndex];
      if (scaleIndex != 0) _logReceiver->Info << ", ";
      _logReceiver->Info << round(scaleEntry.scale) << ": "
                         << FluxDensity::ToNiceString(scaleEntry.rms);
    }
    _logReceiver->Info << "}\n";
  }
}

bool MultiScaleAlgorithm::selectMaximumScale(size_t& scaleWithPeak) {
  // Find max component
  std::map<float, size_t> peakToScaleMap;
  for (size_t i = 0; i != _scaleInfos.size(); ++i) {
    if (_scaleInfos[i].isActive) {
      float maxVal = std::fabs(_scaleInfos[i].maxUnnormalizedImageValue *
                               _scaleInfos[i].biasFactor);
      peakToScaleMap.insert(std::make_pair(maxVal, i));
    }
  }
  if (peakToScaleMap.empty()) {
    scaleWithPeak = size_t(-1);
    return false;
  } else {
    std::map<float, size_t>::const_reverse_iterator mapIter =
        peakToScaleMap.rbegin();
    scaleWithPeak = mapIter->second;
    return true;
  }
}

void MultiScaleAlgorithm::activateScales(size_t scaleWithLastPeak) {
  for (size_t i = 0; i != _scaleInfos.size(); ++i) {
    bool doActivate =
        i == scaleWithLastPeak ||
        /*i == runnerUp ||*/ std::fabs(
            _scaleInfos[i].maxUnnormalizedImageValue) *
                _scaleInfos[i].biasFactor >
            std::fabs(
                _scaleInfos[scaleWithLastPeak].maxUnnormalizedImageValue) *
                (1.0 - _gain) * _scaleInfos[scaleWithLastPeak].biasFactor;
    if (!_scaleInfos[i].isActive && doActivate) {
      _logReceiver->Debug << "Scale " << _scaleInfos[i].scale
                          << " is now significant and is activated.\n";
      _scaleInfos[i].isActive = true;
    } else if (_scaleInfos[i].isActive && !doActivate) {
      _logReceiver->Debug << "Scale " << _scaleInfos[i].scale
                          << " is insignificant and is deactivated.\n";
      _scaleInfos[i].isActive = false;
    }
  }
}

void MultiScaleAlgorithm::measureComponentValues(
    aocommon::UVector<float>& componentValues, size_t scaleIndex,
    ImageSet& imageSet) {
  const ScaleInfo& scale = _scaleInfos[scaleIndex];
  componentValues.resize(imageSet.size());
  _logReceiver->Debug << "Measuring " << scale.maxImageValueX << ','
                      << scale.maxImageValueY << ", scale " << scale.scale
                      << ", integrated=" << scale.maxUnnormalizedImageValue
                      << ":";
  for (size_t i = 0; i != imageSet.size(); ++i) {
    componentValues[i] =
        imageSet[i][scale.maxImageValueX + scale.maxImageValueY * _width];
    _logReceiver->Debug << ' ' << componentValues[i];
  }
  _logReceiver->Debug << '\n';
}

void MultiScaleAlgorithm::addComponentToModel(float* model,
                                              size_t scaleWithPeak,
                                              float componentValue) {
  const size_t x = _scaleInfos[scaleWithPeak].maxImageValueX,
               y = _scaleInfos[scaleWithPeak].maxImageValueY;
  if (_scaleInfos[scaleWithPeak].scale == 0.0)
    model[x + _width * y] += componentValue;
  else
    MultiScaleTransforms::AddShapeComponent(model, _width, _height,
                                            _scaleInfos[scaleWithPeak].scale, x,
                                            y, componentValue, _scaleShape);

  _scaleInfos[scaleWithPeak].nComponentsCleaned++;
  _scaleInfos[scaleWithPeak].totalFluxCleaned += componentValue;

  if (_trackPerScaleMasks) {
    _scaleMasks[scaleWithPeak][x + _width * y] = true;
  }
}

float* MultiScaleAlgorithm::getConvolvedPSF(
    size_t psfIndex, size_t scaleIndex,
    const std::unique_ptr<std::unique_ptr<Image[]>[]>& convolvedPSFs) {
  return convolvedPSFs[psfIndex][scaleIndex].data();
}

void MultiScaleAlgorithm::findPeakDirect(const float* image, float* scratch,
                                         size_t scaleIndex) {
  ScaleInfo& scaleInfo = _scaleInfos[scaleIndex];
  const size_t horBorderSize = round(_width * _cleanBorderRatio),
               vertBorderSize = round(_height * _cleanBorderRatio);
  const float* actualImage;
  if (_rmsFactorImage.empty())
    actualImage = image;
  else {
    for (size_t i = 0; i != _rmsFactorImage.size(); ++i)
      scratch[i] = image[i] * _rmsFactorImage[i];
    actualImage = scratch;
  }

  boost::optional<float> maxValue;
  if (_usePerScaleMasks)
    maxValue = PeakFinder::FindWithMask(
        actualImage, _width, _height, scaleInfo.maxImageValueX,
        scaleInfo.maxImageValueY, _allowNegativeComponents, 0, _height,
        _scaleMasks[scaleIndex].data(), horBorderSize, vertBorderSize);
  else if (_cleanMask == 0)
    maxValue =
        PeakFinder::Find(actualImage, _width, _height, scaleInfo.maxImageValueX,
                         scaleInfo.maxImageValueY, _allowNegativeComponents, 0,
                         _height, horBorderSize, vertBorderSize);
  else
    maxValue = PeakFinder::FindWithMask(
        actualImage, _width, _height, scaleInfo.maxImageValueX,
        scaleInfo.maxImageValueY, _allowNegativeComponents, 0, _height,
        _cleanMask, horBorderSize, vertBorderSize);

  scaleInfo.maxUnnormalizedImageValue = get_optional_value_or(maxValue, 0.0);
  if (_rmsFactorImage.empty())
    scaleInfo.maxNormalizedImageValue = get_optional_value_or(maxValue, 0.0);
  else
    scaleInfo.maxNormalizedImageValue =
        get_optional_value_or(maxValue, 0.0) /
        _rmsFactorImage[scaleInfo.maxImageValueX +
                        scaleInfo.maxImageValueY * _width];
}

static size_t calculateGoodFFTSize(size_t n) {
  size_t bestfac = 2 * n;
  /* NOTE: Starting from f2=2 here instead from f2=1 as usual, because the
                  result needs to be even. */
  for (size_t f2 = 2; f2 < bestfac; f2 *= 2)
    for (size_t f23 = f2; f23 < bestfac; f23 *= 3)
      for (size_t f235 = f23; f235 < bestfac; f235 *= 5)
        for (size_t f2357 = f235; f2357 < bestfac; f2357 *= 7)
          if (f2357 >= n) bestfac = f2357;
  return bestfac;
}

void MultiScaleAlgorithm::getConvolutionDimensions(size_t scaleIndex,
                                                   size_t& width,
                                                   size_t& height) const {
  double scale = _scaleInfos[scaleIndex].scale;
  // The factor of 1.5 comes from some superficial experience with diverging
  // runs. It's supposed to be a balance between diverging runs caused by
  // insufficient padding on one hand, and taking up too much memory on the
  // other. I've seen divergence when padding=1.1, _width=1500, max scale=726
  // and conv width=1650. Divergence occurred on scale 363. Was solved with conv
  // width=2250. 2250 = 1.1*(363*factor + 1500)  --> factor = 1.5 And solved
  // with conv width=2000. 2000 = 1.1*(363*factor + 1500)  --> factor = 0.8
  width = ceil(_convolutionPadding * (scale * 1.5 + _width));
  height = ceil(_convolutionPadding * (scale * 1.5 + _height));
  width = calculateGoodFFTSize(width);
  height = calculateGoodFFTSize(height);
}
