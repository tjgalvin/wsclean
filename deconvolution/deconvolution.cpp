#include "deconvolution.h"

#include "imageset.h"
#include "simpleclean.h"
#include "moresane.h"
#include "iuwtdeconvolution.h"
#include "genericclean.h"

#include "../multiscale/multiscalealgorithm.h"

#include "../casamaskreader.h"
#include "../fitsreader.h"
#include "../image.h"
#include "../ndppp.h"
#include "../rmsimage.h"

#include "../units/fluxdensity.h"

#include "../wsclean/imagingtable.h"
#include "../wsclean/wscleansettings.h"

using namespace aocommon;

Deconvolution::Deconvolution(const class WSCleanSettings& settings)
    : _settings(settings),
      _parallelDeconvolution(settings),
      _autoMaskIsFinished(false),
      _imgWidth(0),  // these are not yet set the in settings obj -- load later
      _imgHeight(0),
      _beamSize(0.0),
      _pixelScaleX(0),
      _pixelScaleY(0) {}

Deconvolution::~Deconvolution() { FreeDeconvolutionAlgorithms(); }

void Deconvolution::Perform(const class ImagingTable& groupTable,
                            bool& reachedMajorThreshold,
                            size_t majorIterationNr) {
  Logger::Info.Flush();
  Logger::Info << " == Cleaning (" << majorIterationNr << ") ==\n";

  ImageSet residualSet(&groupTable, _settings, _imgWidth, _imgHeight),
      modelSet(&groupTable, _settings, _imgWidth, _imgHeight);

  residualSet.LoadAndAverage(*_residualImages);
  modelSet.LoadAndAverage(*_modelImages);

  Image integrated(_imgWidth, _imgHeight);
  residualSet.GetLinearIntegrated(integrated.data());
  double stddev = integrated.StdDevFromMAD();
  Logger::Info << "Estimated standard deviation of background noise: "
               << FluxDensity::ToNiceString(stddev) << '\n';
  if (_settings.autoMask && _autoMaskIsFinished) {
    // When we are in the second phase of automasking, don't use
    // the RMS background anymore
    _parallelDeconvolution.SetRMSFactorImage(Image());
  } else {
    if (!_settings.localRMSImage.empty()) {
      Image rmsImage(_imgWidth, _imgHeight);
      FitsReader reader(_settings.localRMSImage);
      reader.Read(rmsImage.data());
      // Normalize the RMS image
      stddev = rmsImage.Min();
      Logger::Info << "Lowest RMS in image: "
                   << FluxDensity::ToNiceString(stddev) << '\n';
      if (stddev <= 0.0)
        throw std::runtime_error(
            "RMS image can only contain values > 0, but contains values <= "
            "0.0");
      for (double& value : rmsImage) {
        if (value != 0.0) value = stddev / value;
      }
      _parallelDeconvolution.SetRMSFactorImage(std::move(rmsImage));
    } else if (_settings.localRMS) {
      Image rmsImage;
      // TODO this should use full beam parameters
      switch (_settings.localRMSMethod) {
        case WSCleanSettings::RMSWindow:
          RMSImage::Make(rmsImage, integrated, _settings.localRMSWindow,
                         _beamSize, _beamSize, 0.0, _pixelScaleX, _pixelScaleY);
          break;
        case WSCleanSettings::RMSAndMinimumWindow:
          RMSImage::MakeWithNegativityLimit(
              rmsImage, integrated, _settings.localRMSWindow, _beamSize,
              _beamSize, 0.0, _pixelScaleX, _pixelScaleY);
          break;
      }
      // Normalize the RMS image relative to the threshold so that Jy remains
      // Jy.
      stddev = rmsImage.Min();
      Logger::Info << "Lowest RMS in image: "
                   << FluxDensity::ToNiceString(stddev) << '\n';
      for (double& value : rmsImage) {
        if (value != 0.0) value = stddev / value;
      }
      _parallelDeconvolution.SetRMSFactorImage(std::move(rmsImage));
    }
  }
  if (_settings.autoMask && !_autoMaskIsFinished)
    _parallelDeconvolution.SetThreshold(std::max(
        stddev * _settings.autoMaskSigma, _settings.deconvolutionThreshold));
  else if (_settings.autoDeconvolutionThreshold)
    _parallelDeconvolution.SetThreshold(
        std::max(stddev * _settings.autoDeconvolutionThresholdSigma,
                 _settings.deconvolutionThreshold));
  integrated.reset();

  std::vector<aocommon::UVector<double>> psfVecs(residualSet.PSFCount());
  residualSet.LoadAndAveragePSFs(*_psfImages, psfVecs, _psfPolarization);

  aocommon::UVector<const double*> psfs(residualSet.PSFCount());
  for (size_t i = 0; i != psfVecs.size(); ++i) psfs[i] = psfVecs[i].data();

  if (_settings.useMultiscale) {
    if (_settings.autoMask) {
      if (_autoMaskIsFinished)
        _parallelDeconvolution.SetAutoMaskMode(false, true);
      else
        _parallelDeconvolution.SetAutoMaskMode(true, false);
    }
  } else {
    if (_settings.autoMask && _autoMaskIsFinished) {
      if (_autoMask.empty()) {
        _autoMask.resize(_imgWidth * _imgHeight);
        for (size_t imgIndex = 0; imgIndex != modelSet.size(); ++imgIndex) {
          const double* image = modelSet[imgIndex];
          for (size_t i = 0; i != _imgWidth * _imgHeight; ++i) {
            _autoMask[i] = (image[i] == 0.0) ? false : true;
          }
        }
      }
      _parallelDeconvolution.SetCleanMask(_autoMask.data());
    }
  }

  _parallelDeconvolution.ExecuteMajorIteration(residualSet, modelSet, psfs,
                                               reachedMajorThreshold);

  if (!reachedMajorThreshold && _settings.autoMask && !_autoMaskIsFinished) {
    Logger::Info << "Auto-masking threshold reached; continuing next major "
                    "iteration with deeper threshold and mask.\n";
    _autoMaskIsFinished = true;
    reachedMajorThreshold = true;
  }

  if (_settings.majorIterationCount != 0 &&
      majorIterationNr >= _settings.majorIterationCount) {
    reachedMajorThreshold = false;
    Logger::Info << "Maximum number of major iterations was reached: not "
                    "continuing deconvolution.\n";
  }

  residualSet.AssignAndStore(*_residualImages);
  modelSet.InterpolateAndStore(
      *_modelImages, _parallelDeconvolution.FirstAlgorithm().Fitter());
}

void Deconvolution::InitializeDeconvolutionAlgorithm(
    const ImagingTable& groupTable, aocommon::PolarizationEnum psfPolarization,
    double beamSize, size_t threadCount) {
  _imgWidth = _settings.trimmedImageWidth;
  _imgHeight = _settings.trimmedImageHeight;
  _pixelScaleX = _settings.pixelScaleX;
  _pixelScaleY = _settings.pixelScaleY;

  _psfPolarization = psfPolarization;
  _beamSize = beamSize;
  _autoMaskIsFinished = false;
  _autoMask.clear();
  FreeDeconvolutionAlgorithms();
  if (groupTable.SquaredGroupCount() == 0)
    throw std::runtime_error("Nothing to clean");

  if (!std::isfinite(_beamSize)) {
    Logger::Warn << "No proper beam size available in deconvolution!\n";
    _beamSize = 0.0;
  }

  ImagingTable firstSquaredGroup = groupTable.GetSquaredGroup(0);
  size_t squaredCount = firstSquaredGroup.EntryCount();
  _polarizations.clear();
  for (size_t p = 0; p != squaredCount; ++p) {
    if (_polarizations.count(firstSquaredGroup[p].polarization) != 0)
      throw std::runtime_error(
          "Two equal polarizations were given to the deconvolution algorithm "
          "within a single polarized group");
    else
      _polarizations.insert(firstSquaredGroup[p].polarization);
  }

  std::unique_ptr<class DeconvolutionAlgorithm> algorithm;

  if (_settings.useMoreSaneDeconvolution) {
    algorithm.reset(
        new MoreSane(_settings.moreSaneLocation, _settings.moreSaneArgs,
                     _settings.moreSaneSigmaLevels, _settings.prefixName,
                     _parallelDeconvolution.GetFFTWManager()));
  } else if (_settings.useIUWTDeconvolution) {
    IUWTDeconvolution* method =
        new IUWTDeconvolution(_parallelDeconvolution.GetFFTWManager());
    method->SetUseSNRTest(_settings.iuwtSNRTest);
    algorithm.reset(method);
  } else if (_settings.useMultiscale) {
    MultiScaleAlgorithm* msAlgorithm =
        new MultiScaleAlgorithm(_parallelDeconvolution.GetFFTWManager(),
                                beamSize, _pixelScaleX, _pixelScaleY);
    msAlgorithm->SetManualScaleList(_settings.multiscaleScaleList);
    msAlgorithm->SetMultiscaleScaleBias(
        _settings.multiscaleDeconvolutionScaleBias);
    msAlgorithm->SetMaxScales(_settings.multiscaleMaxScales);
    msAlgorithm->SetMultiscaleGain(_settings.multiscaleGain);
    msAlgorithm->SetShape(_settings.multiscaleShapeFunction);
    msAlgorithm->SetTrackComponents(_settings.saveSourceList);
    msAlgorithm->SetConvolutionPadding(_settings.multiscaleConvolutionPadding);
    msAlgorithm->SetUseFastSubMinorLoop(_settings.multiscaleFastSubMinorLoop);
    algorithm.reset(msAlgorithm);
  } else {
    algorithm.reset(new GenericClean(_parallelDeconvolution.GetFFTWManager(),
                                     _settings.useSubMinorOptimization));
  }

  algorithm->SetMaxNIter(_settings.deconvolutionIterationCount);
  algorithm->SetThreshold(_settings.deconvolutionThreshold);
  algorithm->SetGain(_settings.deconvolutionGain);
  algorithm->SetMGain(_settings.deconvolutionMGain);
  algorithm->SetCleanBorderRatio(_settings.deconvolutionBorderRatio);
  algorithm->SetAllowNegativeComponents(_settings.allowNegativeComponents);
  algorithm->SetStopOnNegativeComponents(_settings.stopOnNegativeComponents);
  algorithm->SetThreadCount(threadCount);
  algorithm->SetSpectralFittingMode(_settings.spectralFittingMode,
                                    _settings.spectralFittingTerms);

  ImageSet::CalculateDeconvolutionFrequencies(
      groupTable, _channelFrequencies, _channelWeights,
      _settings.deconvolutionChannelCount);
  algorithm->InitializeFrequencies(_channelFrequencies, _channelWeights);
  _parallelDeconvolution.SetAlgorithm(std::move(algorithm));

  readMask(groupTable);
}

void Deconvolution::readMask(const ImagingTable& groupTable) {
  bool hasMask = false;
  if (!_settings.fitsDeconvolutionMask.empty()) {
    FitsReader maskReader(_settings.fitsDeconvolutionMask, true, true);
    if (maskReader.ImageWidth() != _imgWidth ||
        maskReader.ImageHeight() != _imgHeight)
      throw std::runtime_error(
          "Specified Fits file mask did not have same dimensions as output "
          "image!");
    aocommon::UVector<float> maskData(_imgWidth * _imgHeight);
    if (maskReader.NFrequencies() == 1) {
      Logger::Debug << "Reading mask '" << _settings.fitsDeconvolutionMask
                    << "'...\n";
      maskReader.Read(maskData.data());
    } else if (maskReader.NFrequencies() == _settings.channelsOut) {
      Logger::Debug << "Reading mask '" << _settings.fitsDeconvolutionMask
                    << "' (" << (groupTable.Front().outputChannelIndex + 1)
                    << ")...\n";
      maskReader.ReadIndex(maskData.data(),
                           groupTable.Front().outputChannelIndex);
    } else {
      std::stringstream msg;
      msg << "The number of frequencies in the specified fits mask ("
          << maskReader.NFrequencies()
          << ") does not match the number of requested output channels ("
          << _settings.channelsOut << ")";
      throw std::runtime_error(msg.str());
    }
    _cleanMask.assign(_imgWidth * _imgHeight, false);
    for (size_t i = 0; i != _imgWidth * _imgHeight; ++i)
      _cleanMask[i] = (maskData[i] != 0.0);

    hasMask = true;
  } else if (!_settings.casaDeconvolutionMask.empty()) {
    if (_cleanMask.empty()) {
      Logger::Info << "Reading CASA mask '" << _settings.casaDeconvolutionMask
                   << "'...\n";
      _cleanMask.assign(_imgWidth * _imgHeight, false);
      CasaMaskReader maskReader(_settings.casaDeconvolutionMask);
      if (maskReader.Width() != _imgWidth || maskReader.Height() != _imgHeight)
        throw std::runtime_error(
            "Specified CASA mask did not have same dimensions as output "
            "image!");
      maskReader.Read(_cleanMask.data());
    }

    hasMask = true;
  }

  if (_settings.horizonMask) {
    if (!hasMask) {
      _cleanMask.assign(_imgWidth * _imgHeight, true);
      hasMask = true;
    }

    double fovSq = M_PI * 0.5 - _settings.horizonMaskDistance;
    if (fovSq < 0.0) fovSq = 0.0;
    if (fovSq <= M_PI * 0.5)
      fovSq = std::sin(fovSq);
    else  // a negative horizon distance was given
      fovSq = 1.0 - _settings.horizonMaskDistance;
    fovSq = fovSq * fovSq;
    bool* ptr = _cleanMask.data();

    for (size_t y = 0; y != _imgHeight; ++y) {
      for (size_t x = 0; x != _imgWidth; ++x) {
        double l, m;
        ImageCoordinates::XYToLM(x, y, _pixelScaleX, _pixelScaleY, _imgWidth,
                                 _imgHeight, l, m);
        if (l * l + m * m >= fovSq) *ptr = false;
        ++ptr;
      }
    }

    Logger::Info << "Saving horizon mask...\n";
    Image image(_imgWidth, _imgHeight);
    for (size_t i = 0; i != _imgWidth * _imgHeight; ++i)
      image[i] = _cleanMask[i] ? 1.0 : 0.0;

    FitsWriter writer;
    writer.SetImageDimensions(_imgWidth, _imgHeight, _settings.pixelScaleX,
                              _settings.pixelScaleY);
    writer.Write(_settings.prefixName + "-horizon-mask.fits", image.data());
  }

  if (hasMask) _parallelDeconvolution.SetCleanMask(_cleanMask.data());
}
