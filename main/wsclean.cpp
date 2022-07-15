#include "wsclean.h"

#include "../math/imageoperations.h"

#include "../structures/primarybeam.h"

#include "../gridding/directmsgridder.h"

#include "../io/componentlistwriter.h"
#include "../io/facetreader.h"
#include "../io/imagefilename.h"
#include "../io/imageweightcache.h"
#include "../io/wscfitswriter.h"

#include "../scheduling/griddingtaskmanager.h"

#include "../system/application.h"

#include "../structures/imageweights.h"
#include "../structures/msselection.h"

#include <radler/radler.h>

#include "../idg/averagebeam.h"
#include "../idg/idgmsgridder.h"

#include "../math/renderer.h"

#include "../model/model.h"

#include "../msproviders/contiguousms.h"
#include "../msproviders/msdatadescription.h"

#include "progressbar.h"

#include <aocommon/fits/fitswriter.h>
#include <aocommon/image.h>
#include <aocommon/logger.h>
#include <aocommon/uvector.h>
#include <aocommon/parallelfor.h>
#include <aocommon/units/angle.h>

#include <schaapcommon/facets/facetimage.h>
#include <schaapcommon/fft/resampler.h>
#include <schaapcommon/fft/restoreimage.h>
#include <schaapcommon/fitters/nlplfitter.h>

#include <algorithm>
#include <iostream>
#include <memory>

using aocommon::Image;
using aocommon::Logger;
using aocommon::Polarization;
using aocommon::PolarizationEnum;
using aocommon::units::Angle;

WSClean::WSClean()
    : _globalSelection(),
      _commandLine(),
      _inversionWatch(false),
      _predictingWatch(false),
      _deconvolutionWatch(false),
      _isFirstInversion(true),
      _majorIterationNr(0),
      _psfImages(),
      _modelImages(),
      _residualImages(),
      _deconvolution(),
      _lastStartTime(0.0) {}

WSClean::~WSClean() {}

void WSClean::multiplyImage(double factor, double* image) const {
  if (factor != 1.0) {
    size_t nPix = _settings.trimmedImageWidth * _settings.trimmedImageHeight;
    for (size_t i = 0; i != nPix; ++i) image[i] *= factor;
  }
}

GriddingResult WSClean::loadExistingImage(ImagingTableEntry& entry,
                                          bool isPSF) {
  std::string name;
  if (isPSF) {
    Settings modifiedSettings(_settings);
    modifiedSettings.prefixName = _settings.reusePsfPrefix;
    name =
        ImageFilename::GetPSFPrefix(modifiedSettings, entry.outputChannelIndex,
                                    entry.outputIntervalIndex) +
        "-psf.fits";
  } else {
    Settings modifiedSettings(_settings);
    modifiedSettings.prefixName = _settings.reuseDirtyPrefix;
    name = ImageFilename::GetPrefix(modifiedSettings, entry.polarization,
                                    entry.outputChannelIndex,
                                    entry.outputIntervalIndex, false) +
           "-dirty.fits";
  }
  aocommon::FitsReader reader(name);
  if (reader.ImageWidth() != _settings.trimmedImageWidth ||
      reader.ImageHeight() != _settings.trimmedImageHeight)
    throw std::runtime_error(
        "Image width and height of reused PSF don't match with given settings");
  Image psfImage(reader.ImageWidth(), reader.ImageHeight());
  reader.Read(psfImage.Data());

  GriddingResult result;
  result.images = {std::move(psfImage)};
  result.imageWeight = reader.ReadDoubleKey("WSCIMGWG");
  reader.ReadDoubleKeyIfExists("WSCVWSUM", result.visibilityWeightSum);
  double nVis = 0.0;
  reader.ReadDoubleKeyIfExists("WSCNVIS", nVis);
  result.griddedVisibilityCount = nVis;
  reader.ReadDoubleKeyIfExists("WSCENVIS",
                               result.effectiveGriddedVisibilityCount);
  return result;
}

void WSClean::loadExistingPSF(ImagingTableEntry& entry) {
  Logger::Info << "Loading existing PSF from disk...\n";
  GriddingResult result = loadExistingImage(entry, true);

  if (griddingUsesATerms()) {
    ImageFilename imageName(entry.outputChannelIndex,
                            entry.outputIntervalIndex);

    const std::string beamImageName =
        imageName.GetBeamPrefix(_settings) + ".fits";
    if (!boost::filesystem::exists(beamImageName)) {
      throw std::runtime_error(
          "When reuse-psf is used in combination with the IDG gridder, please "
          "make sure that a corresponding beam image can be found. Expected "
          "file name: " +
          beamImageName);
    }
  }
  const bool writeBeamImage = false;
  imagePSFCallback(entry, result, writeBeamImage);
}

void WSClean::loadExistingDirty(ImagingTableEntry& entry, bool updateBeamInfo) {
  Logger::Info << "Loading existing dirty image from disk...\n";
  GriddingResult result = loadExistingImage(entry, false);
  imageMainCallback(entry, result, updateBeamInfo, true);
}

void WSClean::storeAverageBeam(const ImagingTableEntry& entry,
                               std::unique_ptr<AverageBeam>& averageBeam) {
  if (averageBeam) {
    _scalarBeamImages.SetFitsWriter(
        createWSCFitsWriter(entry, false, false, true).Writer());
    _matrixBeamImages.SetFitsWriter(
        createWSCFitsWriter(entry, false, false, true).Writer());
    averageBeam->Store(_scalarBeamImages, _matrixBeamImages,
                       entry.outputChannelIndex);
  }
}

void WSClean::imagePSF(ImagingTableEntry& entry) {
  Logger::Info.Flush();
  Logger::Info << " == Constructing PSF ==\n";

  GriddingTask task;
  task.operation = GriddingTask::Invert;
  task.imagePSF = true;
  task.polarization = entry.polarization;
  task.subtractModel = false;
  task.verbose = _isFirstInversion;
  task.cache = std::move(_msGridderMetaCache[entry.index]);
  task.storeImagingWeights = _settings.writeImagingWeightSpectrumColumn;
  task.observationInfo = _observationInfo;
  task.facet = entry.facet;
  task.facetIndex = entry.facetIndex;
  task.facetGroupIndex = entry.facetGroupIndex;
  applyFacetPhaseShift(entry, task.observationInfo);
  initializeMSList(entry, task.msList);
  task.imageWeights = initializeImageWeights(entry, task.msList);
  // during PSF imaging, the average beam will never exist, so it is not
  // necessary to set task.averageBeam

  const bool writeBeamImage = true;
  _griddingTaskManager->Run(std::move(task),
                            [this, &entry](GriddingResult& result) {
                              imagePSFCallback(entry, result, writeBeamImage);
                            });
}

void WSClean::imagePSFCallback(ImagingTableEntry& entry, GriddingResult& result,
                               bool writeBeamImage) {
  const size_t channelIndex = entry.outputChannelIndex;
  entry.imageWeight = result.imageWeight;
  entry.normalizationFactor = result.normalizationFactor;
  _infoPerChannel[channelIndex].beamSizeEstimate = result.beamSize;
  _infoPerChannel[channelIndex].weight = entry.imageWeight;
  _infoPerChannel[channelIndex].normalizationFactor = entry.normalizationFactor;
  _infoPerChannel[channelIndex].wGridSize = result.actualWGridSize;
  _infoPerChannel[channelIndex].visibilityCount = result.griddedVisibilityCount;
  _infoPerChannel[channelIndex].effectiveVisibilityCount =
      result.effectiveGriddedVisibilityCount;
  _infoPerChannel[channelIndex].visibilityWeightSum =
      result.visibilityWeightSum;

  if (0 == _facetCount) processFullPSF(result.images[0], entry);

  _lastStartTime = result.startTime;
  _msGridderMetaCache[entry.index] = std::move(result.cache);

  _psfImages.SetFitsWriter(
      createWSCFitsWriter(entry, false, false, false).Writer());
  _psfImages.StoreFacet(result.images[0].Data(),
                        *_settings.polarizations.begin(), channelIndex,
                        entry.facetIndex, entry.facet, false);

  if (writeBeamImage && griddingUsesATerms()) {
    Logger::Info << "Writing IDG beam image...\n";
    ImageFilename imageName(entry.outputChannelIndex,
                            entry.outputIntervalIndex);
    if (!result.averageBeam || result.averageBeam->Empty()) {
      throw std::runtime_error(
          "Trying to write the IDG beam while the beam has not been computed "
          "yet.");
    }
    IdgMsGridder::SaveBeamImage(
        entry, imageName, _settings, _observationInfo.phaseCentreRA,
        _observationInfo.phaseCentreDec, _observationInfo.shiftL,
        _observationInfo.shiftM, *result.averageBeam);
  }

  _isFirstInversion = false;

  storeAverageBeam(entry, result.averageBeam);
}

void WSClean::processFullPSF(Image& image, const ImagingTableEntry& entry) {
  size_t centralIndex =
      _settings.trimmedImageWidth / 2 +
      (_settings.trimmedImageHeight / 2) * _settings.trimmedImageWidth;
  double normFactor;
  if (image[centralIndex] != 0.0)
    normFactor = 1.0 / image[centralIndex];
  else
    normFactor = 0.0;

  const size_t channelIndex = entry.outputChannelIndex;
  _infoPerChannel[channelIndex].psfNormalizationFactor = normFactor;
  image *= normFactor * entry.siCorrection;
  Logger::Debug << "Normalized PSF by factor of " << normFactor << ".\n";

  image.RemoveNans();
  double minPixelScale = std::min(_settings.pixelScaleX, _settings.pixelScaleY);
  double initialFitSize =
      std::max(_infoPerChannel[channelIndex].beamSizeEstimate, minPixelScale);
  double bMaj, bMin, bPA, bTheoretical;
  ImageOperations::DetermineBeamSize(_settings, bMaj, bMin, bPA, bTheoretical,
                                     image.Data(), initialFitSize);
  _infoPerChannel[channelIndex].theoreticBeamSize = bTheoretical;
  _infoPerChannel[channelIndex].beamMaj = bMaj;
  _infoPerChannel[channelIndex].beamMin = bMin;
  _infoPerChannel[channelIndex].beamPA = bPA;

  Logger::Info << "Writing psf image... ";
  if (_settings.isUVImageSaved) {
    saveUVImage(image, entry, false, "uvpsf");
  }

  Logger::Info.Flush();
  const std::string name(
      ImageFilename::GetPSFPrefix(_settings, channelIndex,
                                  entry.outputIntervalIndex) +
      "-psf.fits");
  WSCFitsWriter fitsFile = createWSCFitsWriter(entry, false, false, true);
  fitsFile.WritePSF(name, image.Data());
  Logger::Info << "DONE\n";
}

void WSClean::imageMain(ImagingTableEntry& entry, bool isFirstInversion,
                        bool updateBeamInfo) {
  Logger::Info.Flush();
  Logger::Info << " == Constructing image ==\n";

  GriddingTask task;
  task.operation = GriddingTask::Invert;
  task.imagePSF = false;
  if (_settings.gridderType == GridderType::IDG &&
      _settings.polarizations.size() != 1)
    task.polarization = Polarization::FullStokes;
  else
    task.polarization = entry.polarization;
  task.subtractModel =
      !isFirstInversion || _settings.subtractModel || _settings.continuedRun;
  task.verbose = isFirstInversion && _isFirstInversion;
  task.cache = std::move(_msGridderMetaCache[entry.index]);
  task.storeImagingWeights =
      isFirstInversion && _settings.writeImagingWeightSpectrumColumn;
  initializeMSList(entry, task.msList);
  task.imageWeights = initializeImageWeights(entry, task.msList);
  task.observationInfo = _observationInfo;
  task.facet = entry.facet;
  task.facetIndex = entry.facetIndex;
  task.facetGroupIndex = entry.facetGroupIndex;
  task.averageBeam = AverageBeam::Load(_scalarBeamImages, _matrixBeamImages,
                                       entry.outputChannelIndex);

  applyFacetPhaseShift(entry, task.observationInfo);

  _griddingTaskManager->Run(
      std::move(task),
      [this, &entry, updateBeamInfo, isFirstInversion](GriddingResult& result) {
        imageMainCallback(entry, result, updateBeamInfo, isFirstInversion);
      });
}

void WSClean::imageMainCallback(ImagingTableEntry& entry,
                                GriddingResult& result, bool updateBeamInfo,
                                bool isInitialInversion) {
  size_t joinedChannelIndex = entry.outputChannelIndex;

  _msGridderMetaCache[entry.index] = std::move(result.cache);
  entry.imageWeight = result.imageWeight;
  entry.normalizationFactor = result.normalizationFactor;
  _infoPerChannel[entry.outputChannelIndex].weight = result.imageWeight;
  _infoPerChannel[entry.outputChannelIndex].normalizationFactor =
      result.normalizationFactor;

  _lastStartTime = result.startTime;

  // If no PSF is made, also set the beam size. If the PSF was made, these
  // would already be set after imaging the PSF.
  if (updateBeamInfo) {
    if (_settings.theoreticBeam) {
      _infoPerChannel[entry.outputChannelIndex].beamMaj =
          std::max(result.beamSize, _settings.gaussianTaperBeamSize);
      _infoPerChannel[entry.outputChannelIndex].beamMin =
          std::max(result.beamSize, _settings.gaussianTaperBeamSize);
      _infoPerChannel[entry.outputChannelIndex].beamPA = 0.0;
    } else if (_settings.manualBeamMajorSize != 0.0) {
      _infoPerChannel[entry.outputChannelIndex].beamMaj =
          _settings.manualBeamMajorSize;
      _infoPerChannel[entry.outputChannelIndex].beamMin =
          _settings.manualBeamMinorSize;
      _infoPerChannel[entry.outputChannelIndex].beamPA = _settings.manualBeamPA;
    } else {
      _infoPerChannel[entry.outputChannelIndex].beamMaj =
          std::numeric_limits<double>::quiet_NaN();
      _infoPerChannel[entry.outputChannelIndex].beamMin =
          std::numeric_limits<double>::quiet_NaN();
      _infoPerChannel[entry.outputChannelIndex].beamPA =
          std::numeric_limits<double>::quiet_NaN();
    }
  }

  using PolImagesPair = std::pair<const PolarizationEnum, std::vector<Image>>;
  std::vector<PolImagesPair> imageList;
  if (_settings.gridderType == GridderType::IDG &&
      _settings.polarizations.size() != 1) {
    assert(result.images.size() == _settings.polarizations.size());
    imageList.reserve(result.images.size());
    auto polIter = _settings.polarizations.begin();
    for (size_t polIndex = 0; polIndex != result.images.size();
         ++polIndex, ++polIter)
      imageList.emplace_back(
          *polIter, std::vector<Image>{std::move(result.images[polIndex])});
  } else {
    imageList.emplace_back(entry.polarization, std::move(result.images));
  }

  for (PolImagesPair& polImagePair : imageList) {
    std::vector<Image>& images = polImagePair.second;
    const PolarizationEnum polarization = polImagePair.first;
    for (size_t i = 0; i != images.size(); ++i) {
      // IDG performs normalization on the dirty images, so only normalize if
      // not using IDG
      const double psfFactor =
          _settings.gridderType == GridderType::IDG
              ? 1.0
              : _infoPerChannel[joinedChannelIndex].psfNormalizationFactor;
      if (_settings.gridderType != GridderType::IDG)
        images[i] *= psfFactor * entry.siCorrection;
      const bool isImaginary = i == 1;
      storeAndCombineXYandYX(_residualImages, joinedChannelIndex, entry,
                             polarization, isImaginary, images[i].Data());
    }

    // If facets are used, stitchFacets() performs these actions.
    if (isInitialInversion && 0 == _facetCount) {
      // maxFacetGroupIndex is always 1
      const size_t maxFacetGroupIndex = 1;
      initializeModelImages(entry, polarization, maxFacetGroupIndex);

      _residualImages.SetFitsWriter(
          createWSCFitsWriter(entry, polarization, false, false, false)
              .Writer());
      // If facets are used, stitchFacets() saves the dirty image.
      if (_settings.isDirtySaved) {
        for (size_t imageIndex = 0; imageIndex != entry.imageCount;
             ++imageIndex) {
          const bool isImaginary = (imageIndex == 1);
          WSCFitsWriter writer(createWSCFitsWriter(entry, polarization,
                                                   isImaginary, false, false));
          Image dirtyImage(_settings.trimmedImageWidth,
                           _settings.trimmedImageHeight);
          _residualImages.Load(dirtyImage.Data(), polarization,
                               entry.outputChannelIndex, isImaginary);
          Logger::Info << "Writing dirty image...\n";
          writer.WriteImage("dirty.fits", dirtyImage.Data());
        }
      }
    }
  }

  storeAverageBeam(entry, result.averageBeam);
}

void WSClean::storeAndCombineXYandYX(CachedImageSet& dest,
                                     size_t joinedChannelIndex,
                                     const ImagingTableEntry& entry,
                                     PolarizationEnum polarization,
                                     bool isImaginary, const float* image) {
  if (polarization == Polarization::YX &&
      _settings.polarizations.count(Polarization::XY) != 0) {
    Logger::Info << "Adding XY and YX together...\n";
    Image xyImage;
    if (entry.facet) {
      // Trimmed facet size
      xyImage = Image(entry.facet->GetTrimmedBoundingBox().Width(),
                      entry.facet->GetTrimmedBoundingBox().Height());
    } else {
      // Full image size
      xyImage =
          Image(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
    }

    dest.LoadFacet(xyImage.Data(), Polarization::XY, joinedChannelIndex,
                   entry.facetIndex, entry.facet, isImaginary);
    if (isImaginary) {
      for (size_t i = 0; i != xyImage.Size(); ++i)
        xyImage[i] = (xyImage[i] - image[i]) * 0.5;
    } else {
      for (size_t i = 0; i != xyImage.Size(); ++i)
        xyImage[i] = (xyImage[i] + image[i]) * 0.5;
    }
    dest.StoreFacet(xyImage.Data(), Polarization::XY, joinedChannelIndex,
                    entry.facetIndex, entry.facet, isImaginary);
  } else {
    dest.StoreFacet(image, polarization, joinedChannelIndex, entry.facetIndex,
                    entry.facet, isImaginary);
  }
}

void WSClean::predict(const ImagingTableEntry& entry) {
  Logger::Info.Flush();
  Logger::Info << " == Converting model image to visibilities ==\n";
  size_t width;
  size_t height;
  if (entry.facet) {
    width = entry.facet->GetTrimmedBoundingBox().Width();
    height = entry.facet->GetTrimmedBoundingBox().Height();
  } else {
    width = _settings.trimmedImageWidth;
    height = _settings.trimmedImageHeight;
  }
  const bool isFullStokes = _settings.gridderType == GridderType::IDG &&
                            _settings.polarizations.size() != 1;
  std::vector<PolarizationEnum> polarizations;
  if (isFullStokes)
    polarizations.assign(_settings.polarizations.begin(),
                         _settings.polarizations.end());
  else
    polarizations = {entry.polarization};

  std::vector<Image> modelImages;
  modelImages.reserve(polarizations.size());
  for (PolarizationEnum& polarization : polarizations) {
    modelImages.emplace_back(width, height);
    bool isYX = polarization == Polarization::YX;
    const PolarizationEnum loadPol = isYX ? Polarization::XY : polarization;
    _modelImages.LoadFacet(modelImages.back().Data(), loadPol,
                           entry.outputChannelIndex, entry.facetIndex,
                           entry.facet, false);
    if (Polarization::IsComplex(polarization)) {  // XY or YX
      modelImages.emplace_back(width, height);
      // YX is never stored: it is always combined with XY and stored as XY
      _modelImages.LoadFacet(modelImages.back().Data(), Polarization::XY,
                             entry.outputChannelIndex, entry.facetIndex,
                             entry.facet, true);
      if (isYX) {
        for (float& v : modelImages.back()) v = -v;
      }
    }
  }
  GriddingTask task;
  task.operation = GriddingTask::Predict;
  task.polarization =
      isFullStokes ? Polarization::FullStokes : entry.polarization;
  task.cache = std::move(_msGridderMetaCache[entry.index]);
  task.verbose = false;
  task.storeImagingWeights = false;
  task.modelImages = std::move(modelImages);
  initializeMSList(entry, task.msList);
  task.imageWeights = initializeImageWeights(entry, task.msList);
  task.observationInfo = _observationInfo;
  task.facet = entry.facet;
  task.facetIndex = entry.facetIndex;
  task.facetGroupIndex = entry.facetGroupIndex;
  task.averageBeam = AverageBeam::Load(_scalarBeamImages, _matrixBeamImages,
                                       entry.outputChannelIndex);
  applyFacetPhaseShift(entry, task.observationInfo);
  _griddingTaskManager->Run(
      std::move(task), [this, &entry](GriddingResult& result) {
        _msGridderMetaCache[entry.index] = std::move(result.cache);
      });
}

ObservationInfo WSClean::getObservationInfo() const {
  casacore::MeasurementSet ms(_settings.filenames[0]);
  ObservationInfo observationInfo =
      ReadObservationInfo(ms, _settings.fieldIds[0]);
  if (_settings.hasShift) {
    observationInfo.hasShiftedPhaseCentre = true;
    aocommon::ImageCoordinates::RaDecToLM(
        _settings.shiftRA, _settings.shiftDec, observationInfo.phaseCentreRA,
        observationInfo.phaseCentreDec, observationInfo.shiftL,
        observationInfo.shiftM);
  }
  return observationInfo;
}

void WSClean::applyFacetPhaseShift(const ImagingTableEntry& entry,
                                   ObservationInfo& observationInfo) const {
  if (entry.facet) {
    observationInfo.shiftL -= entry.centreShiftX * _settings.pixelScaleX;
    observationInfo.shiftM += entry.centreShiftY * _settings.pixelScaleY;
    if (entry.centreShiftX != 0.0 || entry.centreShiftY != 0.0) {
      observationInfo.hasShiftedPhaseCentre = true;
    }
  }
}

std::shared_ptr<ImageWeights> WSClean::initializeImageWeights(
    const ImagingTableEntry& entry,
    std::vector<std::unique_ptr<MSDataDescription>>& msList) {
  if (_settings.mfWeighting) {
    return _imageWeightCache->GetMFWeights();
  } else {
    std::shared_ptr<ImageWeights> weights = _imageWeightCache->Get(
        msList, entry.outputChannelIndex, entry.outputIntervalIndex);
    if (_settings.isWeightImageSaved) {
      std::string prefix = ImageFilename::GetPSFPrefix(
          _settings, entry.outputChannelIndex, entry.outputIntervalIndex);
      weights->Save(prefix + "-weights.fits");
    }
    return weights;
  }
}

void WSClean::initializeMFSImageWeights() {
  Logger::Info << "Precalculating MF weights for "
               << _settings.weightMode.ToString() << " weighting...\n";
  std::unique_ptr<ImageWeights> weights = _imageWeightCache->MakeEmptyWeights();
  if (_settings.doReorder) {
    for (const ImagingTable::Group& sqGroup : _imagingTable.SquaredGroups()) {
      const ImagingTableEntry& entry = *sqGroup.front();
      for (size_t msIndex = 0; msIndex != _settings.filenames.size();
           ++msIndex) {
        const ImagingTableEntry::MSInfo& ms = entry.msData[msIndex];
        for (size_t dataDescId = 0;
             dataDescId != _msBands[msIndex].DataDescCount(); ++dataDescId) {
          MSSelection partSelection(_globalSelection);
          const bool hasSelection = partSelection.SelectMsChannels(
              _msBands[msIndex], dataDescId, entry);
          if (hasSelection) {
            const PolarizationEnum pol =
                _settings.gridderType == GridderType::IDG ? getIdgPolarization()
                                                          : entry.polarization;
            PartitionedMS msProvider(_partitionedMSHandles[msIndex],
                                     ms.bands[dataDescId].partIndex, pol,
                                     dataDescId);
            aocommon::BandData selectedBand(_msBands[msIndex][dataDescId]);
            if (partSelection.HasChannelRange()) {
              selectedBand = aocommon::BandData(
                  selectedBand, partSelection.ChannelRangeStart(),
                  partSelection.ChannelRangeEnd());
            }
            weights->Grid(msProvider, selectedBand);
          }
        }
      }
    }
  } else {
    for (size_t i = 0; i != _settings.filenames.size(); ++i) {
      for (size_t d = 0; d != _msBands[i].DataDescCount(); ++d) {
        const PolarizationEnum pol = _settings.gridderType == GridderType::IDG
                                         ? getIdgPolarization()
                                         : *_settings.polarizations.begin();
        ContiguousMS msProvider(_settings.filenames[i],
                                _settings.dataColumnName, _globalSelection, pol,
                                d, _settings.useMPI);
        aocommon::BandData selectedBand = _msBands[i][d];
        if (_globalSelection.HasChannelRange())
          selectedBand = aocommon::BandData(
              selectedBand, _globalSelection.ChannelRangeStart(),
              _globalSelection.ChannelRangeEnd());
        weights->Grid(msProvider, selectedBand);
        Logger::Info << '.';
        Logger::Info.Flush();
      }
    }
  }
  weights->FinishGridding();
  _imageWeightCache->SetMFWeights(std::move(weights));
  if (_settings.isWeightImageSaved)
    _imageWeightCache->GetMFWeights()->Save(_settings.prefixName +
                                            "-weights.fits");
}

void WSClean::performReordering(bool isPredictMode) {
  std::mutex mutex;

  // If there are reordered measurement sets on disk, we have to clean them
  // before writing new ones:
  _partitionedMSHandles.clear();

  _partitionedMSHandles.resize(_settings.filenames.size());
  bool useModel = _settings.deconvolutionMGain != 1.0 || isPredictMode ||
                  _settings.subtractModel || _settings.continuedRun;
  bool initialModelRequired = _settings.subtractModel || _settings.continuedRun;

  if (_settings.parallelReordering != 1) Logger::Info << "Reordering...\n";

  aocommon::ParallelFor<size_t> loop(_settings.parallelReordering);
  loop.Run(0, _settings.filenames.size(), [&](size_t msIndex, size_t) {
    std::vector<PartitionedMS::ChannelRange> channels;
    std::map<PolarizationEnum, size_t> nextIndex;
    for (size_t sqIndex = 0; sqIndex != _imagingTable.SquaredGroupCount();
         ++sqIndex) {
      ImagingTable sqGroup = _imagingTable.GetSquaredGroup(sqIndex);
      for (size_t fgIndex = 0; fgIndex != sqGroup.FacetGroupCount();
           ++fgIndex) {
        ImagingTable facetGroup = sqGroup.GetFacetGroup(fgIndex);
        // The band information is determined from the first facet in the group.
        // After this, all facet entries inside the group are updated.
        const ImagingTableEntry& entry = facetGroup.Front();
        for (size_t d = 0; d != _msBands[msIndex].DataDescCount(); ++d) {
          MSSelection selection(_globalSelection);
          if (selection.SelectMsChannels(_msBands[msIndex], d, entry)) {
            if (entry.polarization == *_settings.polarizations.begin()) {
              PartitionedMS::ChannelRange r;
              r.dataDescId = d;
              r.start = selection.ChannelRangeStart();
              r.end = selection.ChannelRangeEnd();
              channels.push_back(r);
            }
            for (ImagingTableEntry& facetEntry : facetGroup) {
              facetEntry.msData[msIndex].bands[d].partIndex =
                  nextIndex[entry.polarization];
            }
            ++nextIndex[entry.polarization];
          }
        }
      }
    }

    PartitionedMS::Handle partMS = PartitionedMS::Partition(
        _settings.filenames[msIndex], channels, _globalSelection,
        _settings.dataColumnName, useModel, initialModelRequired, _settings);
    std::lock_guard<std::mutex> lock(mutex);
    _partitionedMSHandles[msIndex] = std::move(partMS);
    if (_settings.parallelReordering != 1)
      Logger::Info << "Finished reordering " << _settings.filenames[msIndex]
                   << " [" << msIndex << "]\n";
  });
}

void WSClean::RunClean() {
  _observationInfo = getObservationInfo();
  std::vector<std::shared_ptr<schaapcommon::facets::Facet>> facets =
      FacetReader::ReadFacets(_settings, _observationInfo);
  _facetCount = facets.size();

  schaapcommon::facets::Pixel centerPixel(_settings.trimmedImageWidth / 2,
                                          _settings.trimmedImageHeight / 2);
  const bool hasCenter = std::any_of(
      facets.begin(), facets.end(),
      [&centerPixel](
          const std::shared_ptr<schaapcommon::facets::Facet>& facet) {
        // Point-in-poly test only evaluated if bounding box does
        // contain the centerPixel
        return facet->GetTrimmedBoundingBox().Contains(centerPixel) &&
               facet->Contains(centerPixel);
      });

  // FIXME: raise warning if facets do not cover the entire image, see AST-429

  // Center pixel should be present in one of the facets for the deconvolution
  if (!facets.empty() && _settings.deconvolutionIterationCount > 0 &&
      !hasCenter) {
    throw std::runtime_error(
        "The center pixel of the full image is not found in one of the facets. "
        "Make sure your facet file defines a facet that covers the center "
        "pixel of the main image.");
  }

  _globalSelection = _settings.GetMSSelection();
  MSSelection fullSelection = _globalSelection;

  for (size_t intervalIndex = 0; intervalIndex != _settings.intervalsOut;
       ++intervalIndex) {
    makeImagingTable(intervalIndex);
    if (!facets.empty()) updateFacetsInImagingTable(facets);

    _globalSelection = selectInterval(fullSelection, intervalIndex);

    if (_settings.doReorder) performReordering(false);

    _infoPerChannel.assign(_settings.channelsOut, OutputChannelInfo());

    _msGridderMetaCache.clear();
    _imageWeightCache = createWeightCache();

    if (_settings.mfWeighting) initializeMFSImageWeights();

    _griddingTaskManager = GriddingTaskManager::Make(_settings);

    std::unique_ptr<PrimaryBeam> primaryBeam;
    for (size_t groupIndex = 0;
         groupIndex != _imagingTable.IndependentGroupCount(); ++groupIndex) {
      ImagingTable group = _imagingTable.GetIndependentGroup(groupIndex);
      runIndependentGroup(group, primaryBeam);
    }

    // Needs to be destructed before image allocator, or image allocator will
    // report error caused by leaked memory
    _griddingTaskManager.reset();

    if (_settings.channelsOut > 1) {
      for (PolarizationEnum pol : _settings.polarizations) {
        bool psfWasMade = (_settings.deconvolutionIterationCount > 0 ||
                           _settings.makePSF || _settings.makePSFOnly) &&
                          pol == *_settings.polarizations.begin();

        if (psfWasMade) {
          ImageOperations::MakeMFSImage(_settings, _infoPerChannel, _infoForMFS,
                                        "psf.fits", intervalIndex, pol, false,
                                        true);
          if (_settings.savePsfPb)
            ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                          _infoForMFS, "psf-pb.fits",
                                          intervalIndex, pol, false, true);
        }
        if (griddingUsesATerms()) {
          ImageOperations::MakeMFSImage(_settings, _infoPerChannel, _infoForMFS,
                                        "beam.fits", intervalIndex, pol, false,
                                        false);
        } else if (usesBeam()) {
          // The (complex valued but Hermitian) Mueller matrices are stored with
          // 16 elements:
          constexpr size_t n_matrix_elements = 16;
          for (size_t beam_index = 0; beam_index != n_matrix_elements;
               ++beam_index) {
            ImageOperations::MakeMFSImage(
                _settings, _infoPerChannel, _infoForMFS,
                "beam-" + std::to_string(beam_index) + ".fits", intervalIndex,
                pol, false, false);
          }
        }

        if (!(pol == Polarization::YX &&
              _settings.polarizations.count(Polarization::XY) != 0) &&
            !_settings.makePSFOnly) {
          if (_settings.isDirtySaved)
            ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                          _infoForMFS, "dirty.fits",
                                          intervalIndex, pol, false);
          if (_settings.deconvolutionIterationCount == 0) {
            ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                          _infoForMFS, "image.fits",
                                          intervalIndex, pol, false);
            if (usesBeam())
              ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                            _infoForMFS, "image-pb.fits",
                                            intervalIndex, pol, false);
          } else {
            ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                          _infoForMFS, "residual.fits",
                                          intervalIndex, pol, false);
            ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                          _infoForMFS, "model.fits",
                                          intervalIndex, pol, false);
            ImageOperations::RenderMFSImage(_settings, _infoForMFS,
                                            intervalIndex, pol, false, false);
            if (usesBeam()) {
              ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                            _infoForMFS, "residual-pb.fits",
                                            intervalIndex, pol, false);
              ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                            _infoForMFS, "model-pb.fits",
                                            intervalIndex, pol, false);
              ImageOperations::RenderMFSImage(_settings, _infoForMFS,
                                              intervalIndex, pol, false, true);
            }
          }
          if (Polarization::IsComplex(pol)) {
            if (_settings.isDirtySaved)
              ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                            _infoForMFS, "dirty.fits",
                                            intervalIndex, pol, true);
            if (_settings.deconvolutionIterationCount == 0) {
              ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                            _infoForMFS, "image.fits",
                                            intervalIndex, pol, true);
              if (usesBeam())
                ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                              _infoForMFS, "image-pb.fits",
                                              intervalIndex, pol, true);
            } else {
              ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                            _infoForMFS, "residual.fits",
                                            intervalIndex, pol, true);
              ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                            _infoForMFS, "model.fits",
                                            intervalIndex, pol, true);
              ImageOperations::RenderMFSImage(_settings, _infoForMFS,
                                              intervalIndex, pol, true, false);
              if (usesBeam()) {
                ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                              _infoForMFS, "residual-pb.fits",
                                              intervalIndex, pol, true);
                ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                              _infoForMFS, "model-pb.fits",
                                              intervalIndex, pol, true);
                ImageOperations::RenderMFSImage(_settings, _infoForMFS,
                                                intervalIndex, pol, true, true);
              }
            }
          }
        }
      }
    }

    // This will erase the temporary files
    _partitionedMSHandles.clear();
  }
}

std::unique_ptr<ImageWeightCache> WSClean::createWeightCache() {
  std::unique_ptr<ImageWeightCache> cache(new ImageWeightCache(
      _settings.weightMode, _settings.paddedImageWidth,
      _settings.paddedImageHeight, _settings.pixelScaleX, _settings.pixelScaleY,
      _settings.minUVInLambda, _settings.maxUVInLambda,
      _settings.rankFilterLevel, _settings.rankFilterSize,
      _settings.useWeightsAsTaper, _settings.threadCount));
  cache->SetTaperInfo(
      _settings.gaussianTaperBeamSize, _settings.tukeyTaperInLambda,
      _settings.tukeyInnerTaperInLambda, _settings.edgeTaperInLambda,
      _settings.edgeTukeyTaperInLambda);
  return cache;
}

void WSClean::RunPredict() {
  // When facets are used, the initialization of the imaging table and the
  // facets depend on eachother. We therefore use this approach:
  // 1. Count the number of facets and store in _facetCount.
  // 2. Create the imaging table using _facetCount and set the facet index in
  //    the imaging table entries. Each interval loop iteration creates a new
  //    imaging table.
  // 3. Read the image size and pixel scale from the input fits file
  //    corresponding to the first imaging table entry. This way, the user does
  //    not have to specify these values on the command line.
  // 4. In the first interval loop iteration, update the settings using the
  //    values from the input fits file. In subsequent iterations, check if the
  //    image size and pixel scale match the existing settings.
  // 5. In the first interval loop iteration, create the facets using the new
  //    settings. In subsequent iterations, the settings do not change so
  //    recreating the facets is not needed.
  // 6. Set the facets and related properties in the imaging table entries,
  //    using the existing facet index in the entries.

  assert(!_deconvolution.has_value());
  _observationInfo = getObservationInfo();
  std::vector<std::shared_ptr<schaapcommon::facets::Facet>> facets;
  _facetCount = FacetReader::CountFacets(_settings.facetRegionFilename);

  _globalSelection = _settings.GetMSSelection();
  MSSelection fullSelection = _globalSelection;

  for (size_t intervalIndex = 0; intervalIndex != _settings.intervalsOut;
       ++intervalIndex) {
    makeImagingTable(intervalIndex);

    if (_settings.predictionChannels != 0) {
      // TODO
    }

    _infoPerChannel.assign(_settings.channelsOut, OutputChannelInfo());
    _msGridderMetaCache.clear();

    _globalSelection = selectInterval(fullSelection, intervalIndex);

    if (_settings.doReorder) performReordering(true);

    if (_facetCount != 0) {
      std::string prefix =
          ImageFilename::GetPrefix(_settings, _imagingTable[0].polarization,
                                   _imagingTable[0].outputChannelIndex,
                                   _imagingTable[0].outputIntervalIndex, false);
      const std::string suffix =
          (_settings.applyFacetBeam || !_settings.facetSolutionFiles.empty())
              ? "-model-pb.fits"
              : "-model.fits";
      aocommon::FitsReader reader(prefix + suffix);
      overrideImageSettings(reader);
      if (intervalIndex == 0) {
        facets = FacetReader::ReadFacets(_settings, _observationInfo);
        // FIXME: raise warning if facets do not cover the entire image, see
        // AST-429
      }

      updateFacetsInImagingTable(facets);
    }

    _griddingTaskManager = GriddingTaskManager::Make(_settings);
    predictGroup(_imagingTable);
    _griddingTaskManager.reset();
  }
}

double WSClean::minTheoreticalBeamSize(const ImagingTable& table) const {
  double beam = 0.0;
  for (const ImagingTableEntry& e : table) {
    const OutputChannelInfo& info = _infoPerChannel[e.outputChannelIndex];
    if (std::isfinite(info.theoreticBeamSize) &&
        (info.theoreticBeamSize < beam || beam == 0.0))
      beam = info.theoreticBeamSize;
  }
  return beam;
}

void WSClean::runIndependentGroup(ImagingTable& groupTable,
                                  std::unique_ptr<PrimaryBeam>& primaryBeam) {
  WSCFitsWriter modelWriter(
      createWSCFitsWriter(groupTable.Front(), false, true, false));
  _modelImages.Initialize(modelWriter.Writer(), _settings.polarizations.size(),
                          _settings.channelsOut, _facetCount,
                          _settings.prefixName + "-model");
  WSCFitsWriter writer(
      createWSCFitsWriter(groupTable.Front(), false, false, false));
  _residualImages.Initialize(writer.Writer(), _settings.polarizations.size(),
                             _settings.channelsOut, _facetCount,
                             _settings.prefixName + "-residual");
  if (groupTable.Front().polarization == *_settings.polarizations.begin()) {
    _psfImages.Initialize(writer.Writer(), 1, groupTable.SquaredGroups().size(),
                          _facetCount, _settings.prefixName + "-psf");
    _scalarBeamImages.Initialize(writer.Writer(), 1,
                                 groupTable.SquaredGroups().size(), _facetCount,
                                 _settings.prefixName + "-scalar-beam");
    _matrixBeamImages.Initialize(writer.Writer(), 2,
                                 groupTable.SquaredGroups().size(), _facetCount,
                                 _settings.prefixName + "-matrix-beam");
  }

  // In the case of IDG we have to directly ask for all four polarizations.
  const bool requestPolarizationsAtOnce =
      _settings.gridderType == GridderType::IDG &&
      _settings.polarizations.size() > 1;

  // In case XY/YX polarizations are requested, we should not parallelize over
  // those since they need to be combined after imaging, and this currently
  // requires XY before YX.
  const bool parallelizePolarizations =
      _settings.polarizations.count(Polarization::XY) == 0 &&
      _settings.polarizations.count(Polarization::YX) == 0;

  _inversionWatch.Start();
  const bool doMakePSF = _settings.deconvolutionIterationCount > 0 ||
                         _settings.makePSF || _settings.makePSFOnly;
  for (ImagingTableEntry& entry : groupTable) {
    const bool isFirstPol =
        entry.polarization == *_settings.polarizations.begin();
    if (doMakePSF && isFirstPol) {
      if (_settings.reusePsf)
        loadExistingPSF(entry);
      else
        imagePSF(entry);
    }
  }
  _griddingTaskManager->Finish();
  if (doMakePSF) stitchFacets(groupTable, _psfImages, false, true);

  if (!_settings.makePSFOnly) {
    runFirstInversions(groupTable, primaryBeam, requestPolarizationsAtOnce,
                       parallelizePolarizations);
  }

  _inversionWatch.Pause();

  if (!_settings.makePSFOnly) {
    runMajorIterations(groupTable, primaryBeam, requestPolarizationsAtOnce,
                       parallelizePolarizations);
  }

  Logger::Info << "Inversion: " << _inversionWatch.ToString()
               << ", prediction: " << _predictingWatch.ToString()
               << ", deconvolution: " << _deconvolutionWatch.ToString() << '\n';
}

void WSClean::saveRestoredImagesForGroup(
    const ImagingTable& table,
    std::unique_ptr<PrimaryBeam>& primaryBeam) const {
  const ImagingTableEntry tableEntry = table.Front();
  assert(tableEntry.facetIndex == 0);

  // Restore model to residual and save image
  size_t currentChannelIndex = tableEntry.outputChannelIndex;

  PolarizationEnum curPol = tableEntry.polarization;
  for (size_t imageIter = 0; imageIter != tableEntry.imageCount; ++imageIter) {
    bool isImaginary = (imageIter == 1);
    WSCFitsWriter writer(
        createWSCFitsWriter(tableEntry, isImaginary, false, true));
    Image restoredImage(_settings.trimmedImageWidth,
                        _settings.trimmedImageHeight);
    _residualImages.Load(restoredImage.Data(), curPol, currentChannelIndex,
                         isImaginary);

    if (_settings.deconvolutionIterationCount != 0)
      writer.WriteImage("residual.fits", restoredImage.Data());

    if (_settings.isUVImageSaved)
      saveUVImage(restoredImage, tableEntry, isImaginary, "uv");

    Image modelImage(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
    _modelImages.Load(modelImage.Data(), curPol, currentChannelIndex,
                      isImaginary);
    double beamMaj = _infoPerChannel[currentChannelIndex].beamMaj;
    double beamMin, beamPA;
    std::string beamStr;
    if (std::isfinite(beamMaj)) {
      beamMin = _infoPerChannel[currentChannelIndex].beamMin;
      beamPA = _infoPerChannel[currentChannelIndex].beamPA;
      beamStr = "(beam=" + Angle::ToNiceString(beamMin) + "-" +
                Angle::ToNiceString(beamMaj) +
                ", PA=" + Angle::ToNiceString(beamPA) + ")";
    } else {
      beamStr = "(beam is neither fitted nor estimated -- using delta scales!)";
      beamMaj = 0.0;
      beamMin = 0.0;
      beamPA = 0.0;
    }
    Logger::Info << "Rendering sources to restored image " + beamStr + "... ";
    Logger::Info.Flush();
    schaapcommon::fft::RestoreImage(
        restoredImage.Data(), modelImage.Data(), _settings.trimmedImageWidth,
        _settings.trimmedImageHeight, beamMaj, beamMin, beamPA,
        _settings.pixelScaleX, _settings.pixelScaleY, _settings.threadCount);
    Logger::Info << "DONE\n";
    modelImage.Reset();

    Logger::Info << "Writing restored image... ";
    Logger::Info.Flush();
    writer.WriteImage("image.fits", restoredImage.Data());
    Logger::Info << "DONE\n";
    restoredImage.Reset();

    // H5 corrections will be applied on the beam images
    const bool applyH5OnBeamImages =
        (_settings.applyFacetBeam && !_settings.facetSolutionFiles.empty());
    ImageFilename imageName =
        ImageFilename(currentChannelIndex, tableEntry.outputIntervalIndex);

    // This conditional ensures that the images are available before applying
    // the primarybeam
    if (curPol == *_settings.polarizations.rbegin()) {
      if (griddingUsesATerms()) {
        IdgMsGridder::SavePBCorrectedImages(writer.Writer(), imageName, "image",
                                            _settings);
        if (_settings.savePsfPb)
          IdgMsGridder::SavePBCorrectedImages(writer.Writer(), imageName, "psf",
                                              _settings);
        if (_settings.deconvolutionIterationCount != 0) {
          IdgMsGridder::SavePBCorrectedImages(writer.Writer(), imageName,
                                              "residual", _settings);
          IdgMsGridder::SavePBCorrectedImages(writer.Writer(), imageName,
                                              "model", _settings);
        }
      } else if (_settings.applyPrimaryBeam || _settings.applyFacetBeam) {
        primaryBeam->CorrectImages(writer.Writer(), imageName, "image", table,
                                   _msGridderMetaCache, applyH5OnBeamImages);
        if (_settings.savePsfPb)
          primaryBeam->CorrectImages(writer.Writer(), imageName, "psf", table,
                                     _msGridderMetaCache, applyH5OnBeamImages);
        if (_settings.deconvolutionIterationCount != 0) {
          primaryBeam->CorrectImages(writer.Writer(), imageName, "residual",
                                     table, _msGridderMetaCache,
                                     applyH5OnBeamImages);
          primaryBeam->CorrectImages(writer.Writer(), imageName, "model", table,
                                     _msGridderMetaCache, applyH5OnBeamImages);
        }
      }
    }

    // Apply the H5 solutions to the facets. In case a H5 solution file is
    // provided, but no primary beam correction was applied.
    // This can be done on a per-polarization basis (as long as we do not
    // fully support a Full Jones correction).
    if (!_settings.facetSolutionFiles.empty() && !applyH5OnBeamImages) {
      correctImagesH5(writer.Writer(), table, imageName, "image");
      if (_settings.savePsfPb)
        correctImagesH5(writer.Writer(), table, imageName, "psf");
      if (_settings.deconvolutionIterationCount != 0) {
        correctImagesH5(writer.Writer(), table, imageName, "residual");
        correctImagesH5(writer.Writer(), table, imageName, "model");
      }
    }
  }
}

void WSClean::writeFirstResidualImages(const ImagingTable& groupTable) const {
  Logger::Info << "Writing first iteration image(s)...\n";
  Image ptr(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
  for (const ImagingTableEntry& entry : groupTable) {
    size_t ch = entry.outputChannelIndex;
    if (entry.polarization == Polarization::YX) {
      _residualImages.Load(ptr.Data(), Polarization::XY, ch, true);
      WSCFitsWriter writer(
          createWSCFitsWriter(entry, Polarization::XY, true, false, false));
      writer.WriteImage("first-residual.fits", ptr.Data());
    } else {
      _residualImages.Load(ptr.Data(), entry.polarization, ch, false);
      WSCFitsWriter writer(createWSCFitsWriter(entry, false, false, false));
      writer.WriteImage("first-residual.fits", ptr.Data());
    }
  }
}

void WSClean::writeModelImages(const ImagingTable& groupTable) const {
  Logger::Info << "Writing model image...\n";
  Image ptr(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
  for (const ImagingTableEntry& entry : groupTable) {
    size_t ch = entry.outputChannelIndex;
    if (entry.polarization == Polarization::YX) {
      _modelImages.Load(ptr.Data(), Polarization::XY, ch, true);
      WSCFitsWriter writer(
          createWSCFitsWriter(entry, Polarization::XY, true, true, true));
      writer.WriteImage("model.fits", ptr.Data());
    } else {
      _modelImages.Load(ptr.Data(), entry.polarization, ch, false);
      WSCFitsWriter writer(createWSCFitsWriter(entry, false, true, true));
      writer.WriteImage("model.fits", ptr.Data());
    }
  }
}

void WSClean::partitionModelIntoFacets(const ImagingTable& table,
                                       bool isPredictOnly) {
  if (_facetCount != 0) {
    Logger::Info << "Clipping model image into facets...\n";
    // Allocate full image
    Image fullImage(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
    // Initialize FacetImage with properties of stitched image, always
    // stitch facets for 1 spectral term.
    schaapcommon::facets::FacetImage facetImage(
        _settings.trimmedImageWidth, _settings.trimmedImageHeight, 1);
    for (size_t facetGroupIndex = 0; facetGroupIndex != table.FacetGroupCount();
         ++facetGroupIndex) {
      const ImagingTable clipGroup = table.GetFacetGroup(facetGroupIndex);
      const size_t imageCount = clipGroup.Front().imageCount;
      _modelImages.Load(fullImage.Data(), clipGroup.Front().polarization,
                        clipGroup.Front().outputChannelIndex, false);
      for (size_t imageIndex = 0; imageIndex != imageCount; ++imageIndex) {
        partitionSingleGroup(clipGroup, imageIndex, _modelImages, fullImage,
                             facetImage, isPredictOnly);
      }
    }
  }
}

void WSClean::partitionSingleGroup(const ImagingTable& facetGroup,
                                   size_t imageIndex,
                                   CachedImageSet& imageCache,
                                   const Image& fullImage,
                                   schaapcommon::facets::FacetImage& facetImage,
                                   bool isPredictOnly) {
  const bool isImaginary = (imageIndex == 1);
  for (const ImagingTableEntry& facetEntry : facetGroup) {
    facetImage.SetFacet(*facetEntry.facet, true);
    facetImage.CopyToFacet({fullImage.Data()});
    if (!isPredictOnly) {
      if (_settings.applyFacetBeam || !_settings.facetSolutionFiles.empty()) {
        // Apply average direction dependent correction before storing
        const long double m =
            _msGridderMetaCache[facetEntry.index]->correctionSum /
            facetEntry.imageWeight;
        facetImage *= 1.0f / std::sqrt(m);
      }
    }
    imageCache.StoreFacet(facetImage.Data(0), facetEntry.polarization,
                          facetEntry.outputChannelIndex, facetEntry.facetIndex,
                          facetEntry.facet, isImaginary);
  }
}

void WSClean::initializeModelImages(const ImagingTableEntry& entry,
                                    PolarizationEnum polarization,
                                    size_t maxFacetGroupIndex) {
  _modelImages.SetFitsWriter(
      createWSCFitsWriter(entry, polarization, false, true, false).Writer());

  if (_settings.continuedRun) {
    readExistingModelImages(entry, polarization, maxFacetGroupIndex);
  } else {
    // Set model to zero: already done if this is YX of XY/YX imaging combi
    if (!(polarization == Polarization::YX &&
          _settings.polarizations.count(Polarization::XY) != 0)) {
      Image modelImage(_settings.trimmedImageWidth,
                       _settings.trimmedImageHeight, 0.0f);
      _modelImages.Store(modelImage.Data(), polarization,
                         entry.outputChannelIndex, false);
      if (Polarization::IsComplex(polarization))
        _modelImages.Store(modelImage.Data(), polarization,
                           entry.outputChannelIndex, true);
    }
  }
}

void WSClean::readExistingModelImages(const ImagingTableEntry& entry,
                                      PolarizationEnum polarization,
                                      size_t maxFacetGroupIndex) {
  // load image(s) from disk and store them in the model-image cache.
  for (size_t i = 0; i != entry.imageCount; ++i) {
    std::string prefix = ImageFilename::GetPrefix(
        _settings, polarization, entry.outputChannelIndex,
        entry.outputIntervalIndex, i == 1);

    const std::string suffix =
        (_settings.applyFacetBeam || !_settings.facetSolutionFiles.empty() ||
         griddingUsesATerms())
            ? "-model-pb.fits"
            : "-model.fits";
    aocommon::FitsReader reader(prefix + suffix);
    Logger::Info << "Reading " << reader.Filename() << "...\n";

    const bool resetGridder = overrideImageSettings(reader);

    // TODO check phase centre

    // FIXME: resetGridder in conjuncion with overrideImageSettings makes sure
    // that the image dimensions are set and passed to the _griddingTaskManager
    // only once. This probably can be simplified?
    if (resetGridder) {
      // Do not reset model column for a continuedRun
      if (!_settings.continuedRun) {
        resetModelColumns(entry);
      }
      _griddingTaskManager = GriddingTaskManager::Make(_settings);
      _griddingTaskManager->Start(getMaxNrMSProviders() *
                                  (maxFacetGroupIndex + 1));
    }

    if (!_imageWeightCache) {
      // The construction of the weight cache is delayed in prediction mode,
      // because only now the image size and scale is known.
      _imageWeightCache = createWeightCache();
      if (_settings.mfWeighting) initializeMFSImageWeights();
    }

    aocommon::FitsWriter writer(reader);
    _modelImages.SetFitsWriter(writer);

    Image buffer(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
    reader.Read(buffer.Data());
    for (size_t j = 0;
         j != _settings.trimmedImageWidth * _settings.trimmedImageHeight; ++j) {
      if (!std::isfinite(buffer[j]))
        throw std::runtime_error(
            "The input image contains non-finite values -- can't predict "
            "from an image with non-finite values");
    }
    _modelImages.Store(buffer.Data(), polarization, entry.outputChannelIndex,
                       i == 1);
  }
}

bool WSClean::overrideImageSettings(const aocommon::FitsReader& reader) {
  bool resetGridder = false;
  if (_settings.trimmedImageWidth == 0 && _settings.trimmedImageHeight == 0) {
    _settings.trimmedImageWidth = reader.ImageWidth();
    _settings.trimmedImageHeight = reader.ImageHeight();
    _settings.RecalculatePaddedDimensions();
    resetGridder = true;
  } else if (reader.ImageWidth() != _settings.trimmedImageWidth ||
             reader.ImageHeight() != _settings.trimmedImageHeight) {
    std::ostringstream msg;
    msg << "Inconsistent image size: dimensions of input image did not "
           "match, input: "
        << reader.ImageWidth() << " x " << reader.ImageHeight()
        << ", specified: " << _settings.trimmedImageWidth << " x "
        << _settings.trimmedImageHeight;
    throw std::runtime_error(msg.str());
  }

  if (reader.PixelSizeX() == 0.0 || reader.PixelSizeY() == 0.0)
    Logger::Warn
        << "Warning: input fits file misses the pixel size keywords.\n";
  else if (_settings.pixelScaleX == 0 && _settings.pixelScaleY == 0) {
    _settings.pixelScaleX = reader.PixelSizeX();
    _settings.pixelScaleY = reader.PixelSizeY();
    Logger::Debug << "Using pixel size of "
                  << Angle::ToNiceString(_settings.pixelScaleX) << " x "
                  << Angle::ToNiceString(_settings.pixelScaleY) << ".\n";
    resetGridder = true;
  }
  // Check if image corresponds with image dimensions of the settings
  // Here I require the pixel scale to be accurate enough so that the image
  // is at most 1/10th pixel larger/smaller.
  else if (std::fabs(reader.PixelSizeX() - _settings.pixelScaleX) *
                   _settings.trimmedImageWidth >
               0.1 * _settings.pixelScaleX ||
           std::fabs(reader.PixelSizeY() - _settings.pixelScaleY) *
                   _settings.trimmedImageHeight >
               0.1 * _settings.pixelScaleY) {
    std::ostringstream msg;
    msg << "Inconsistent pixel size: pixel size of input image did not "
           "match. Input: "
        << reader.PixelSizeX() << " x " << reader.PixelSizeY()
        << ", specified: " << _settings.pixelScaleX << " x "
        << _settings.pixelScaleY;
    throw std::runtime_error(msg.str());
  }
  if (_settings.pixelScaleX == 0.0 || _settings.pixelScaleY == 0.0) {
    throw std::runtime_error(
        "Could not determine proper pixel size. The input image did not "
        "provide proper pixel size values, and no or an invalid -scale was "
        "provided to WSClean");
  }
  return resetGridder;
}

void WSClean::predictGroup(const ImagingTable& groupTable) {
  const bool gridPolarizationsAtOnce =
      _settings.gridderType == GridderType::IDG &&
      _settings.polarizations.size() != 1;

  resetModelColumns(groupTable);
  _predictingWatch.Start();
  _griddingTaskManager->Start(getMaxNrMSProviders() *
                              (groupTable.MaxFacetGroupIndex() + 1));

  for (size_t groupIndex = 0; groupIndex != groupTable.IndependentGroupCount();
       ++groupIndex) {
    const ImagingTable independentGroup =
        groupTable.GetIndependentGroup(groupIndex);

    // Initialize the model images before entering the gridding loop. This is
    // necessary because in IDG mode, predicting Stokes I will require all
    // model images to have been initialized.
    for (size_t facetGroupIndex = 0;
         facetGroupIndex != independentGroup.FacetGroupCount();
         ++facetGroupIndex) {
      const ImagingTable facetGroup =
          independentGroup.GetFacetGroup(facetGroupIndex);
      // For facet-based prediction: facetGroup contains only a list of facets
      // from the same (full) image. The meta data for the full model image can
      // be inferred from the first entry in the facetGroup table
      _modelImages.Initialize(
          createWSCFitsWriter(facetGroup.Front(), false, true, false).Writer(),
          _settings.polarizations.size(), _settings.channelsOut, _facetCount,
          _settings.prefixName + "-model");

      readExistingModelImages(facetGroup.Front(),
                              facetGroup.Front().polarization,
                              groupTable.MaxFacetGroupIndex());
      partitionModelIntoFacets(facetGroup, true);
    }

    for (size_t facetGroupIndex = 0;
         facetGroupIndex != independentGroup.FacetGroupCount();
         ++facetGroupIndex) {
      const ImagingTable facetGroup =
          independentGroup.GetFacetGroup(facetGroupIndex);

      for (const auto& entry : facetGroup) {
        if (!gridPolarizationsAtOnce ||
            entry.polarization == *_settings.polarizations.begin()) {
          predict(entry);
        }
      }  // facets
    }    // facet groups of different polarizations
  }      // independent groups (channels)

  _griddingTaskManager->Finish();
  _predictingWatch.Pause();

  Logger::Info << "Inversion: " << _inversionWatch.ToString()
               << ", prediction: " << _predictingWatch.ToString()
               << ", cleaning: " << _deconvolutionWatch.ToString() << '\n';
}

void WSClean::initializeMSList(
    const ImagingTableEntry& entry,
    std::vector<std::unique_ptr<MSDataDescription>>& msList) {
  const PolarizationEnum pol = _settings.gridderType == GridderType::IDG
                                   ? getIdgPolarization()
                                   : entry.polarization;

  msList.clear();
  for (size_t msIndex = 0; msIndex != _settings.filenames.size(); ++msIndex) {
    for (size_t dataDescId = 0; dataDescId != _msBands[msIndex].DataDescCount();
         ++dataDescId) {
      MSSelection selection(_globalSelection);
      if (selection.SelectMsChannels(_msBands[msIndex], dataDescId, entry)) {
        std::unique_ptr<MSDataDescription> dataDescription;
        if (_settings.doReorder)
          dataDescription = MSDataDescription::ForPartitioned(
              _partitionedMSHandles[msIndex], selection,
              entry.msData[msIndex].bands[dataDescId].partIndex, pol,
              dataDescId, _settings.useMPI);
        else
          dataDescription = MSDataDescription::ForContiguous(
              _settings.filenames[msIndex], _settings.dataColumnName, selection,
              pol, dataDescId, _settings.useMPI);
        msList.emplace_back(std::move(dataDescription));
      }
    }
  }
}

void WSClean::resetModelColumns(const ImagingTable& groupTable) {
  if (groupTable.FacetCount() > 1) {
    for (const ImagingTable::Group& facetGroup : groupTable.FacetGroups()) {
      resetModelColumns(*facetGroup.front());
    }
  }
}

void WSClean::resetModelColumns(const ImagingTableEntry& entry) {
  std::vector<std::unique_ptr<MSDataDescription>> msList;
  initializeMSList(entry, msList);
  for (auto& ms : msList) {
    ms->GetProvider()->ResetModelColumn();
  }
}

void WSClean::runFirstInversions(ImagingTable& groupTable,
                                 std::unique_ptr<PrimaryBeam>& primaryBeam,
                                 bool requestPolarizationsAtOnce,
                                 bool parallelizePolarizations) {
  const size_t facetCount = groupTable.FacetCount();
  for (size_t facetIndex = 0; facetIndex < facetCount; ++facetIndex) {
    ImagingTable facetTable = groupTable.GetFacet(facetIndex);

    if (requestPolarizationsAtOnce) {
      for (ImagingTableEntry& entry : facetTable) {
        if (entry.polarization == *_settings.polarizations.begin())
          runSingleFirstInversion(entry, primaryBeam);
      }
    } else if (parallelizePolarizations) {
      for (ImagingTableEntry& entry : facetTable) {
        runSingleFirstInversion(entry, primaryBeam);
      }
    } else {
      bool hasMore;
      size_t sqIndex = 0;
      do {
        hasMore = false;
        // Run the inversion for one entry out of each squared group
        for (const ImagingTable::Group& sqGroup : facetTable.SquaredGroups()) {
          if (sqIndex < sqGroup.size()) {
            hasMore = true;
            runSingleFirstInversion(*sqGroup[sqIndex], primaryBeam);
          }
        }
        ++sqIndex;
        _griddingTaskManager->Finish();
      } while (hasMore);
    }
  }
  if (requestPolarizationsAtOnce) {
    _griddingTaskManager->Finish();
    groupTable.AssignGridDataFromPolarization(*_settings.polarizations.begin());
  } else if (parallelizePolarizations) {
    _griddingTaskManager->Finish();
  }
  stitchFacets(groupTable, _residualImages, _settings.isDirtySaved, false);
}

void WSClean::runSingleFirstInversion(
    ImagingTableEntry& entry, std::unique_ptr<PrimaryBeam>& primaryBeam) {
  const bool isLastPol =
      entry.polarization == *_settings.polarizations.rbegin();
  const bool doMakePSF = _settings.deconvolutionIterationCount > 0 ||
                         _settings.makePSF || _settings.makePSFOnly;

  if (isLastPol) {
    ImageFilename imageName =
        ImageFilename(entry.outputChannelIndex, entry.outputIntervalIndex);
    if (_settings.applyPrimaryBeam || _settings.applyFacetBeam) {
      std::vector<std::unique_ptr<MSDataDescription>> msList;
      initializeMSList(entry, msList);
      std::shared_ptr<ImageWeights> weights =
          initializeImageWeights(entry, msList);
      primaryBeam.reset(new PrimaryBeam(_settings));
      for (std::unique_ptr<MSDataDescription>& description : msList)
        primaryBeam->AddMS(std::move(description));
      primaryBeam->SetPhaseCentre(
          _observationInfo.phaseCentreRA, _observationInfo.phaseCentreDec,
          _observationInfo.shiftL, _observationInfo.shiftM);
      // Only generate beam images for facetIndex == 0 in facet group
      if (entry.facetIndex == 0) {
        primaryBeam->MakeBeamImages(imageName, entry, std::move(weights));
      }
    }
  }

  if (_settings.reuseDirty)
    loadExistingDirty(entry, !doMakePSF);
  else
    imageMain(entry, true, !doMakePSF);

  _isFirstInversion = false;
}

void WSClean::runMajorIterations(ImagingTable& groupTable,
                                 std::unique_ptr<PrimaryBeam>& primaryBeam,
                                 bool requestPolarizationsAtOnce,
                                 bool parallelizePolarizations) {
  std::unique_ptr<radler::WorkTable> deconvolution_table =
      groupTable.GetFacet(0).CreateDeconvolutionTable(
          _settings.deconvolutionChannelCount, _psfImages, _modelImages,
          _residualImages);

  _deconvolution.emplace(_settings.GetRadlerSettings(),
                         std::move(deconvolution_table),
                         minTheoreticalBeamSize(groupTable));

  if (_settings.deconvolutionIterationCount > 0) {
    // Start major cleaning loop
    _majorIterationNr = 1;
    bool reachedMajorThreshold = false;
    do {
      _deconvolutionWatch.Start();
      _deconvolution->Perform(reachedMajorThreshold, _majorIterationNr);
      _deconvolutionWatch.Pause();

      if (_majorIterationNr == 1 && _settings.deconvolutionMGain != 1.0 &&
          _settings.isFirstResidualSaved)
        writeFirstResidualImages(groupTable);
      const bool isFinished = !reachedMajorThreshold;
      if (isFinished) {
        writeModelImages(groupTable);
      }

      if (_settings.deconvolutionMGain != 1.0) {
        partitionModelIntoFacets(groupTable, false);
        if (requestPolarizationsAtOnce) {
          resetModelColumns(groupTable);
          _predictingWatch.Start();
          _griddingTaskManager->Start(getMaxNrMSProviders() *
                                      (groupTable.MaxFacetGroupIndex() + 1));
          // Iterate over polarizations, channels & facets
          for (const ImagingTableEntry& entry : groupTable) {
            // Only request one polarization for each facet/channel. The
            // gridder will grid all polarizations
            if (entry.polarization == *_settings.polarizations.begin())
              predict(entry);
          }
          _griddingTaskManager->Finish();

          _predictingWatch.Pause();
          _inversionWatch.Start();

          for (ImagingTableEntry& entry : groupTable) {
            if (entry.polarization == *_settings.polarizations.begin())
              imageMain(entry, false, false);
          }
          _griddingTaskManager->Finish();
          _inversionWatch.Pause();
        } else if (parallelizePolarizations) {
          resetModelColumns(groupTable);
          _predictingWatch.Start();
          _griddingTaskManager->Start(getMaxNrMSProviders() *
                                      (groupTable.MaxFacetGroupIndex() + 1));
          for (const ImagingTable::Group& sqGroup :
               groupTable.SquaredGroups()) {
            for (const ImagingTable::EntryPtr& entry : sqGroup) {
              predict(*entry);
            }
          }
          _griddingTaskManager->Finish();
          _predictingWatch.Pause();

          _inversionWatch.Start();
          for (const ImagingTable::Group& sqGroup :
               groupTable.SquaredGroups()) {
            for (const ImagingTable::EntryPtr& entry : sqGroup) {
              imageMain(*entry, false, false);
            }  // end of polarization & facets loop
          }    // end of joined channels loop
          _griddingTaskManager->Finish();
          _inversionWatch.Pause();
        } else {  // only parallelize channels
          resetModelColumns(groupTable);
          _predictingWatch.Start();
          _griddingTaskManager->Start(getMaxNrMSProviders() *
                                      (groupTable.MaxFacetGroupIndex() + 1));
          bool hasMore;
          size_t sqIndex = 0;
          do {
            hasMore = false;
            for (const ImagingTable::Group& sqGroup :
                 groupTable.SquaredGroups()) {
              if (sqIndex < sqGroup.size()) {
                hasMore = true;
                predict(*sqGroup[sqIndex]);
              }
            }
            ++sqIndex;
            _griddingTaskManager->Finish();
          } while (hasMore);
          _predictingWatch.Pause();

          _inversionWatch.Start();
          sqIndex = 0;
          do {
            hasMore = false;
            for (const ImagingTable::Group& sqGroup :
                 groupTable.SquaredGroups()) {
              if (sqIndex < sqGroup.size()) {
                hasMore = true;
                imageMain(*sqGroup[sqIndex], false, false);
              }
            }
            ++sqIndex;
            _griddingTaskManager->Finish();
          } while (hasMore);
          _inversionWatch.Pause();
        }
        stitchFacets(groupTable, _residualImages, false, false);
      }

      ++_majorIterationNr;
    } while (reachedMajorThreshold);

    --_majorIterationNr;
    Logger::Info << _majorIterationNr << " major iterations were performed.\n";
  }

  for (size_t facetGroupIndex = 0;
       facetGroupIndex != groupTable.FacetGroupCount(); ++facetGroupIndex) {
    const ImagingTable facetGroup = groupTable.GetFacetGroup(facetGroupIndex);
    saveRestoredImagesForGroup(facetGroup, primaryBeam);
  }

  if (_settings.saveSourceList) {
    std::unique_ptr<radler::WorkTable> deconvolution_table =
        groupTable.CreateDeconvolutionTable(_settings.deconvolutionChannelCount,
                                            _psfImages, _modelImages,
                                            _residualImages);
    ComponentListWriter componentListWriter(_settings,
                                            std::move(deconvolution_table));
    componentListWriter.SaveSourceList(*_deconvolution,
                                       _observationInfo.phaseCentreRA,
                                       _observationInfo.phaseCentreDec);
    if (usesBeam()) {
      componentListWriter.SavePbCorrectedSourceList(
          *_deconvolution, _observationInfo.phaseCentreRA,
          _observationInfo.phaseCentreDec);
    }
  }

  _deconvolution->FreeDeconvolutionAlgorithms();
}

MSSelection WSClean::selectInterval(MSSelection& fullSelection,
                                    size_t intervalIndex) {
  if (_settings.intervalsOut == 1)
    return fullSelection;
  else {
    size_t tS, tE;
    if (fullSelection.HasInterval()) {
      tS = fullSelection.IntervalStart();
      tE = fullSelection.IntervalEnd();
    } else {
      casacore::MeasurementSet ms(_settings.filenames[0]);
      Logger::Info << "Counting number of scans... ";
      Logger::Info.Flush();
      casacore::ScalarColumn<double> timeColumn(
          ms, casacore::MS::columnName(casacore::MSMainEnums::TIME));
      double time = timeColumn(0);
      size_t timestepIndex = 1;
      for (size_t row = 0; row != ms.nrow(); ++row) {
        if (time != timeColumn(row)) {
          ++timestepIndex;
          time = timeColumn(row);
        }
      }
      Logger::Info << "DONE (" << timestepIndex << ")\n";
      tS = 0;
      tE = timestepIndex;
      // Store the full interval in the selection, so that it doesn't need to
      // be determined again.
      fullSelection.SetInterval(tS, tE);
    }
    if (_settings.intervalsOut > tE - tS) {
      std::ostringstream str;
      str << "Invalid interval selection: " << _settings.intervalsOut
          << " intervals requested, but measurement set has only " << tE - tS
          << " intervals.";
      throw std::runtime_error(str.str());
    }
    MSSelection newSelection(fullSelection);
    newSelection.SetInterval(
        tS + (tE - tS) * intervalIndex / _settings.intervalsOut,
        tS + (tE - tS) * (intervalIndex + 1) / _settings.intervalsOut);
    return newSelection;
  }
}

void WSClean::saveUVImage(const Image& image, const ImagingTableEntry& entry,
                          bool isImaginary, const std::string& prefix) const {
  Image realUV(_settings.trimmedImageWidth, _settings.trimmedImageHeight),
      imagUV(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
  schaapcommon::fft::Resampler fft(
      _settings.trimmedImageWidth, _settings.trimmedImageHeight,
      _settings.trimmedImageWidth, _settings.trimmedImageHeight, 1);
  fft.SingleFT(image.Data(), realUV.Data(), imagUV.Data());
  // Factors of 2 involved: because of SingleFT()
  // (also one from the fact that normF excludes a factor of two?)
  realUV *=
      _infoPerChannel[entry.outputChannelIndex].normalizationFactor /
      sqrt(0.5 * _settings.trimmedImageWidth * _settings.trimmedImageHeight);
  imagUV *=
      _infoPerChannel[entry.outputChannelIndex].normalizationFactor /
      sqrt(0.5 * _settings.trimmedImageWidth * _settings.trimmedImageHeight);
  WSCFitsWriter writer(createWSCFitsWriter(entry, isImaginary, false, false));
  writer.WriteUV(prefix + "-real.fits", realUV.Data());
  writer.WriteUV(prefix + "-imag.fits", imagUV.Data());
}

void WSClean::stitchFacets(const ImagingTable& table,
                           CachedImageSet& imageCache, bool writeDirty,
                           bool isPSF) {
  if (_facetCount != 0) {
    Logger::Info << "Stitching facets onto full image...\n";
    // Allocate full image
    Image fullImage(_settings.trimmedImageWidth, _settings.trimmedImageHeight,
                    0.0f);
    // Initialize FacetImage with properties of stitched image, always
    // stitch facets for 1 spectral term.
    schaapcommon::facets::FacetImage facetImage(
        _settings.trimmedImageWidth, _settings.trimmedImageHeight, 1);
    for (size_t facetGroupIndex = 0; facetGroupIndex != table.FacetGroupCount();
         ++facetGroupIndex) {
      const ImagingTable stitchGroup = table.GetFacetGroup(facetGroupIndex);

      // The PSF is only once imaged for all polarizations
      if (!isPSF || stitchGroup.Front().polarization ==
                        *_settings.polarizations.begin()) {
        const size_t imageCount = stitchGroup.Front().imageCount;
        for (size_t imageIndex = 0; imageIndex != imageCount; ++imageIndex) {
          stitchSingleGroup(stitchGroup, imageIndex, imageCache, writeDirty,
                            isPSF, fullImage, facetImage,
                            table.MaxFacetGroupIndex());
        }
      }
    }
  }
}

void WSClean::stitchSingleGroup(const ImagingTable& facetGroup,
                                size_t imageIndex, CachedImageSet& imageCache,
                                bool writeDirty, bool isPSF, Image& fullImage,
                                schaapcommon::facets::FacetImage& facetImage,
                                size_t maxFacetGroupIndex) {
  const bool isImaginary = (imageIndex == 1);
  fullImage = 0.0f;
  for (const ImagingTableEntry& facetEntry : facetGroup) {
    facetImage.SetFacet(*facetEntry.facet, true);
    imageCache.LoadFacet(facetImage.Data(0), facetEntry.polarization,
                         facetEntry.outputChannelIndex, facetEntry.facetIndex,
                         facetEntry.facet, isImaginary);

    if (!isPSF &&
        (_settings.applyFacetBeam || !_settings.facetSolutionFiles.empty())) {
      const long double m =
          _msGridderMetaCache[facetEntry.index]->correctionSum /
          facetEntry.imageWeight;
      facetImage *= 1.0f / std::sqrt(m);
    }

    // TODO with our current stitching implementation, facets should always be
    // directly copied to the full image, not added. The facets should not
    // overlap though.
    facetImage.AddToImage({fullImage.Data()});
  }
  if (writeDirty) {
    initializeModelImages(facetGroup.Front(), facetGroup.Front().polarization,
                          maxFacetGroupIndex);
    _residualImages.SetFitsWriter(
        createWSCFitsWriter(facetGroup.Front(), false, false, true).Writer());
    WSCFitsWriter writer(
        createWSCFitsWriter(facetGroup.Front(), isImaginary, false, true));
    Logger::Info << "Writing dirty image...\n";
    writer.WriteImage("dirty.fits", fullImage.Data());
  }

  if (isPSF) {
    const ImagingTableEntry& entry = facetGroup.Front();
    processFullPSF(fullImage, entry);
  }

  const size_t channelIndex = facetGroup.Front().outputChannelIndex;
  const PolarizationEnum polarization = facetGroup.Front().polarization;
  imageCache.Store(fullImage.Data(), polarization, channelIndex, isImaginary);
}

void WSClean::makeImagingTable(size_t outputIntervalIndex) {
  std::set<aocommon::ChannelInfo> channelSet;
  _msBands.assign(_settings.filenames.size(), aocommon::MultiBandData());
  for (size_t i = 0; i != _settings.filenames.size(); ++i) {
    casacore::MeasurementSet ms(_settings.filenames[i]);
    _msBands[i] = aocommon::MultiBandData(ms);
    std::set<size_t> dataDescIds = _msBands[i].GetUsedDataDescIds(ms);
    if (dataDescIds.size() != _msBands[i].DataDescCount()) {
      Logger::Debug << dataDescIds.size() << "/" << _msBands[i].DataDescCount()
                    << " spws are used of " << _settings.filenames[i] << '\n';
    }

    // Apply user selection: remove unselected spws
    if (!_settings.spectralWindows.empty()) {
      for (std::set<size_t>::iterator d = dataDescIds.begin();
           d != dataDescIds.end();) {
        if (_settings.spectralWindows.find(_msBands[i].GetBandIndex(*d)) ==
            _settings.spectralWindows.end())
          d = dataDescIds.erase(d);
        else
          ++d;
      }
    }
    // accumulate channel info
    for (const size_t dataDescId : dataDescIds) {
      bool increasing = true;
      if (_msBands[i][dataDescId].ChannelCount() >= 2) {
        increasing = _msBands[i][dataDescId].Channel(1) >
                     _msBands[i][dataDescId].Channel(0);
      }
      channelSet.insert(_msBands[i][dataDescId].Channel(0));
      for (size_t ch = 1; ch != _msBands[i][dataDescId].ChannelCount(); ++ch) {
        bool chanIncreasing = _msBands[i][dataDescId].Channel(ch) >
                              _msBands[i][dataDescId].Channel(ch - 1);
        if (chanIncreasing != increasing)
          throw std::runtime_error(
              "Your measurement set has an incorrect frequency axis: the "
              "channels do neither only increase nor only decrease in "
              "frequency");
        if (_msBands[i][dataDescId].Channel(ch) ==
            _msBands[i][dataDescId].Channel(ch - 1))
          throw std::runtime_error(
              "Your measurement set has an incorrect frequency axis: two "
              "adjacent channels had the same frequency. Channels should "
              "either strictly increase or strictly decrease in frequency.");
        channelSet.insert(_msBands[i][dataDescId].Channel(ch));
      }
    }
  }
  if (channelSet.size() < _settings.channelsOut) {
    std::ostringstream str;
    str << "Parameter '-channels-out' was set to an invalid value: "
        << _settings.channelsOut
        << " output channels requested, but combined in all specified "
           "measurement sets, there are only "
        << channelSet.size() << " unique channels.";
    throw std::runtime_error(str.str());
  }
  std::vector<aocommon::ChannelInfo> inputChannelFrequencies(channelSet.begin(),
                                                             channelSet.end());
  Logger::Debug << "Total nr of channels found in measurement sets: "
                << inputChannelFrequencies.size() << '\n';

  _imagingTable.Clear();

  ImagingTableEntry templateEntry;
  templateEntry.joinedGroupIndex = 0;
  templateEntry.squaredDeconvolutionIndex = 0;

  // for(size_t interval=0; interval!=_settings.intervalsOut; ++interval)
  //{
  for (size_t outChannelIndex = 0; outChannelIndex != _settings.channelsOut;
       ++outChannelIndex) {
    makeImagingTableEntry(inputChannelFrequencies, outputIntervalIndex,
                          outChannelIndex, templateEntry);
    templateEntry.outputChannelIndex = outChannelIndex;

    if (_settings.joinedFrequencyDeconvolution) {
      templateEntry.joinedGroupIndex = 0;
    }
    addPolarizationsToImagingTable(templateEntry);
  }
  //}
  _imagingTable.Update();
  _imagingTable.Print();
}

void WSClean::makeImagingTableEntry(
    const std::vector<aocommon::ChannelInfo>& channels, size_t outIntervalIndex,
    size_t outChannelIndex, ImagingTableEntry& entry) {
  size_t startCh, endCh;
  if (_settings.endChannel != 0) {
    if (_settings.endChannel > channels.size())
      throw std::runtime_error(
          "Bad channel selection -- more channels selected than available");
    startCh = _settings.startChannel;
    endCh = _settings.endChannel;
  } else {
    startCh = 0;
    endCh = channels.size();
  }
  std::vector<aocommon::ChannelInfo> groupChannels(channels.begin() + startCh,
                                                   channels.begin() + endCh);

  if (_settings.divideChannelFrequencies.empty()) {
    makeImagingTableEntryChannelSettings(groupChannels, outIntervalIndex,
                                         outChannelIndex, _settings.channelsOut,
                                         entry);
  } else {
    // We need to separately divide the channels into groups as specified and
    // call the freq division for the group corresponding with the
    // outChannelIndex.
    const size_t nSplits = _settings.divideChannelFrequencies.size();
    for (size_t i = 0; i != nSplits + 1; ++i) {
      const size_t outChannelStart = _settings.channelsOut * i / (nSplits + 1);
      const size_t outChannelEnd =
          _settings.channelsOut * (i + 1) / (nSplits + 1);
      if (outChannelIndex >= outChannelStart &&
          outChannelIndex < outChannelEnd) {
        double splitFreqLow =
            (i == 0) ? 0.0 : _settings.divideChannelFrequencies[i - 1];
        double splitFreqHigh = (i == nSplits)
                                   ? std::numeric_limits<double>::max()
                                   : _settings.divideChannelFrequencies[i];
        std::vector<aocommon::ChannelInfo> splittedChannels;
        for (const aocommon::ChannelInfo& channel : groupChannels) {
          if (channel.Frequency() >= splitFreqLow &&
              channel.Frequency() < splitFreqHigh)
            splittedChannels.emplace_back(channel);
        }
        size_t nOutChannels = outChannelEnd - outChannelStart;
        makeImagingTableEntryChannelSettings(splittedChannels, outIntervalIndex,
                                             outChannelIndex - outChannelStart,
                                             nOutChannels, entry);
      }
    }
  }

  if (_settings.spectralCorrection.empty())
    entry.siCorrection = 1.0;
  else {
    double bandwidthCentre =
        0.5 * (channels.front().Frequency() + channels.back().Frequency());
    double chCentralFrequency =
        0.5 * (entry.lowestFrequency + entry.highestFrequency);
    double chFlux = schaapcommon::fitters::NonLinearPowerLawFitter::Evaluate(
        chCentralFrequency, _settings.spectralCorrection,
        _settings.spectralCorrectionFrequency);
    double midFlux = schaapcommon::fitters::NonLinearPowerLawFitter::Evaluate(
        bandwidthCentre, _settings.spectralCorrection,
        _settings.spectralCorrectionFrequency);
    entry.siCorrection = midFlux / chFlux;
    if (outChannelIndex == 0)
      Logger::Debug << "SI correction for first channel: " << entry.siCorrection
                    << '\n';
    if (outChannelIndex + 1 == _settings.channelsOut)
      Logger::Debug << "SI correction for last channel: " << entry.siCorrection
                    << '\n';
  }

  entry.msData.resize(_settings.filenames.size());
  for (size_t msIndex = 0; msIndex != _settings.filenames.size(); ++msIndex) {
    entry.msData[msIndex].bands.resize(_msBands[msIndex].DataDescCount());
  }
}

void WSClean::makeImagingTableEntryChannelSettings(
    const std::vector<aocommon::ChannelInfo>& channels, size_t outIntervalIndex,
    size_t outChannelIndex, size_t nOutChannels, ImagingTableEntry& entry) {
  size_t chLowIndex, chHighIndex;
  if (_settings.divideChannelsByGaps) {
    std::multimap<double, size_t> gaps;
    for (size_t i = 1; i != channels.size(); ++i) {
      double left = channels[i - 1].Frequency();
      double right = channels[i].Frequency();
      gaps.emplace(right - left, i);
    }
    std::vector<size_t> orderedGaps;
    auto iter = gaps.rbegin();
    for (size_t i = 0; i != nOutChannels - 1; ++i) {
      if (iter == gaps.rend())
        throw std::runtime_error(
            "Channel gap division leads to invalid selection");
      orderedGaps.push_back(iter->second);
      ++iter;
    }
    std::sort(orderedGaps.begin(), orderedGaps.end());
    if (outChannelIndex == 0)
      chLowIndex = 0;
    else
      chLowIndex = orderedGaps[outChannelIndex - 1];
    if (outChannelIndex + 1 == nOutChannels)
      chHighIndex = channels.size() - 1;
    else
      chHighIndex = orderedGaps[outChannelIndex] - 1;
  } else {
    chLowIndex = outChannelIndex * channels.size() / nOutChannels;
    chHighIndex = (outChannelIndex + 1) * channels.size() / nOutChannels - 1;
    if (chLowIndex == chHighIndex + 1)
      throw std::runtime_error(
          "Too many output channels requested: output channel " +
          std::to_string(outChannelIndex) +
          " would be empty. Number of output channels requested: " +
          std::to_string(_settings.channelsOut) +
          ". Number of channels in the measurement set(s) available (after "
          "applying channel range selections and splits): " +
          std::to_string(channels.size()));
  }
  if (channels[chLowIndex].Frequency() > channels[chHighIndex].Frequency())
    std::swap(chLowIndex, chHighIndex);
  entry.inputChannelCount = chHighIndex + 1 - chLowIndex;
  entry.lowestFrequency = channels[chLowIndex].Frequency();
  entry.highestFrequency = channels[chHighIndex].Frequency();
  entry.bandStartFrequency =
      entry.lowestFrequency - channels[chLowIndex].Width() * 0.5;
  entry.bandEndFrequency =
      entry.highestFrequency + channels[chHighIndex].Width() * 0.5;
  entry.outputIntervalIndex = outIntervalIndex;
}

void WSClean::addPolarizationsToImagingTable(ImagingTableEntry& templateEntry) {
  for (PolarizationEnum p : _settings.polarizations) {
    templateEntry.polarization = p;
    if (p == Polarization::XY)
      templateEntry.imageCount = 2;
    else if (p == Polarization::YX)
      templateEntry.imageCount = 0;
    else
      templateEntry.imageCount = 1;

    addFacetsToImagingTable(templateEntry);

    if (!_settings.joinedPolarizationDeconvolution) {
      ++templateEntry.joinedGroupIndex;
      ++templateEntry.squaredDeconvolutionIndex;
    }
  }

  if (_settings.joinedPolarizationDeconvolution) {
    ++templateEntry.joinedGroupIndex;
    ++templateEntry.squaredDeconvolutionIndex;
  }
}

void WSClean::addFacetsToImagingTable(ImagingTableEntry& templateEntry) {
  // Create a single entry (with facetIndex == 0) when facets are not used.
  const size_t facet_entry_count = std::max(_facetCount, std::size_t(1));
  for (size_t f = 0; f != facet_entry_count; ++f) {
    auto entry = std::make_unique<ImagingTableEntry>(templateEntry);
    entry->facetIndex = f;
    entry->facet.reset();  // updateFacetsInImagingTable will set the facet.
    _imagingTable.AddEntry(std::move(entry));
  }
  ++templateEntry.facetGroupIndex;
}

void WSClean::updateFacetsInImagingTable(
    const std::vector<std::shared_ptr<schaapcommon::facets::Facet>>& facets) {
  for (ImagingTableEntry& entry : _imagingTable) {
    assert(entry.facetIndex < facets.size());
    entry.facet = facets[entry.facetIndex];

    // Calculate phase center delta for entry
    entry.centreShiftX = entry.facet->GetUntrimmedBoundingBox().Centre().x -
                         _settings.trimmedImageWidth / 2;
    entry.centreShiftY = entry.facet->GetUntrimmedBoundingBox().Centre().y -
                         _settings.trimmedImageHeight / 2;
  }
}

WSCFitsWriter WSClean::createWSCFitsWriter(const ImagingTableEntry& entry,
                                           bool isImaginary, bool isModel,
                                           bool isFullImage) const {
  ObservationInfo observationInfo = _observationInfo;
  if (!isFullImage) {
    applyFacetPhaseShift(entry, observationInfo);
  }

  return WSCFitsWriter(entry, isImaginary, _settings, _deconvolution,
                       observationInfo, _majorIterationNr, _commandLine,
                       _infoPerChannel[entry.outputChannelIndex], isModel,
                       _lastStartTime);
}

WSCFitsWriter WSClean::createWSCFitsWriter(const ImagingTableEntry& entry,
                                           PolarizationEnum polarization,
                                           bool isImaginary, bool isModel,
                                           bool isFullImage) const {
  ObservationInfo observationInfo = _observationInfo;
  if (!isFullImage) {
    applyFacetPhaseShift(entry, observationInfo);
  }

  return WSCFitsWriter(entry, polarization, isImaginary, _settings,
                       _deconvolution, observationInfo, _majorIterationNr,
                       _commandLine, _infoPerChannel[entry.outputChannelIndex],
                       isModel, _lastStartTime);
}

void WSClean::correctImagesH5(aocommon::FitsWriter& writer,
                              const ImagingTable& table,
                              const ImageFilename& imageName,
                              const std::string& filenameKind) const {
  const PolarizationEnum pol = table.Front().polarization;

  if (pol == Polarization::StokesI || pol == Polarization::XX ||
      pol == Polarization::YY) {
    ImageFilename iName(imageName);
    iName.SetPolarization(pol);
    std::string prefix;
    if (filenameKind == "psf")
      prefix = iName.GetPSFPrefix(_settings);
    else
      prefix = iName.GetPrefix(_settings);

    aocommon::FitsReader reader(prefix + "-" + filenameKind + ".fits");
    Image image(reader.ImageWidth(), reader.ImageHeight());
    reader.Read(image.Data());

    schaapcommon::facets::FacetImage facetImage(
        _settings.trimmedImageWidth, _settings.trimmedImageHeight, 1);
    std::vector<float*> imagePtr{image.Data()};
    for (const ImagingTableEntry& entry : table) {
      facetImage.SetFacet(*entry.facet, true);
      // Requires std::map::at in order to comply with constness of member
      // function
      const long double m = _msGridderMetaCache.at(entry.index)->correctionSum /
                            entry.imageWeight;
      facetImage.MultiplyImageInsideFacet(imagePtr, 1.0f / std::sqrt(m));
    }

    // Always write to -pb.fits
    writer.Write(prefix + "-" + filenameKind + "-pb.fits", image.Data());
  } else {
    throw std::runtime_error(
        "H5 correction is requested, but this is not supported "
        "when imaging a single polarization that is not Stokes I, XX, or YY.");
  }
}
