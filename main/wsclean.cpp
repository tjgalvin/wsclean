#include "wsclean.h"

#include "../math/imageoperations.h"

#include "../structures/primarybeam.h"

#include "../gridding/directmsgridder.h"
#include "../gridding/measurementsetgridder.h"

#include "../io/facetreader.h"
#include "../io/imagefilename.h"
#include "../io/imageweightcache.h"
#include "../io/logger.h"
#include "../io/wscfitswriter.h"

#include "../scheduling/griddingtaskmanager.h"

#include "../units/angle.h"

#include "../system/application.h"

#include "../structures/image.h"
#include "../structures/imageweights.h"
#include "../structures/msselection.h"

#include "../deconvolution/deconvolutionalgorithm.h"
#include "../deconvolution/imageset.h"

#include "../idg/idgmsgridder.h"

#include "../math/fftresampler.h"
#include "../math/modelrenderer.h"
#include "../math/nlplfitter.h"

#include "../model/model.h"

#include "../msproviders/contiguousms.h"
#include "../msproviders/msdatadescription.h"

#include "progressbar.h"

#include <schaapcommon/facets/facetimage.h>

#include <aocommon/fits/fitswriter.h>
#include <aocommon/uvector.h>
#include <aocommon/parallelfor.h>

#include <schaapcommon/facets/facetimage.h>

#include <iostream>
#include <memory>

std::string commandLine;

WSClean::WSClean()
    : _globalSelection(),
      _commandLine(),
      _inversionWatch(false),
      _predictingWatch(false),
      _deconvolutionWatch(false),
      _isFirstInversion(true),
      _majorIterationNr(0),
      _deconvolution(_settings),
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
  size_t channelIndex = entry.outputChannelIndex;
  Settings modifiedSettings(_settings);
  modifiedSettings.prefixName = _settings.reusePsfPrefix;
  std::string name;
  if (isPSF)
    name = ImageFilename::GetPSFPrefix(modifiedSettings, channelIndex,
                                       entry.outputIntervalIndex) +
           "-psf.fits";
  else
    name = ImageFilename::GetPrefix(modifiedSettings, entry.polarization,
                                    channelIndex, entry.outputIntervalIndex,
                                    false) +
           "-dirty.fits";
  FitsReader reader(name);
  if (reader.ImageWidth() != _settings.trimmedImageWidth ||
      reader.ImageHeight() != _settings.trimmedImageHeight)
    throw std::runtime_error(
        "Image width and height of reused PSF don't match with given settings");
  ImageF psfImage(reader.ImageWidth(), reader.ImageHeight());
  reader.Read(psfImage.data());

  GriddingResult result;
  result.imageRealResult = std::move(psfImage);
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
  imagePSFCallback(entry, result);
}

void WSClean::loadExistingDirty(ImagingTableEntry& entry, bool updateBeamInfo) {
  Logger::Info << "Loading existing dirty image from disk...\n";
  GriddingResult result = loadExistingImage(entry, false);
  imageMainCallback(entry, result, updateBeamInfo, true);
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
  applyFacetPhaseShift(entry, task.observationInfo);
  initializeMSList(entry, task.msList);
  task.imageWeights = initializeImageWeights(entry, task.msList);

  _griddingTaskManager->Run(std::move(task),
                            [this, &entry](GriddingResult& result) {
                              imagePSFCallback(entry, result);
                            });
}

void WSClean::imagePSFCallback(ImagingTableEntry& entry,
                               GriddingResult& result) {
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

  if (_facets.empty()) processFullPSF(result.imageRealResult, entry);

  _lastStartTime = result.startTime;
  _msGridderMetaCache[entry.index] = std::move(result.cache);

  _psfImages.SetFitsWriter(
      createWSCFitsWriter(entry, false, false, false).Writer());
  _psfImages.StoreFacet(result.imageRealResult.data(),
                        *_settings.polarizations.begin(), channelIndex,
                        entry.facetIndex, entry.facet, false);

  bool isLastPol = entry.polarization == *_settings.polarizations.rbegin();
  if (isLastPol &&
      (_settings.gridWithBeam || !_settings.atermConfigFilename.empty())) {
    Logger::Info << "Writing IDG beam image...\n";
    ImageFilename imageName(entry.outputChannelIndex,
                            entry.outputIntervalIndex);
    IdgMsGridder::SaveBeamImage(
        entry, imageName, _settings, _observationInfo.phaseCentreRA,
        _observationInfo.phaseCentreDec, _observationInfo.shiftL,
        _observationInfo.shiftM, *_msGridderMetaCache[entry.index]);
  }

  _isFirstInversion = false;
}

void WSClean::processFullPSF(ImageF& image, const ImagingTableEntry& entry) {
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

  DeconvolutionAlgorithm::RemoveNaNsInPSF(
      image.data(), _settings.trimmedImageWidth, _settings.trimmedImageHeight);

  double minPixelScale = std::min(_settings.pixelScaleX, _settings.pixelScaleY);
  double initialFitSize =
      std::max(_infoPerChannel[channelIndex].beamSizeEstimate, minPixelScale);
  double bMaj, bMin, bPA, bTheoretical;
  ImageOperations::DetermineBeamSize(_settings, bMaj, bMin, bPA, bTheoretical,
                                     image.data(), initialFitSize);
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
  fitsFile.WritePSF(name, image.data());
  Logger::Info << "DONE\n";
}

void WSClean::imageMain(ImagingTableEntry& entry, bool isFirstInversion,
                        bool updateBeamInfo) {
  Logger::Info.Flush();
  Logger::Info << " == Constructing image ==\n";

  GriddingTask task;
  task.operation = GriddingTask::Invert;
  task.imagePSF = false;
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

  result.imageRealResult *=
      _infoPerChannel[joinedChannelIndex].psfNormalizationFactor *
      entry.siCorrection;
  storeAndCombineXYandYX(_residualImages, joinedChannelIndex, entry, false,
                         result.imageRealResult.data());
  if (aocommon::Polarization::IsComplex(entry.polarization)) {
    result.imageImaginaryResult *=
        _infoPerChannel[joinedChannelIndex].psfNormalizationFactor *
        entry.siCorrection;
    storeAndCombineXYandYX(_residualImages, joinedChannelIndex, entry, true,
                           result.imageImaginaryResult.data());
  }
  _msGridderMetaCache[entry.index] = std::move(result.cache);
  entry.imageWeight = result.imageWeight;
  entry.normalizationFactor = result.normalizationFactor;
  _infoPerChannel[entry.outputChannelIndex].weight = result.imageWeight;
  _infoPerChannel[entry.outputChannelIndex].normalizationFactor =
      result.normalizationFactor;

  _lastStartTime = result.startTime;

  // If no PSF is made, also set the beam size. If the PSF was made, these would
  // already be set after imaging the PSF.
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

  // If !_facets.empty(), these actions are performed in stitchFacets
  if (isInitialInversion && _facets.empty()) {
    initializeModelImages(entry);

    _residualImages.SetFitsWriter(
        createWSCFitsWriter(entry, false, false, false).Writer());
    // If !_facets.empty(), dirty image is saved in stitchFacets
    if (_settings.isDirtySaved) {
      for (size_t imageIndex = 0; imageIndex != entry.imageCount;
           ++imageIndex) {
        const bool isImaginary = (imageIndex == 1);
        WSCFitsWriter writer(
            createWSCFitsWriter(entry, isImaginary, false, false));
        ImageF dirtyImage(_settings.trimmedImageWidth,
                          _settings.trimmedImageHeight);
        _residualImages.Load(dirtyImage.data(), entry.polarization,
                             entry.outputChannelIndex, isImaginary);
        Logger::Info << "Writing dirty image...\n";
        writer.WriteImage("dirty.fits", dirtyImage.data());
      }
    }
  }
}

void WSClean::storeAndCombineXYandYX(CachedImageSet& dest,
                                     size_t joinedChannelIndex,
                                     const ImagingTableEntry& entry,
                                     bool isImaginary, const float* image) {
  if (entry.polarization == aocommon::Polarization::YX &&
      _settings.polarizations.count(aocommon::Polarization::XY) != 0) {
    Logger::Info << "Adding XY and YX together...\n";
    ImageF xyImage;
    if (entry.facet) {
      // Trimmed facet size
      xyImage = ImageF(entry.facet->GetTrimmedBoundingBox().Width(),
                       entry.facet->GetTrimmedBoundingBox().Height());
    } else {
      // Full image size
      xyImage =
          ImageF(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
    }

    dest.LoadFacet(xyImage.data(), aocommon::Polarization::XY,
                   joinedChannelIndex, entry.facetIndex, entry.facet,
                   isImaginary);
    if (isImaginary) {
      for (size_t i = 0; i != xyImage.size(); ++i)
        xyImage[i] = (xyImage[i] - image[i]) * 0.5;
    } else {
      for (size_t i = 0; i != xyImage.size(); ++i)
        xyImage[i] = (xyImage[i] + image[i]) * 0.5;
    }
    dest.StoreFacet(xyImage.data(), aocommon::Polarization::XY,
                    joinedChannelIndex, entry.facetIndex, entry.facet,
                    isImaginary);
  } else {
    dest.StoreFacet(image, entry.polarization, joinedChannelIndex,
                    entry.facetIndex, entry.facet, isImaginary);
  }
}

void WSClean::predict(const ImagingTableEntry& entry) {
  Logger::Info.Flush();
  Logger::Info << " == Converting model image to visibilities ==\n";
  ImageF modelImageReal(_settings.trimmedImageWidth,
                        _settings.trimmedImageHeight),
      modelImageImaginary;

  if (entry.facet)
    throw std::runtime_error("Predicting facets is not implemented");

  if (entry.polarization == aocommon::Polarization::YX) {
    _modelImages.Load(modelImageReal.data(), aocommon::Polarization::XY,
                      entry.outputChannelIndex, false);
    modelImageImaginary =
        ImageF(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
    _modelImages.Load(modelImageImaginary.data(), aocommon::Polarization::XY,
                      entry.outputChannelIndex, true);
    const size_t size =
        _settings.trimmedImageWidth * _settings.trimmedImageHeight;
    for (size_t i = 0; i != size; ++i)
      modelImageImaginary[i] = -modelImageImaginary[i];
  } else {
    _modelImages.Load(modelImageReal.data(), entry.polarization,
                      entry.outputChannelIndex, false);
    if (aocommon::Polarization::IsComplex(entry.polarization)) {
      modelImageImaginary =
          ImageF(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
      _modelImages.Load(modelImageImaginary.data(), entry.polarization,
                        entry.outputChannelIndex, true);
    }
  }

  GriddingTask task;
  task.operation = GriddingTask::Predict;
  task.polarization = entry.polarization;
  task.addToModel = false;
  task.cache = std::move(_msGridderMetaCache[entry.index]);
  task.verbose = false;
  task.storeImagingWeights = false;
  task.modelImageReal = std::move(modelImageReal);
  task.modelImageImaginary = std::move(modelImageImaginary);
  initializeMSList(entry, task.msList);
  task.imageWeights = initializeImageWeights(entry, task.msList);
  task.observationInfo = _observationInfo;
  task.facet = entry.facet;
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
          partSelection.SetBandId(dataDescId);
          const bool hasSelection =
              selectChannels(partSelection, msIndex, dataDescId, entry);
          if (hasSelection) {
            aocommon::PolarizationEnum pol =
                _settings.useIDG ? aocommon::Polarization::Instrumental
                                 : entry.polarization;
            PartitionedMS msProvider(_partitionedMSHandles[msIndex],
                                     ms.bands[dataDescId].partIndex, pol,
                                     dataDescId);
            weights->Grid(msProvider, partSelection);
          }
        }
      }
    }
  } else {
    for (size_t i = 0; i != _settings.filenames.size(); ++i) {
      for (size_t d = 0; d != _msBands[i].DataDescCount(); ++d) {
        aocommon::PolarizationEnum pol =
            _settings.useIDG ? aocommon::Polarization::Instrumental
                             : *_settings.polarizations.begin();
        ContiguousMS msProvider(_settings.filenames[i],
                                _settings.dataColumnName, _globalSelection, pol,
                                d);
        weights->Grid(msProvider, _globalSelection);
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
  _partitionedMSHandles.resize(_settings.filenames.size());
  bool useModel = _settings.deconvolutionMGain != 1.0 || isPredictMode ||
                  _settings.subtractModel || _settings.continuedRun;
  bool initialModelRequired = _settings.subtractModel || _settings.continuedRun;

  if (_settings.parallelReordering != 1) Logger::Info << "Reordering...\n";

  aocommon::ParallelFor<size_t> loop(_settings.parallelReordering);
  loop.Run(0, _settings.filenames.size(), [&](size_t i, size_t) {
    std::vector<PartitionedMS::ChannelRange> channels;
    std::map<aocommon::PolarizationEnum, size_t> nextIndex;
    for (size_t sqIndex = 0; sqIndex != _imagingTable.SquaredGroupCount();
         ++sqIndex) {
      ImagingTable sqGroup = _imagingTable.GetSquaredGroup(sqIndex);
      for (size_t fgIndex = 0; fgIndex != sqGroup.FacetGroupCount();
           ++fgIndex) {
        ImagingTable facetGroup = sqGroup.GetFacetGroup(fgIndex);
        // The band information is determined from the first facet in the group.
        // After this, all facet entries inside the group are updated.
        const ImagingTableEntry& entry = facetGroup.Front();
        for (size_t d = 0; d != _msBands[i].DataDescCount(); ++d) {
          MSSelection selection(_globalSelection);
          if (selectChannels(selection, i, d, entry)) {
            if (entry.polarization == *_settings.polarizations.begin()) {
              PartitionedMS::ChannelRange r;
              r.dataDescId = d;
              r.start = selection.ChannelRangeStart();
              r.end = selection.ChannelRangeEnd();
              channels.push_back(r);
            }
            for (ImagingTableEntry& facetEntry : facetGroup) {
              facetEntry.msData[i].bands[d].partIndex =
                  nextIndex[entry.polarization];
            }
            ++nextIndex[entry.polarization];
          }
        }
      }
    }

    PartitionedMS::Handle partMS = PartitionedMS::Partition(
        _settings.filenames[i], channels, _globalSelection,
        _settings.dataColumnName, useModel, initialModelRequired, _settings);
    std::lock_guard<std::mutex> lock(mutex);
    _partitionedMSHandles[i] = std::move(partMS);
    if (_settings.parallelReordering != 1)
      Logger::Info << "Finished reordering " << _settings.filenames[i] << " ["
                   << i << "]\n";
  });
}

void WSClean::RunClean() {
  _observationInfo = getObservationInfo();
  _facets = FacetReader::ReadFacets(_settings.facetRegionFilename);
  for (schaapcommon::facets::Facet& facet : _facets) {
    const size_t alignment = 4;
    facet.CalculatePixels(
        _observationInfo.phaseCentreRA, _observationInfo.phaseCentreDec,
        _settings.pixelScaleX, _settings.pixelScaleY,
        _settings.trimmedImageWidth, _settings.trimmedImageHeight,
        _observationInfo.shiftL, _observationInfo.shiftM,
        _settings.imagePadding, alignment, _settings.useIDG);
  }

  _globalSelection = _settings.GetMSSelection();
  MSSelection fullSelection = _globalSelection;

  for (size_t intervalIndex = 0; intervalIndex != _settings.intervalsOut;
       ++intervalIndex) {
    makeImagingTable(intervalIndex);

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
      for (std::set<aocommon::PolarizationEnum>::const_iterator pol =
               _settings.polarizations.begin();
           pol != _settings.polarizations.end(); ++pol) {
        bool psfWasMade = (_settings.deconvolutionIterationCount > 0 ||
                           _settings.makePSF || _settings.makePSFOnly) &&
                          pol == _settings.polarizations.begin();

        if (psfWasMade) {
          ImageOperations::MakeMFSImage(_settings, _infoPerChannel, _infoForMFS,
                                        "psf.fits", intervalIndex, *pol, false,
                                        true);
          if (_settings.savePsfPb)
            ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                          _infoForMFS, "psf-pb.fits",
                                          intervalIndex, *pol, false, true);
        }

        if (!(*pol == aocommon::Polarization::YX &&
              _settings.polarizations.count(aocommon::Polarization::XY) != 0) &&
            !_settings.makePSFOnly) {
          if (_settings.isDirtySaved)
            ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                          _infoForMFS, "dirty.fits",
                                          intervalIndex, *pol, false);
          if (_settings.deconvolutionIterationCount == 0) {
            ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                          _infoForMFS, "image.fits",
                                          intervalIndex, *pol, false);
            if (_settings.applyPrimaryBeam ||
                (_settings.gridWithBeam ||
                 !_settings.atermConfigFilename.empty()))
              ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                            _infoForMFS, "image-pb.fits",
                                            intervalIndex, *pol, false);
          } else {
            ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                          _infoForMFS, "residual.fits",
                                          intervalIndex, *pol, false);
            ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                          _infoForMFS, "model.fits",
                                          intervalIndex, *pol, false);
            ImageOperations::RenderMFSImage(_settings, _infoForMFS,
                                            intervalIndex, *pol, false, false);
            if (_settings.applyPrimaryBeam ||
                (_settings.gridWithBeam ||
                 !_settings.atermConfigFilename.empty())) {
              ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                            _infoForMFS, "residual-pb.fits",
                                            intervalIndex, *pol, false);
              ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                            _infoForMFS, "model-pb.fits",
                                            intervalIndex, *pol, false);
              ImageOperations::RenderMFSImage(_settings, _infoForMFS,
                                              intervalIndex, *pol, false, true);
            }
          }
          if (aocommon::Polarization::IsComplex(*pol)) {
            if (_settings.isDirtySaved)
              ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                            _infoForMFS, "dirty.fits",
                                            intervalIndex, *pol, true);
            if (_settings.deconvolutionIterationCount == 0) {
              ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                            _infoForMFS, "image.fits",
                                            intervalIndex, *pol, true);
              if (_settings.applyPrimaryBeam ||
                  (_settings.gridWithBeam ||
                   !_settings.atermConfigFilename.empty()))
                ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                              _infoForMFS, "image-pb.fits",
                                              intervalIndex, *pol, true);
            } else {
              ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                            _infoForMFS, "residual.fits",
                                            intervalIndex, *pol, true);
              ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                            _infoForMFS, "model.fits",
                                            intervalIndex, *pol, true);
              ImageOperations::RenderMFSImage(_settings, _infoForMFS,
                                              intervalIndex, *pol, true, false);
              if (_settings.applyPrimaryBeam ||
                  (_settings.gridWithBeam ||
                   !_settings.atermConfigFilename.empty())) {
                ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                              _infoForMFS, "residual-pb.fits",
                                              intervalIndex, *pol, true);
                ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                              _infoForMFS, "model-pb.fits",
                                              intervalIndex, *pol, true);
                ImageOperations::RenderMFSImage(
                    _settings, _infoForMFS, intervalIndex, *pol, true, true);
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
      _settings.useWeightsAsTaper));
  cache->SetTaperInfo(
      _settings.gaussianTaperBeamSize, _settings.tukeyTaperInLambda,
      _settings.tukeyInnerTaperInLambda, _settings.edgeTaperInLambda,
      _settings.edgeTukeyTaperInLambda);
  return std::move(cache);
}

void WSClean::RunPredict() {
  _observationInfo = getObservationInfo();

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

    _griddingTaskManager = GriddingTaskManager::Make(_settings);

    for (const ImagingTable::Group& group : _imagingTable.SquaredGroups()) {
      predictGroup(group);
    }

    // Needs to be destructed before image allocator, or image allocator will
    // report error caused by leaked memory
    _griddingTaskManager.reset();
  }
}

bool WSClean::selectChannels(MSSelection& selection, size_t msIndex,
                             size_t dataDescId,
                             const ImagingTableEntry& entry) {
  const BandData& band = _msBands[msIndex][dataDescId];
  double firstCh = band.ChannelFrequency(0),
         lastCh = band.ChannelFrequency(band.ChannelCount() - 1);
  // Some mses have decreasing (i.e. reversed) channel frequencies in them
  bool isReversed = false;
  if (firstCh > lastCh) {
    std::swap(firstCh, lastCh);
    isReversed = true;
    Logger::Debug << "Warning: MS has reversed channel frequencies.\n";
  }
  if (band.ChannelCount() != 0 && entry.lowestFrequency <= lastCh &&
      entry.highestFrequency >= firstCh) {
    size_t newStart, newEnd;
    if (isReversed) {
      BandData::const_reverse_iterator lowPtr, highPtr;
      lowPtr =
          std::lower_bound(band.rbegin(), band.rend(), entry.lowestFrequency);
      highPtr = std::lower_bound(lowPtr, band.rend(), entry.highestFrequency);

      if (highPtr == band.rend()) --highPtr;
      newStart = band.ChannelCount() - 1 - (highPtr - band.rbegin());
      newEnd = band.ChannelCount() - (lowPtr - band.rbegin());
    } else {
      const double *lowPtr, *highPtr;
      lowPtr =
          std::lower_bound(band.begin(), band.end(), entry.lowestFrequency);
      highPtr = std::lower_bound(lowPtr, band.end(), entry.highestFrequency);

      if (highPtr == band.end()) --highPtr;
      newStart = lowPtr - band.begin();
      newEnd = highPtr - band.begin() + 1;
    }

    selection.SetChannelRange(newStart, newEnd);
    return true;
  } else {
    return false;
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
                          _settings.channelsOut, _facets.size(),
                          _settings.prefixName + "-model");
  WSCFitsWriter writer(
      createWSCFitsWriter(groupTable.Front(), false, false, false));
  _residualImages.Initialize(writer.Writer(), _settings.polarizations.size(),
                             _settings.channelsOut, _facets.size(),
                             _settings.prefixName + "-residual");
  if (groupTable.Front().polarization == *_settings.polarizations.begin())
    _psfImages.Initialize(writer.Writer(), 1, groupTable.SquaredGroups().size(),
                          _facets.size(), _settings.prefixName + "-psf");

  // In the case of IDG we have to directly ask for all four polarizations. This
  // can't be parallelized in the current structure.
  bool parallelizeChannels =
      !_settings.useIDG || _settings.polarizations.size() == 1;
  // In case XY/YX polarizations are requested, we should not parallelize over
  // those since they need to be combined after imaging, and this currently
  // requires XY before YX.
  bool parallelizePolarizations =
      _settings.polarizations.count(aocommon::Polarization::XY) == 0 &&
      _settings.polarizations.count(aocommon::Polarization::YX) == 0;

  const std::string rootPrefix = _settings.prefixName;

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

  const size_t facetCount = groupTable.FacetCount();
  for (size_t facetIndex = 0; facetIndex < facetCount; ++facetIndex) {
    ImagingTable facetTable = groupTable.GetFacet(facetIndex);

    if (parallelizePolarizations) {
      for (ImagingTableEntry& entry : facetTable) {
        runFirstInversion(entry, primaryBeam);
      }
    } else {
      bool hasMore;
      size_t sqIndex = 0;
      do {
        hasMore = false;
        for (const ImagingTable::Group& sqGroup : facetTable.SquaredGroups()) {
          if (sqIndex < sqGroup.size()) {
            hasMore = true;
            runFirstInversion(*sqGroup[sqIndex], primaryBeam);
          }
        }
        ++sqIndex;
        _griddingTaskManager->Finish();
      } while (hasMore);
    }
  }
  if (parallelizePolarizations) {
    _griddingTaskManager->Finish();
  }
  stitchFacets(groupTable, _residualImages, _settings.isDirtySaved, false);
  _inversionWatch.Pause();

  _deconvolution.InitializeDeconvolutionAlgorithm(
      groupTable.GetFacet(0), *_settings.polarizations.begin(),
      minTheoreticalBeamSize(groupTable), _settings.threadCount);

  if (!_settings.makePSFOnly) {
    if (_settings.deconvolutionIterationCount > 0) {
      // Start major cleaning loop
      _majorIterationNr = 1;
      bool reachedMajorThreshold = false;
      do {
        _deconvolution.InitializeImages(_residualImages, _modelImages,
                                        _psfImages);
        _deconvolutionWatch.Start();
        _deconvolution.Perform(groupTable.GetFacet(0), reachedMajorThreshold,
                               _majorIterationNr);
        _deconvolutionWatch.Pause();

        if (_majorIterationNr == 1 && _settings.deconvolutionMGain != 1.0 &&
            _settings.isFirstResidualSaved)
          writeFirstResidualImages(groupTable);

        if (!reachedMajorThreshold) writeModelImages(groupTable);

        if (_settings.deconvolutionMGain != 1.0) {
          if (parallelizeChannels && parallelizePolarizations) {
            _predictingWatch.Start();
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
              }  // end of polarization loop
            }    // end of joined channels loop
            _griddingTaskManager->Finish();
            _inversionWatch.Pause();
          } else if (parallelizePolarizations) {
            for (const ImagingTable::Group& sqGroup :
                 groupTable.SquaredGroups()) {
              _predictingWatch.Start();
              for (const ImagingTable::EntryPtr& entry : sqGroup) {
                predict(*entry);
              }
              _griddingTaskManager->Finish();
              _predictingWatch.Pause();

              _inversionWatch.Start();
              for (const ImagingTable::EntryPtr& entry : sqGroup) {
                imageMain(*entry, false, false);
              }  // end of polarization loop
              _griddingTaskManager->Finish();
              _inversionWatch.Pause();
            }       // end of joined channels loop
          } else {  // only parallize channels
            _predictingWatch.Start();
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
      Logger::Info << _majorIterationNr
                   << " major iterations were performed.\n";
    }

    for (const ImagingTableEntry& joinedEntry : groupTable) {
      if (joinedEntry.facetIndex == 0) {
        saveRestoredImagesForGroup(joinedEntry, primaryBeam);
      }
    }

    if (_settings.saveSourceList) {
      _deconvolution.SaveSourceList(groupTable, _observationInfo.phaseCentreRA,
                                    _observationInfo.phaseCentreDec);
      if (_settings.applyPrimaryBeam || _settings.gridWithBeam ||
          !_settings.atermConfigFilename.empty()) {
        _deconvolution.SavePBSourceList(groupTable,
                                        _observationInfo.phaseCentreRA,
                                        _observationInfo.phaseCentreDec);
      }
    }
  }

  _deconvolution.FreeDeconvolutionAlgorithms();

  Logger::Info << "Inversion: " << _inversionWatch.ToString()
               << ", prediction: " << _predictingWatch.ToString()
               << ", deconvolution: " << _deconvolutionWatch.ToString() << '\n';

  _settings.prefixName = rootPrefix;
}

void WSClean::saveRestoredImagesForGroup(
    const ImagingTableEntry& tableEntry,
    std::unique_ptr<PrimaryBeam>& primaryBeam) const {
  // Restore model to residual and save image
  size_t currentChannelIndex = tableEntry.outputChannelIndex;

  aocommon::PolarizationEnum curPol = tableEntry.polarization;
  for (size_t imageIter = 0; imageIter != tableEntry.imageCount; ++imageIter) {
    bool isImaginary = (imageIter == 1);
    WSCFitsWriter writer(
        createWSCFitsWriter(tableEntry, isImaginary, false, true));
    ImageF restoredImage(_settings.trimmedImageWidth,
                         _settings.trimmedImageHeight);
    _residualImages.Load(restoredImage.data(), curPol, currentChannelIndex,
                         isImaginary);

    if (_settings.deconvolutionIterationCount != 0)
      writer.WriteImage("residual.fits", restoredImage.data());

    if (_settings.isUVImageSaved)
      saveUVImage(restoredImage, tableEntry, isImaginary, "uv");

    ImageF modelImage(_settings.trimmedImageWidth,
                      _settings.trimmedImageHeight);
    _modelImages.Load(modelImage.data(), curPol, currentChannelIndex,
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
    ModelRenderer::Restore(
        restoredImage.data(), modelImage.data(), _settings.trimmedImageWidth,
        _settings.trimmedImageHeight, beamMaj, beamMin, beamPA,
        _settings.pixelScaleX, _settings.pixelScaleY);
    Logger::Info << "DONE\n";
    modelImage.reset();

    Logger::Info << "Writing restored image... ";
    Logger::Info.Flush();
    writer.WriteImage("image.fits", restoredImage.data());
    Logger::Info << "DONE\n";
    restoredImage.reset();

    if (curPol == *_settings.polarizations.rbegin()) {
      ImageFilename imageName =
          ImageFilename(currentChannelIndex, tableEntry.outputIntervalIndex);
      if (_settings.gridWithBeam || !_settings.atermConfigFilename.empty()) {
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
      } else if (_settings.applyPrimaryBeam) {
        primaryBeam->CorrectImages(writer.Writer(), imageName, "image");
        if (_settings.savePsfPb)
          primaryBeam->CorrectImages(writer.Writer(), imageName, "psf");
        if (_settings.deconvolutionIterationCount != 0) {
          primaryBeam->CorrectImages(writer.Writer(), imageName, "residual");
          primaryBeam->CorrectImages(writer.Writer(), imageName, "model");
        }
      }
    }
  }
}

void WSClean::writeFirstResidualImages(const ImagingTable& groupTable) const {
  Logger::Info << "Writing first iteration image(s)...\n";
  ImageF ptr(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
  for (const ImagingTableEntry& entry : groupTable) {
    size_t ch = entry.outputChannelIndex;
    if (entry.polarization == aocommon::Polarization::YX) {
      _residualImages.Load(ptr.data(), aocommon::Polarization::XY, ch, true);
      WSCFitsWriter writer(createWSCFitsWriter(
          entry, aocommon::Polarization::XY, true, false, false));
      writer.WriteImage("first-residual.fits", ptr.data());
    } else {
      _residualImages.Load(ptr.data(), entry.polarization, ch, false);
      WSCFitsWriter writer(createWSCFitsWriter(entry, false, false, false));
      writer.WriteImage("first-residual.fits", ptr.data());
    }
  }
}

void WSClean::writeModelImages(const ImagingTable& groupTable) const {
  Logger::Info << "Writing model image...\n";
  ImageF ptr(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
  for (const ImagingTableEntry& entry : groupTable) {
    size_t ch = entry.outputChannelIndex;
    if (entry.polarization == aocommon::Polarization::YX) {
      _modelImages.Load(ptr.data(), aocommon::Polarization::XY, ch, true);
      WSCFitsWriter writer(createWSCFitsWriter(
          entry, aocommon::Polarization::XY, true, true, true));
      writer.WriteImage("model.fits", ptr.data());
    } else {
      _modelImages.Load(ptr.data(), entry.polarization, ch, false);
      WSCFitsWriter writer(createWSCFitsWriter(entry, false, true, true));
      writer.WriteImage("model.fits", ptr.data());
    }
  }
}

void WSClean::initializeModelImages(const ImagingTableEntry& entry) {
  _modelImages.SetFitsWriter(
      createWSCFitsWriter(entry, false, true, false).Writer());

  if (_settings.continuedRun) {
    readExistingModelImages(entry);
  } else {
    // Set model to zero: already done if this is YX of XY/YX imaging combi
    if (!(entry.polarization == aocommon::Polarization::YX &&
          _settings.polarizations.count(aocommon::Polarization::XY) != 0)) {
      ImageF modelImage(_settings.trimmedImageWidth,
                        _settings.trimmedImageHeight, 0.0f);
      _modelImages.Store(modelImage.data(), entry.polarization,
                         entry.outputChannelIndex, false);
      if (aocommon::Polarization::IsComplex(entry.polarization))
        _modelImages.Store(modelImage.data(), entry.polarization,
                           entry.outputChannelIndex, true);
    }
  }
}

void WSClean::readExistingModelImages(const ImagingTableEntry& entry) {
  // load image(s) from disk and store them in the model-image cache.
  for (size_t i = 0; i != entry.imageCount; ++i) {
    std::string prefix = ImageFilename::GetPrefix(
        _settings, entry.polarization, entry.outputChannelIndex,
        entry.outputIntervalIndex, i == 1);
    FitsReader reader(prefix + "-model.fits");
    Logger::Info << "Reading " << reader.Filename() << "...\n";
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
    // Here I require the pixel scale to be accurate enough so that the image is
    // at most 1/10th pixel larger/smaller.
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

    // TODO check phase centre

    if (resetGridder)
      _griddingTaskManager = GriddingTaskManager::Make(_settings);

    if (!_imageWeightCache) {
      // The construction of the weight cache is delayed in prediction mode,
      // because only now the image size and scale is known.
      _imageWeightCache = createWeightCache();
      if (_settings.mfWeighting) initializeMFSImageWeights();
    }

    FitsWriter writer(reader);
    _modelImages.SetFitsWriter(writer);

    ImageF buffer(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
    reader.Read(buffer.data());
    for (size_t j = 0;
         j != _settings.trimmedImageWidth * _settings.trimmedImageHeight; ++j) {
      if (!std::isfinite(buffer[j]))
        throw std::runtime_error(
            "The input image contains non-finite values -- can't predict from "
            "an image with non-finite values");
    }
    _modelImages.Store(buffer.data(), entry.polarization,
                       entry.outputChannelIndex, i == 1);
  }
}

void WSClean::predictGroup(const ImagingTable::Group& imagingGroup) {
  _modelImages.Initialize(
      createWSCFitsWriter(*imagingGroup.front(), false, true, false).Writer(),
      _settings.polarizations.size(), 1, _facets.size(),
      _settings.prefixName + "-model");

  const std::string rootPrefix = _settings.prefixName;

  _predictingWatch.Start();
  for (const ImagingTable::EntryPtr& entry : imagingGroup) {
    const bool calculatePixelPositions = _settings.trimmedImageWidth == 0;
    readExistingModelImages(*entry);

    if (calculatePixelPositions) {
      for (schaapcommon::facets::Facet& facet : _facets) {
        const size_t alignment = 4;
        facet.CalculatePixels(
            _observationInfo.phaseCentreRA, _observationInfo.phaseCentreDec,
            _settings.pixelScaleX, _settings.pixelScaleY,
            _settings.trimmedImageWidth, _settings.trimmedImageHeight,
            _observationInfo.shiftL, _observationInfo.shiftM,
            _settings.imagePadding, alignment, _settings.useIDG);
      }
    }

    predict(*entry);
  }  // end of polarization loop
  _predictingWatch.Pause();

  _griddingTaskManager->Finish();

  Logger::Info << "Inversion: " << _inversionWatch.ToString()
               << ", prediction: " << _predictingWatch.ToString()
               << ", cleaning: " << _deconvolutionWatch.ToString() << '\n';

  _settings.prefixName = rootPrefix;
}

void WSClean::initializeMSList(
    const ImagingTableEntry& entry,
    std::vector<std::unique_ptr<MSDataDescription>>& msList) {
  aocommon::PolarizationEnum pol = _settings.useIDG
                                       ? aocommon::Polarization::Instrumental
                                       : entry.polarization;

  msList.clear();
  for (size_t i = 0; i != _settings.filenames.size(); ++i) {
    for (size_t d = 0; d != _msBands[i].DataDescCount(); ++d) {
      MSSelection selection(_globalSelection);
      if (selectChannels(selection, i, d, entry)) {
        std::unique_ptr<MSDataDescription> dataDescription;
        if (_settings.doReorder)
          dataDescription = MSDataDescription::ForPartitioned(
              _partitionedMSHandles[i], selection,
              entry.msData[i].bands[d].partIndex, pol, d);
        else
          dataDescription = MSDataDescription::ForContiguous(
              _settings.filenames[i], _settings.dataColumnName, selection, pol,
              d);
        msList.emplace_back(std::move(dataDescription));
      }
    }
  }
}

void WSClean::runFirstInversion(
    ImagingTableEntry& entry, std::unique_ptr<class PrimaryBeam>& primaryBeam) {
  const bool isLastPol =
      entry.polarization == *_settings.polarizations.rbegin();
  const bool doMakePSF = _settings.deconvolutionIterationCount > 0 ||
                         _settings.makePSF || _settings.makePSFOnly;

  if (isLastPol) {
    ImageFilename imageName =
        ImageFilename(entry.outputChannelIndex, entry.outputIntervalIndex);
    if (_settings.applyPrimaryBeam) {
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
      primaryBeam->MakeBeamImages(imageName, entry, std::move(weights));
    }
  }

  if (!_settings.makePSFOnly) {
    if (_settings.reuseDirty)
      loadExistingDirty(entry, !doMakePSF);
    else
      imageMain(entry, true, !doMakePSF);

    _isFirstInversion = false;
  }
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
      // Store the full interval in the selection, so that it doesn't need to be
      // determined again.
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

void WSClean::saveUVImage(const ImageF& image, const ImagingTableEntry& entry,
                          bool isImaginary, const std::string& prefix) const {
  ImageF realUV(_settings.trimmedImageWidth, _settings.trimmedImageHeight),
      imagUV(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
  FFTResampler fft(_settings.trimmedImageWidth, _settings.trimmedImageHeight,
                   _settings.trimmedImageWidth, _settings.trimmedImageHeight, 1,
                   true);
  fft.SingleFT(image.data(), realUV.data(), imagUV.data());
  // Factors of 2 involved: because of SingleFT()
  // (also one from the fact that normF excludes a factor of two?)
  realUV *=
      _infoPerChannel[entry.outputChannelIndex].normalizationFactor /
      sqrt(0.5 * _settings.trimmedImageWidth * _settings.trimmedImageHeight);
  imagUV *=
      _infoPerChannel[entry.outputChannelIndex].normalizationFactor /
      sqrt(0.5 * _settings.trimmedImageWidth * _settings.trimmedImageHeight);
  WSCFitsWriter writer(createWSCFitsWriter(entry, isImaginary, false, false));
  writer.WriteUV(prefix + "-real.fits", realUV.data());
  writer.WriteUV(prefix + "-imag.fits", imagUV.data());
}

void WSClean::stitchFacets(const ImagingTable& table,
                           CachedImageSet& imageCache, bool writeDirty,
                           bool isPSF) {
  if (!_facets.empty()) {
    Logger::Info << "Stitching facets onto full image...\n";
    // Allocate full image
    ImageF fullImage(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
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
                            isPSF, fullImage, facetImage);
        }
      }
    }
  }
}

void WSClean::stitchSingleGroup(const ImagingTable& facetGroup,
                                size_t imageIndex, CachedImageSet& imageCache,
                                bool writeDirty, bool isPSF, ImageF& fullImage,
                                schaapcommon::facets::FacetImage& facetImage) {
  const bool isImaginary = (imageIndex == 1);
  fullImage = 0.0f;
  for (const ImagingTableEntry& facetEntry : facetGroup) {
    facetImage.SetFacet(*facetEntry.facet, true);
    imageCache.LoadFacet(facetImage.Data(0), facetEntry.polarization,
                         facetEntry.outputChannelIndex, facetEntry.facetIndex,
                         facetEntry.facet, isImaginary);
    // TODO with our current stitching implementation, facets should always be
    // directly copied to the full image, not added. The facets should not
    // overlap though.
    facetImage.AddToImage({fullImage.data()});
  }
  if (writeDirty) {
    initializeModelImages(facetGroup.Front());
    _residualImages.SetFitsWriter(
        createWSCFitsWriter(facetGroup.Front(), false, false, true).Writer());
    WSCFitsWriter writer(
        createWSCFitsWriter(facetGroup.Front(), isImaginary, false, true));
    Logger::Info << "Writing dirty image...\n";
    writer.WriteImage("dirty.fits", fullImage.data());
  }

  if (isPSF) {
    const ImagingTableEntry& entry = facetGroup.Front();
    processFullPSF(fullImage, entry);
  }

  const size_t channelIndex = facetGroup.Front().outputChannelIndex;
  const aocommon::PolarizationEnum polarization =
      facetGroup.Front().polarization;
  imageCache.Store(fullImage.data(), polarization, channelIndex, isImaginary);
}

void WSClean::makeImagingTable(size_t outputIntervalIndex) {
  std::set<aocommon::ChannelInfo> channelSet;
  _msBands.assign(_settings.filenames.size(), MultiBandData());
  for (size_t i = 0; i != _settings.filenames.size(); ++i) {
    casacore::MeasurementSet ms(_settings.filenames[i]);
    _msBands[i] = MultiBandData(ms.spectralWindow(), ms.dataDescription());
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
      size_t outChannelStart = _settings.channelsOut * i / (nSplits + 1),
             outChannelEnd = _settings.channelsOut * (i + 1) / (nSplits + 1);
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
    double chFlux = NonLinearPowerLawFitter::Evaluate(
        chCentralFrequency, _settings.spectralCorrection,
        _settings.spectralCorrectionFrequency);
    double midFlux = NonLinearPowerLawFitter::Evaluate(
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

void WSClean::addFacetsToImagingTable(ImagingTableEntry& templateEntry) {
  if (_facets.empty()) {
    std::unique_ptr<ImagingTableEntry> entry(
        new ImagingTableEntry(templateEntry));
    entry->facetIndex = 0;
    entry->facet = nullptr;
    _imagingTable.AddEntry(std::move(entry));
  } else {
    for (size_t f = 0; f != _facets.size(); ++f) {
      std::unique_ptr<ImagingTableEntry> entry(
          new ImagingTableEntry(templateEntry));
      entry->facetIndex = f;
      entry->facet = &_facets[f];

      // Calculate phase center delta for entry
      entry->centreShiftX = _facets[f].GetUntrimmedBoundingBox().Centre().x -
                            _settings.trimmedImageWidth / 2;
      entry->centreShiftY = _facets[f].GetUntrimmedBoundingBox().Centre().y -
                            _settings.trimmedImageHeight / 2;
      _imagingTable.AddEntry(std::move(entry));
    }
  }
  ++templateEntry.facetGroupIndex;
}

void WSClean::addPolarizationsToImagingTable(ImagingTableEntry& templateEntry) {
  for (aocommon::PolarizationEnum p : _settings.polarizations) {
    templateEntry.polarization = p;
    if (p == aocommon::Polarization::XY)
      templateEntry.imageCount = 2;
    else if (p == aocommon::Polarization::YX)
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

WSCFitsWriter WSClean::createWSCFitsWriter(
    const ImagingTableEntry& entry, aocommon::PolarizationEnum polarization,
    bool isImaginary, bool isModel, bool isFullImage) const {
  ObservationInfo observationInfo = _observationInfo;
  if (!isFullImage) {
    applyFacetPhaseShift(entry, observationInfo);
  }

  return WSCFitsWriter(entry, polarization, isImaginary, _settings,
                       _deconvolution, observationInfo, _majorIterationNr,
                       _commandLine, _infoPerChannel[entry.outputChannelIndex],
                       isModel, _lastStartTime);
}
