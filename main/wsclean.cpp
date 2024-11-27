#include "wsclean.h"

#include "../math/imageoperations.h"

#include "../gridding/directmsgridder.h"

#include "../io/componentlistwriter.h"
#include "../io/facetreader.h"
#include "../io/imagefilename.h"
#include "../io/imageweightcache.h"
#include "../io/wscfitswriter.h"

#include "../scheduling/griddingtaskmanager.h"

#include "../system/application.h"

#include "../structures/facetutil.h"
#include "../structures/imageweights.h"
#include "../structures/msselection.h"
#include "../structures/primarybeam.h"

#include <radler/radler.h>

#include "../idg/averagebeam.h"
#include "../idg/idgmsgridder.h"

#include "../math/renderer.h"
#include "../math/tophatconvolution.h"

#include "../model/model.h"

#include "../msproviders/contiguousms.h"
#include "../msproviders/msdatadescription.h"

#include "progressbar.h"

#include <aocommon/fits/fitswriter.h>
#include <aocommon/image.h>
#include <aocommon/logger.h>
#include <aocommon/uvector.h>
#include <aocommon/units/angle.h>

#include <schaapcommon/facets/facetimage.h>
#include <schaapcommon/math/resampler.h>
#include <schaapcommon/math/restoreimage.h>
#include <schaapcommon/fitters/nlplfitter.h>

#include <algorithm>
#include <iostream>
#include <memory>

using aocommon::Image;
using aocommon::Logger;
using aocommon::Polarization;
using aocommon::PolarizationEnum;
using aocommon::units::Angle;

namespace wsclean {

WSClean::WSClean()
    : _globalSelection(),
      _commandLine(),
      _inversionWatch(false),
      _predictingWatch(false),
      _deconvolutionWatch(false),
      _isFirstInversionTask(true),
      _majorIterationNr(0),
      _psfImages(),
      _modelImages(),
      _residualImages(),
      _deconvolution(),
      _ddPsfCount(0),
      _lastStartTime(0.0) {}

WSClean::~WSClean() = default;

void WSClean::loadExistingImage(ImagingTableEntry& entry, bool isPSF) {
  if (isPSF)
    Logger::Info << "Loading existing PSF from disk...\n";
  else
    Logger::Info << "Loading existing dirty image from disk...\n";

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
  Image image(reader.ImageWidth(), reader.ImageHeight());
  reader.Read(image.Data());

  const size_t channel_index = entry.outputChannelIndex;
  OutputChannelInfo& channel_info = _infoPerChannel[channel_index];

  entry.imageWeight = reader.ReadDoubleKey("WSCIMGWG");
  channel_info.weight = entry.imageWeight;
  reader.ReadDoubleKeyIfExists("WSCVWSUM", channel_info.visibilityWeightSum);
  double n_visibilities = channel_info.visibilityCount;
  reader.ReadDoubleKeyIfExists("WSCNVIS", n_visibilities);
  channel_info.visibilityCount = static_cast<size_t>(n_visibilities);
  reader.ReadDoubleKeyIfExists("WSCENVIS",
                               channel_info.effectiveVisibilityCount);
  reader.ReadDoubleKeyIfExists("WSCNORMF", channel_info.normalizationFactor);
  entry.normalizationFactor = channel_info.normalizationFactor;

  if (isPSF) {
    processFullPSF(image, entry);

    _psfImages.SetWSCFitsWriter(createWSCFitsWriter(entry, false, false));
    _psfImages.Store(image.Data(), *_settings.polarizations.begin(),
                     channel_index, false);
  } else {
    // maxFacetGroupIndex is always 1
    const size_t maxFacetGroupIndex = 1;
    initializeModelImages(entry, *_settings.polarizations.begin(),
                          maxFacetGroupIndex);
    _residualImages.SetWSCFitsWriter(createWSCFitsWriter(
        entry, *_settings.polarizations.begin(), false, false));
    for (size_t imageIndex = 0; imageIndex != entry.imageCount; ++imageIndex) {
      const bool isImaginary = (imageIndex == 1);
      WSCFitsWriter writer(createWSCFitsWriter(
          entry, *_settings.polarizations.begin(), isImaginary, false));
      _residualImages.Store(image, *_settings.polarizations.begin(),
                            channel_index, isImaginary);
      if (_settings.isDirtySaved) {
        Logger::Info << "Writing dirty image...\n";
        writer.WriteImage("dirty.fits", image);
      }
    }
  }
}

void WSClean::storeAverageBeam(const ImagingTableEntry& entry,
                               std::unique_ptr<AverageBeam>& averageBeam) {
  if (averageBeam) {
    _scalarBeamImages.SetWSCFitsWriter(
        createWSCFitsWriter(entry, false, false));
    _matrixBeamImages.SetWSCFitsWriter(
        createWSCFitsWriter(entry, false, false));
    averageBeam->Store(_scalarBeamImages, _matrixBeamImages,
                       entry.outputChannelIndex);
  }
}

void WSClean::ImagePsf(ImagingTable::Group&& facet_group) {
  std::vector<GriddingTask> tasks = _griddingTaskFactory->CreatePsfTasks(
      facet_group, *_imageWeightCache, _settings.compound_tasks,
      _isFirstInversionTask);

  if (!_settings.compound_tasks) {
    for (size_t i = 0; i < facet_group.size(); ++i) {
      Logger::Info.Flush();
      Logger::Info << " == Constructing PSF ==\n";

      _griddingTaskManager->Run(
          std::move(tasks[i]),
          [this, entry = facet_group[i]](GriddingResult& result) {
            ImagePsfCallback({std::move(entry)}, result);
          });
    }
  } else {
    assert(tasks.size() == 1 &&
           tasks.front().facets.size() == facet_group.size());

    Logger::Info.Flush();
    Logger::Info << " == Constructing PSF group ==\n";

    _griddingTaskManager->Run(
        std::move(tasks.front()),
        [this, facet_group = std::move(facet_group)](GriddingResult& result) {
          ImagePsfCallback(std::move(facet_group), result);
        });
  }
}

void WSClean::ImagePsfCallback(ImagingTable::Group facet_group,
                               GriddingResult& result) {
  assert(facet_group.size() == result.facets.size());

  const size_t channel_index = facet_group.front()->outputChannelIndex;
  OutputChannelInfo& channel_info = _infoPerChannel[channel_index];
  _lastStartTime = result.startTime;

  channel_info.beamSizeEstimate = result.beamSize;
  channel_info.visibilityCount = result.griddedVisibilityCount;
  channel_info.visibilityWeightSum = result.visibilityWeightSum;

  for (size_t i = 0; i < facet_group.size(); ++i) {
    ImagingTableEntry& entry = *facet_group[i];
    GriddingResult::FacetData& facet_result = result.facets[i];

    entry.imageWeight = facet_result.imageWeight;
    entry.normalizationFactor = facet_result.normalizationFactor;
    channel_info.weight = entry.imageWeight;
    channel_info.normalizationFactor = entry.normalizationFactor;
    channel_info.wGridSize = facet_result.actualWGridSize;
    channel_info.effectiveVisibilityCount =
        facet_result.effectiveGriddedVisibilityCount;

    if (entry.isDdPsf) {
      channel_info.averageDdPsfCorrection[entry.facetIndex] =
          facet_result.averageCorrection;
    } else {
      channel_info.averageFacetCorrection[entry.facetIndex] =
          facet_result.averageCorrection;
      channel_info.averageBeamFacetCorrection[entry.facetIndex] =
          facet_result.averageBeamCorrection;
    }

    _griddingTaskFactory->SetMetaDataCacheEntry(entry,
                                                std::move(facet_result.cache));

    if (entry.isDdPsf || _facetCount == 0)
      processFullPSF(facet_result.images[0], entry);

    _psfImages.SetWSCFitsWriter(createWSCFitsWriter(entry, false, false));
    _psfImages.StoreFacet(facet_result.images[0],
                          *_settings.polarizations.begin(), channel_index,
                          entry.facetIndex, entry.facet, false);
  }
}

void WSClean::processFullPSF(Image& image, const ImagingTableEntry& entry) {
  Settings settings(_settings);
  settings.trimmedImageWidth = image.Width();
  settings.trimmedImageHeight = image.Height();
  const size_t centralIndex =
      settings.trimmedImageWidth / 2 +
      (settings.trimmedImageHeight / 2) * settings.trimmedImageWidth;
  const size_t channelIndex = entry.outputChannelIndex;

  // When imaging with dd psfs, facet corrections are applied, so correct for
  // this. di psfs do not receive facet corrections and do not require this
  // multiplication.
  if (entry.isDdPsf && _settings.UseFacetCorrections()) {
    // This is by itself an ineffective correction, because PSFs are normalized
    // to unity in the centre, so this scaling is lost. It is still useful
    // because by correcting this we get a PSF normalization factor that is more
    // sensible (see ~10 lines below this), and this may be displayed to the
    // user. It's also used to scale the normal images in the same way.
    //
    // Nevertheless, for the PSF we do not correct polarization leakage, because
    // we only make a Stokes I beam, hence we cannot do a full Mueller matrix
    // correction as we do for facets.
    //
    // One could imagine scenarios where the I and the Q,U,V PSFs are different.
    // For example when most of the time either the X dipole or the Y dipole is
    // active, and only a fraction of the time both dipoles are active. Then
    // Stokes I still has good uv-coverage, while Q,U,V have poor uv-coverage.
    // In practice that won't happen. The final deconvolution result is not
    // sensitive to small errors in the PSF anyway.
    const double factor = _infoPerChannel[channelIndex]
                              .averageDdPsfCorrection[entry.facetIndex]
                              .GetStokesIValue();
    image *= 1.0 / std::sqrt(factor);
  }

  double normFactor;
  if (image[centralIndex] != 0.0)
    normFactor = 1.0 / image[centralIndex];
  else
    normFactor = 0.0;

  image *= normFactor * entry.siCorrection;
  Logger::Debug << "Normalized PSF by factor of " << normFactor << ".\n";

  image.RemoveNans();
  double minPixelScale = std::min(settings.pixelScaleX, settings.pixelScaleY);
  double initialFitSize =
      std::max(_infoPerChannel[channelIndex].beamSizeEstimate, minPixelScale);
  double bMaj, bMin, bPA, bTheoretical;
  ImageOperations::DetermineBeamSize(settings, bMaj, bMin, bPA, bTheoretical,
                                     image, initialFitSize);
  // Create a temporary copy of the output channel info
  // and put the fitting and normalization result in the copy
  OutputChannelInfo channel_info(_infoPerChannel[channelIndex]);
  channel_info.theoreticBeamSize = bTheoretical;
  channel_info.beamMaj = bMaj;
  channel_info.beamMin = bMin;
  channel_info.beamPA = bPA;
  channel_info.psfNormalizationFactor = normFactor;

  // If this entry is the main psf, or if it is the dd-psf at the centre
  // then use the fitting result for the main image
  if (!entry.isDdPsf ||
      entry.facet->Contains(schaapcommon::facets::PixelPosition(
          _settings.trimmedImageWidth / 2, _settings.trimmedImageHeight / 2))) {
    _infoPerChannel[channelIndex] = channel_info;
    _infoPerChannel[channelIndex].centralPsfIndex = entry.facetIndex;
    _infoForMFS.centralPsfIndex = entry.facetIndex;
  }

  Logger::Info << "Writing psf image... ";
  if (settings.isUVImageSaved) {
    saveUVImage(image, entry, channel_info, false, "uvpsf");
  }

  Logger::Info.Flush();
  const std::string name(
      (entry.isDdPsf ? ImageFilename::GetPSFPrefix(settings, channelIndex,
                                                   entry.outputIntervalIndex,
                                                   entry.facetIndex)
                     : ImageFilename::GetPSFPrefix(settings, channelIndex,
                                                   entry.outputIntervalIndex)) +
      "-psf.fits");
  WSCFitsWriter fits_writer =
      createWSCFitsWriter(entry, channel_info, false, false);
  if (entry.isDdPsf) {
    fits_writer.WriteFullNameImage(name, image, *entry.facet);
  } else {
    fits_writer.WriteFullNameImage(name, image);
  }
  Logger::Info << "DONE\n";
}

void WSClean::ImageMain(ImagingTable::Group& facet_group, bool isFirstInversion,
                        bool updateBeamInfo) {
  std::vector<std::unique_ptr<AverageBeam>> average_beams;
  if (!isFirstInversion && griddingUsesATerms()) {
    average_beams.reserve(facet_group.size());
    for (const std::shared_ptr<ImagingTableEntry>& entry : facet_group) {
      average_beams.push_back(AverageBeam::Load(
          _scalarBeamImages, _matrixBeamImages, entry->outputChannelIndex));
    }
  }

  std::vector<GriddingTask> tasks = _griddingTaskFactory->CreateInvertTasks(
      facet_group, *_imageWeightCache, _settings.compound_tasks,
      _isFirstInversionTask, isFirstInversion, std::move(average_beams));

  if (!_settings.compound_tasks) {
    for (size_t i = 0; i < facet_group.size(); ++i) {
      Logger::Info.Flush();
      Logger::Info << " == Constructing image ==\n";

      _griddingTaskManager->Run(
          std::move(tasks[i]), [this, entry = facet_group[i], updateBeamInfo,
                                isFirstInversion](GriddingResult& result) {
            ImageMainCallback({entry}, result, updateBeamInfo,
                              isFirstInversion);
          });
    }
  } else {
    assert(tasks.size() == 1 &&
           tasks.front().facets.size() == facet_group.size());

    Logger::Info.Flush();
    Logger::Info << " == Constructing image group ==\n";

    // Copy 'facet_group' by value into the callback function, and then
    // move that value when invoking ImageMainCallback.
    _griddingTaskManager->Run(
        std::move(tasks.front()), [this, facet_group, updateBeamInfo,
                                   isFirstInversion](GriddingResult& result) {
          ImageMainCallback(std::move(facet_group), result, updateBeamInfo,
                            isFirstInversion);
        });
  }

  _isFirstInversionTask = false;
}

void WSClean::ImageMainCallback(ImagingTable::Group facet_group,
                                GriddingResult& result, bool updateBeamInfo,
                                bool isFirstInversion) {
  assert(facet_group.size() == result.facets.size());

  const size_t channel_index = facet_group.front()->outputChannelIndex;
  OutputChannelInfo& channel_info = _infoPerChannel[channel_index];
  _lastStartTime = result.startTime;

  // If no PSF is made, also set the beam size. If the PSF was made, these
  // would already be set after imaging the PSF.
  if (updateBeamInfo) {
    if (_settings.theoreticBeam) {
      channel_info.beamMaj =
          std::max(result.beamSize, _settings.gaussianTaperBeamSize);
      channel_info.beamMin =
          std::max(result.beamSize, _settings.gaussianTaperBeamSize);
      channel_info.beamPA = 0.0;
    } else if (_settings.manualBeamMajorSize != 0.0) {
      channel_info.beamMaj = _settings.manualBeamMajorSize;
      channel_info.beamMin = _settings.manualBeamMinorSize;
      channel_info.beamPA = _settings.manualBeamPA;
    } else {
      channel_info.beamMaj = std::numeric_limits<double>::quiet_NaN();
      channel_info.beamMin = std::numeric_limits<double>::quiet_NaN();
      channel_info.beamPA = std::numeric_limits<double>::quiet_NaN();
    }
  }

  for (size_t i = 0; i < facet_group.size(); ++i) {
    ImagingTableEntry& entry = *facet_group[i];
    GriddingResult::FacetData& facet_result = result.facets[i];

    _griddingTaskFactory->SetMetaDataCacheEntry(entry,
                                                std::move(facet_result.cache));
    entry.imageWeight = facet_result.imageWeight;
    entry.normalizationFactor = facet_result.normalizationFactor;
    channel_info.weight = facet_result.imageWeight;
    channel_info.normalizationFactor = facet_result.normalizationFactor;

    if (!facet_result.averageCorrection.IsZero()) {
      channel_info.averageFacetCorrection[entry.facetIndex] =
          facet_result.averageCorrection;
      channel_info.averageBeamFacetCorrection[entry.facetIndex] =
          facet_result.averageBeamCorrection;
    }

    using PolImagesPair = std::pair<const PolarizationEnum, std::vector<Image>>;
    std::vector<PolImagesPair> imageList;
    if (_settings.gridderType == GridderType::IDG &&
        _settings.polarizations.size() != 1) {
      assert(facet_result.images.size() == _settings.polarizations.size());
      imageList.reserve(facet_result.images.size());
      auto polIter = _settings.polarizations.begin();
      for (size_t polIndex = 0; polIndex != facet_result.images.size();
           ++polIndex, ++polIter)
        imageList.emplace_back(
            *polIter,
            std::vector<Image>{std::move(facet_result.images[polIndex])});
    } else {
      imageList.emplace_back(entry.polarization,
                             std::move(facet_result.images));
    }

    for (PolImagesPair& polImagePair : imageList) {
      std::vector<Image>& images = polImagePair.second;
      const PolarizationEnum polarization = polImagePair.first;
      for (size_t i = 0; i != images.size(); ++i) {
        // IDG performs normalization on the dirty images, so only normalize if
        // not using IDG
        const double psfFactor = _settings.gridderType == GridderType::IDG
                                     ? 1.0
                                     : channel_info.psfNormalizationFactor;
        images[i] *= psfFactor * entry.siCorrection;
        const bool isImaginary = i == 1;
        storeAndCombineXYandYX(_residualImages, channel_index, entry,
                               polarization, isImaginary, images[i]);
      }

      // If facets are used, stitchFacets() performs these actions.
      if (isFirstInversion && 0 == _facetCount) {
        // maxFacetGroupIndex is always 1
        const size_t maxFacetGroupIndex = 1;
        initializeModelImages(entry, polarization, maxFacetGroupIndex);
        _residualImages.SetWSCFitsWriter(
            createWSCFitsWriter(entry, polarization, false, false));
        // If facets are used, stitchFacets() saves the dirty image.
        if (_settings.isDirtySaved) {
          for (size_t imageIndex = 0; imageIndex != entry.imageCount;
               ++imageIndex) {
            const bool isImaginary = (imageIndex == 1);
            WSCFitsWriter writer(
                createWSCFitsWriter(entry, polarization, isImaginary, false));
            Image dirtyImage(_settings.trimmedImageWidth,
                             _settings.trimmedImageHeight);
            _residualImages.Load(dirtyImage.Data(), polarization,
                                 entry.outputChannelIndex, isImaginary);
            Logger::Info << "Writing dirty image...\n";
            writer.WriteImage("dirty.fits", dirtyImage);
          }
        }
      }
    }

    storeAverageBeam(entry, facet_result.averageBeam);

    if (isFirstInversion && griddingUsesATerms()) {
      Logger::Info << "Writing IDG beam image...\n";
      ImageFilename imageName(entry.outputChannelIndex,
                              entry.outputIntervalIndex);
      if (!facet_result.averageBeam || facet_result.averageBeam->Empty()) {
        throw std::runtime_error(
            "Trying to write the IDG beam while the beam has not been computed "
            "yet.");
      }
      IdgMsGridder::SaveBeamImage(entry, imageName, _settings,
                                  _observationInfo.phaseCentreRA,
                                  _observationInfo.phaseCentreDec, _l_shift,
                                  _m_shift, *facet_result.averageBeam);
    }
  }
}

void WSClean::storeAndCombineXYandYX(CachedImageSet& dest,
                                     size_t joinedChannelIndex,
                                     const ImagingTableEntry& entry,
                                     PolarizationEnum polarization,
                                     bool isImaginary, const Image& image) {
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
    dest.StoreFacet(xyImage, Polarization::XY, joinedChannelIndex,
                    entry.facetIndex, entry.facet, isImaginary);
  } else {
    dest.StoreFacet(image, polarization, joinedChannelIndex, entry.facetIndex,
                    entry.facet, isImaginary);
  }
}

void WSClean::Predict(const ImagingTable::Group& facet_group) {
  const bool is_full_stokes = _settings.gridderType == GridderType::IDG &&
                              _settings.polarizations.size() != 1;
  std::vector<std::vector<Image>> model_images;
  std::vector<std::unique_ptr<AverageBeam>> average_beams;
  model_images.reserve(facet_group.size());
  average_beams.reserve(facet_group.size());

  for (const std::shared_ptr<ImagingTableEntry>& entry : facet_group) {
    Logger::Info.Flush();
    Logger::Info << " == Converting model image to visibilities ==\n";

    size_t width;
    size_t height;
    if (entry->facet) {
      width = entry->facet->GetTrimmedBoundingBox().Width();
      height = entry->facet->GetTrimmedBoundingBox().Height();
    } else {
      width = _settings.trimmedImageWidth;
      height = _settings.trimmedImageHeight;
    }

    std::vector<PolarizationEnum> polarizations;
    if (is_full_stokes)
      polarizations.assign(_settings.polarizations.begin(),
                           _settings.polarizations.end());
    else
      polarizations = {entry->polarization};

    std::vector<Image> entry_model_images;
    entry_model_images.reserve(polarizations.size());
    for (PolarizationEnum& polarization : polarizations) {
      entry_model_images.emplace_back(width, height);
      bool is_yx = polarization == Polarization::YX;
      const PolarizationEnum loadPol = is_yx ? Polarization::XY : polarization;
      _modelImages.LoadFacet(entry_model_images.back().Data(), loadPol,
                             entry->outputChannelIndex, entry->facetIndex,
                             entry->facet, false);
      if (Polarization::IsComplex(polarization)) {  // XY or YX
        entry_model_images.emplace_back(width, height);
        // YX is never stored: it is always combined with XY and stored as XY
        _modelImages.LoadFacet(entry_model_images.back().Data(),
                               Polarization::XY, entry->outputChannelIndex,
                               entry->facetIndex, entry->facet, true);
        if (is_yx) {
          for (float& v : entry_model_images.back()) v = -v;
        }
      }
    }
    model_images.push_back(std::move(entry_model_images));

    average_beams.push_back(AverageBeam::Load(
        _scalarBeamImages, _matrixBeamImages, entry->outputChannelIndex));
  }

  std::vector<GriddingTask> tasks = _griddingTaskFactory->CreatePredictTasks(
      facet_group, *_imageWeightCache, _settings.compound_tasks,
      std::move(model_images), std::move(average_beams));

  if (!_settings.compound_tasks) {
    assert(tasks.size() == facet_group.size());
    for (size_t i = 0; i < facet_group.size(); ++i) {
      _griddingTaskManager->Run(
          std::move(tasks[i]),
          [this, entry = facet_group[i]](GriddingResult& result) {
            GriddingResult::FacetData& facet_result = result.facets.front();
            _griddingTaskFactory->SetMetaDataCacheEntry(
                *entry, std::move(facet_result.cache));
          });
    }
  } else {
    assert(tasks.size() == 1 &&
           tasks.front().facets.size() == facet_group.size());

    _griddingTaskManager->Run(
        std::move(tasks.front()), [this, facet_group](GriddingResult& result) {
          assert(facet_group.size() == result.facets.size());

          for (size_t i = 0; i < facet_group.size(); ++i) {
            _griddingTaskFactory->SetMetaDataCacheEntry(
                *facet_group[i], std::move(result.facets[i].cache));
          }
        });
  }
}

ObservationInfo WSClean::getObservationInfo() const {
  casacore::MeasurementSet ms(_settings.filenames[0]);
  ObservationInfo observationInfo =
      ReadObservationInfo(ms, _settings.fieldIds[0]);
  return observationInfo;
}

std::pair<double, double> WSClean::getLMShift() const {
  double l_shift = 0.0;
  double m_shift = 0.0;
  if (_settings.hasShift) {
    aocommon::ImageCoordinates::RaDecToLM(
        _settings.shiftRA, _settings.shiftDec, _observationInfo.phaseCentreRA,
        _observationInfo.phaseCentreDec, l_shift, m_shift);
  }
  return std::make_pair(l_shift, m_shift);
}

void WSClean::RunClean() {
  _observationInfo = getObservationInfo();
  std::tie(_l_shift, _m_shift) = getLMShift();

  std::vector<std::shared_ptr<schaapcommon::facets::Facet>> facets =
      FacetReader::ReadFacets(
          _settings.facetRegionFilename, _settings.trimmedImageWidth,
          _settings.trimmedImageHeight, _settings.pixelScaleX,
          _settings.pixelScaleY, _observationInfo.phaseCentreRA,
          _observationInfo.phaseCentreDec, _l_shift, _m_shift,
          _settings.imagePadding, _settings.gridderType == GridderType::IDG,
          _settings.GetFeatherSize());
  _facetCount = facets.size();

  std::vector<std::shared_ptr<schaapcommon::facets::Facet>> dd_psfs;
  if (_settings.ddPsfGridWidth > 1 || _settings.ddPsfGridHeight > 1) {
    const schaapcommon::facets::Facet::InitializationData facet_data =
        CreateFacetInitializationData(
            _settings.trimmedImageWidth, _settings.trimmedImageHeight,
            _settings.pixelScaleX, _settings.pixelScaleY,
            _observationInfo.phaseCentreRA, _observationInfo.phaseCentreDec,
            _l_shift, _m_shift, _settings.imagePadding,
            _settings.gridderType == GridderType::IDG, 0);
    dd_psfs = CreateFacetGrid(facet_data, _settings.ddPsfGridWidth,
                              _settings.ddPsfGridHeight);
  }
  _ddPsfCount = dd_psfs.size();

  schaapcommon::facets::PixelPosition centerPixel(
      _settings.trimmedImageWidth / 2, _settings.trimmedImageHeight / 2);
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

  MSSelection fullSelection = _settings.GetMSSelection();

  for (size_t intervalIndex = 0; intervalIndex != _settings.intervalsOut;
       ++intervalIndex) {
    makeImagingTable(intervalIndex);
    if (!facets.empty()) updateFacetsInImagingTable(facets, false);
    if (!dd_psfs.empty()) updateFacetsInImagingTable(dd_psfs, true);

    _globalSelection = selectInterval(fullSelection, intervalIndex);

    _msHelper =
        std::make_unique<MsHelper>(_settings, _globalSelection, _msBands);

    // Read reordered files if the option reuse-reordered is
    // specified.
    if (_settings.reuseReorder)
      _msHelper->ReuseReorderedFiles(_imagingTable);
    else if (_settings.doReorder)
      _msHelper->PerformReordering(_imagingTable, false);

    _infoPerChannel.assign(_settings.channelsOut,
                           OutputChannelInfo(std::max<size_t>(1, _facetCount),
                                             std::max<size_t>(1, _ddPsfCount)));

    _imageWeightCache = createWeightCache();

    _image_weight_initializer = std::make_unique<ImageWeightInitializer>(
        _settings, _globalSelection, _msBands,
        _msHelper->GetReorderedMsHandles());
    if (_settings.mfWeighting)
      _image_weight_initializer->InitializeMf(_imagingTable,
                                              *_imageWeightCache);
    _griddingTaskFactory = std::make_unique<GriddingTaskFactory>(
        *_msHelper, *_image_weight_initializer, _observationInfo, _l_shift,
        _m_shift, _imagingTable.EntryCount());
    _griddingTaskManager = GriddingTaskManager::Make(_settings);
    std::unique_ptr<PrimaryBeam> primaryBeam;
    for (size_t groupIndex = 0;
         groupIndex != _imagingTable.IndependentGroupCount(); ++groupIndex) {
      ImagingTable group = _imagingTable.GetIndependentGroup(groupIndex);
      runIndependentGroup(group, primaryBeam);
    }

    _griddingTaskManager.reset();
    _griddingTaskFactory.reset();
    _image_weight_initializer.reset();
    // Resetting the MsHelper will destroy its reordered ms handles and
    // thereby clear the temporary files if -save-reordered is not present.
    _msHelper.reset();

    if (_settings.channelsOut > 1) {
      for (PolarizationEnum pol : _settings.polarizations) {
        bool psfWasMade = (_settings.deconvolutionIterationCount > 0 ||
                           _settings.makePSF || _settings.makePSFOnly) &&
                          pol == *_settings.polarizations.begin();

        if (psfWasMade) {
          const size_t n_directions = _ddPsfCount ? _ddPsfCount : 1;
          Logger::Debug << "Using PSF index " << _infoForMFS.centralPsfIndex
                        << " for MF restoring beam fitting.\n";
          for (size_t direction = 0; direction != n_directions; ++direction) {
            std::optional<size_t> dd_psf_dir =
                _ddPsfCount ? std::optional<size_t>(direction) : std::nullopt;
            OutputChannelInfo mfs_info = _infoForMFS;
            ImageOperations::MakeMFSImage(_settings, _infoPerChannel, mfs_info,
                                          "psf", intervalIndex, pol,
                                          ImageFilenameType::Psf, dd_psf_dir);
            if (_settings.savePsfPb)
              ImageOperations::MakeMFSImage(
                  _settings, _infoPerChannel, mfs_info, "psf-pb", intervalIndex,
                  pol, ImageFilenameType::Psf, dd_psf_dir);
            if (direction == _infoForMFS.centralPsfIndex || _ddPsfCount == 0) {
              _infoForMFS = mfs_info;
            }
          }
        }
        if (griddingUsesATerms()) {
          ImageOperations::MakeMFSImage(_settings, _infoPerChannel, _infoForMFS,
                                        "", intervalIndex, pol,
                                        ImageFilenameType::Beam);
          // When faceting without beam, no beam images are stored, so skip
          // making the MFS beams in this case. In all other cases with beam, do
          // make them:
        } else if (usesBeam() && (_settings.facetSolutionFiles.empty() ||
                                  _settings.applyFacetBeam)) {
          // The (complex valued but Hermitian) Mueller matrices are stored with
          // 16 elements:
          constexpr size_t n_matrix_elements = 16;
          for (size_t beam_index = 0; beam_index != n_matrix_elements;
               ++beam_index) {
            ImageOperations::MakeMFSImage(
                _settings, _infoPerChannel, _infoForMFS,
                std::to_string(beam_index), intervalIndex, pol,
                ImageFilenameType::Beam);
          }
        }

        if (!(pol == Polarization::YX &&
              _settings.polarizations.count(Polarization::XY) != 0) &&
            !_settings.makePSFOnly) {
          if (_settings.isDirtySaved)
            ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                          _infoForMFS, "dirty", intervalIndex,
                                          pol, ImageFilenameType::Normal);
          if (_settings.deconvolutionIterationCount == 0) {
            ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                          _infoForMFS, "image", intervalIndex,
                                          pol, ImageFilenameType::Normal);
            if (usesBeam())
              ImageOperations::MakeMFSImage(
                  _settings, _infoPerChannel, _infoForMFS, "image-pb",
                  intervalIndex, pol, ImageFilenameType::Normal);
          } else {
            ImageOperations::MakeMFSImage(
                _settings, _infoPerChannel, _infoForMFS, "residual",
                intervalIndex, pol, ImageFilenameType::Normal);
            ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                          _infoForMFS, "model", intervalIndex,
                                          pol, ImageFilenameType::Normal);
            ImageOperations::RenderMFSImage(_settings, _infoForMFS,
                                            intervalIndex, pol, false, false);
            if (usesBeam()) {
              ImageOperations::MakeMFSImage(
                  _settings, _infoPerChannel, _infoForMFS, "residual-pb",
                  intervalIndex, pol, ImageFilenameType::Normal);
              ImageOperations::MakeMFSImage(
                  _settings, _infoPerChannel, _infoForMFS, "model-pb",
                  intervalIndex, pol, ImageFilenameType::Normal);
              ImageOperations::RenderMFSImage(_settings, _infoForMFS,
                                              intervalIndex, pol, false, true);
            }
          }
          if (Polarization::IsComplex(pol)) {
            if (_settings.isDirtySaved)
              ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                            _infoForMFS, "dirty", intervalIndex,
                                            pol, ImageFilenameType::Imaginary);
            if (_settings.deconvolutionIterationCount == 0) {
              ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                            _infoForMFS, "image", intervalIndex,
                                            pol, ImageFilenameType::Imaginary);
              if (usesBeam())
                ImageOperations::MakeMFSImage(
                    _settings, _infoPerChannel, _infoForMFS, "image-pb",
                    intervalIndex, pol, ImageFilenameType::Imaginary);
            } else {
              ImageOperations::MakeMFSImage(
                  _settings, _infoPerChannel, _infoForMFS, "residual",
                  intervalIndex, pol, ImageFilenameType::Imaginary);
              ImageOperations::MakeMFSImage(_settings, _infoPerChannel,
                                            _infoForMFS, "model", intervalIndex,
                                            pol, ImageFilenameType::Imaginary);
              ImageOperations::RenderMFSImage(_settings, _infoForMFS,
                                              intervalIndex, pol, true, false);
              if (usesBeam()) {
                ImageOperations::MakeMFSImage(
                    _settings, _infoPerChannel, _infoForMFS, "residual-pb",
                    intervalIndex, pol, ImageFilenameType::Imaginary);
                ImageOperations::MakeMFSImage(
                    _settings, _infoPerChannel, _infoForMFS, "model-pb",
                    intervalIndex, pol, ImageFilenameType::Imaginary);
                ImageOperations::RenderMFSImage(_settings, _infoForMFS,
                                                intervalIndex, pol, true, true);
              }
            }
          }
        }
      }
    }
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
  std::tie(_l_shift, _m_shift) = getLMShift();

  MSSelection fullSelection = _settings.GetMSSelection();

  for (size_t intervalIndex = 0; intervalIndex != _settings.intervalsOut;
       ++intervalIndex) {
    makeImagingTable(intervalIndex);
    _globalSelection = selectInterval(fullSelection, intervalIndex);

    _infoPerChannel.assign(
        _settings.channelsOut,
        OutputChannelInfo(std::max<size_t>(1, _facetCount), 0));

    _msHelper =
        std::make_unique<MsHelper>(_settings, _globalSelection, _msBands);
    if (_settings.reuseReorder)
      _msHelper->ReuseReorderedFiles(_imagingTable);
    else if (_settings.doReorder)
      _msHelper->PerformReordering(_imagingTable, true);

    _image_weight_initializer = std::make_unique<ImageWeightInitializer>(
        _settings, _globalSelection, _msBands,
        _msHelper->GetReorderedMsHandles());

    _griddingTaskFactory = std::make_unique<GriddingTaskFactory>(
        *_msHelper, *_image_weight_initializer, _observationInfo, _l_shift,
        _m_shift, _imagingTable.EntryCount());

    if (_facetCount != 0) {
      std::string prefix =
          ImageFilename::GetPrefix(_settings, _imagingTable[0].polarization,
                                   _imagingTable[0].outputChannelIndex,
                                   _imagingTable[0].outputIntervalIndex, false);
      aocommon::FitsReader reader(prefix + PredictModelFileSuffix());
      overrideImageSettings(reader);
      if (intervalIndex == 0) {
        facets = FacetReader::ReadFacets(
            _settings.facetRegionFilename, _settings.trimmedImageWidth,
            _settings.trimmedImageHeight, _settings.pixelScaleX,
            _settings.pixelScaleY, _observationInfo.phaseCentreRA,
            _observationInfo.phaseCentreDec, _l_shift, _m_shift,
            _settings.imagePadding, _settings.gridderType == GridderType::IDG,
            _settings.GetFeatherSize());

        // FIXME: raise warning if facets do not cover the entire image, see
        // AST-429
      }

      updateFacetsInImagingTable(facets, false);
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
      createWSCFitsWriter(groupTable.Front(), false, true));
  _modelImages.Initialize(modelWriter, _settings.polarizations.size(),
                          _settings.channelsOut, _facetCount,
                          _settings.prefixName + "-model");
  WSCFitsWriter writer(createWSCFitsWriter(groupTable.Front(), false, false));
  _residualImages.Initialize(writer, _settings.polarizations.size(),
                             _settings.channelsOut, _facetCount,
                             _settings.prefixName + "-residual");

  if (groupTable.Front().polarization == *_settings.polarizations.begin()) {
    const size_t image_count = _ddPsfCount ? _ddPsfCount : _facetCount;
    _psfImages.Initialize(writer, 1, groupTable.SquaredGroups().size(),
                          image_count, _settings.prefixName + "-psf");
    _scalarBeamImages.Initialize(writer, 1, groupTable.SquaredGroups().size(),
                                 image_count,
                                 _settings.prefixName + "-scalar-beam");
    _matrixBeamImages.Initialize(writer, 2, groupTable.SquaredGroups().size(),
                                 image_count,
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
  const bool doMakeDdPsf = doMakePSF && (_settings.ddPsfGridHeight > 1 ||
                                         _settings.ddPsfGridWidth > 1);

  if (doMakePSF) {
    ImagingTable::Groups facet_groups =
        groupTable.FacetGroups([&](const ImagingTableEntry& entry) {
          const bool is_first_polarization =
              entry.polarization == *_settings.polarizations.begin();
          return (entry.isDdPsf == doMakeDdPsf) && is_first_polarization;
        });
    if (_settings.reusePsf) {
      for (ImagingTable::Group& facet_group : facet_groups) {
        loadExistingImage(*facet_group.front(), true);
      }
    } else {
      for (ImagingTable::Group& facet_group : facet_groups) {
        ImagePsf(std::move(facet_group));
      }
      _isFirstInversionTask = false;
    }

    _griddingTaskManager->Finish();

    if (!doMakeDdPsf && !_settings.reusePsf)
      stitchFacets(groupTable, _psfImages, false, true, false);
  }

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
    const ImagingTable::Group& group,
    std::unique_ptr<PrimaryBeam>& primaryBeam) const {
  const ImagingTableEntry& tableEntry = *group.front();
  assert(tableEntry.facetIndex == 0);

  // Restore model to residual and save image
  const size_t channel_index = tableEntry.outputChannelIndex;
  const OutputChannelInfo& channel_info = _infoPerChannel[channel_index];

  PolarizationEnum curPol = tableEntry.polarization;
  for (size_t imageIter = 0; imageIter != tableEntry.imageCount; ++imageIter) {
    bool isImaginary = (imageIter == 1);
    WSCFitsWriter writer(createWSCFitsWriter(tableEntry, isImaginary, false));
    Image restoredImage(_settings.trimmedImageWidth,
                        _settings.trimmedImageHeight);
    _residualImages.Load(restoredImage.Data(), curPol, channel_index,
                         isImaginary);

    if (_settings.deconvolutionIterationCount != 0)
      writer.WriteImage("residual.fits", restoredImage);

    if (_settings.isUVImageSaved)
      saveUVImage(restoredImage, tableEntry, isImaginary, "uv");

    Image modelImage(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
    _modelImages.Load(modelImage.Data(), curPol, channel_index, isImaginary);
    double beamMaj = channel_info.beamMaj;
    double beamMin, beamPA;
    std::string beamStr;
    if (std::isfinite(beamMaj)) {
      beamMin = channel_info.beamMin;
      beamPA = channel_info.beamPA;
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
    schaapcommon::math::RestoreImage(
        restoredImage.Data(), modelImage.Data(), _settings.trimmedImageWidth,
        _settings.trimmedImageHeight, beamMaj, beamMin, beamPA,
        _settings.pixelScaleX, _settings.pixelScaleY);
    Logger::Info << "DONE\n";
    modelImage.Reset();

    Logger::Info << "Writing restored image... ";
    Logger::Info.Flush();
    writer.WriteImage("image.fits", restoredImage);
    Logger::Info << "DONE\n";
    restoredImage.Reset();

    const bool has_facet_solutions = !_settings.facetSolutionFiles.empty();
    ImageFilename imageName =
        ImageFilename(channel_index, tableEntry.outputIntervalIndex);

    // Checking for the _last_ polarization ensures that all images are
    // available before applying the primarybeam
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
      } else if (_settings.applyPrimaryBeam || _settings.applyFacetBeam ||
                 has_facet_solutions) {
        if (has_facet_solutions)
          primaryBeam->CorrectBeamForFacetGain(imageName, group, channel_info);
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
  Image ptr(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
  for (const ImagingTableEntry& entry : groupTable) {
    size_t ch = entry.outputChannelIndex;
    if (entry.polarization == Polarization::YX) {
      _residualImages.Load(ptr.Data(), Polarization::XY, ch, true);
      WSCFitsWriter writer(
          createWSCFitsWriter(entry, Polarization::XY, true, false));
      writer.WriteImage("first-residual.fits", ptr);
    } else {
      _residualImages.Load(ptr.Data(), entry.polarization, ch, false);
      WSCFitsWriter writer(createWSCFitsWriter(entry, false, false));
      writer.WriteImage("first-residual.fits", ptr);
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
          createWSCFitsWriter(entry, Polarization::XY, true, true));
      writer.WriteImage("model.fits", ptr);
    } else {
      _modelImages.Load(ptr.Data(), entry.polarization, ch, false);
      WSCFitsWriter writer(createWSCFitsWriter(entry, false, true));
      writer.WriteImage("model.fits", ptr);
    }
  }
}

void WSClean::partitionModelIntoFacets(const ImagingTable::Groups& facetGroups,
                                       bool isPredictOnly) {
  if (_facetCount != 0) {
    Logger::Info << "Clipping model image into facets...\n";
    // Allocate full image
    Image fullImage(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
    // Initialize FacetImage with properties of stitched image, always
    // stitch facets for 1 spectral term.
    schaapcommon::facets::FacetImage facetImage(
        _settings.trimmedImageWidth, _settings.trimmedImageHeight, 1);
    for (const ImagingTable::Group& facetGroup : facetGroups) {
      const size_t imageCount = facetGroup.front()->imageCount;
      _modelImages.Load(fullImage.Data(), facetGroup.front()->polarization,
                        facetGroup.front()->outputChannelIndex, false);
      for (size_t imageIndex = 0; imageIndex != imageCount; ++imageIndex) {
        partitionSingleGroup(facetGroup, imageIndex, _modelImages, fullImage,
                             facetImage, isPredictOnly);
      }
    }
  }
}

void WSClean::partitionSingleGroup(const ImagingTable::Group& facetGroup,
                                   size_t imageIndex,
                                   CachedImageSet& imageCache,
                                   const Image& fullImage,
                                   schaapcommon::facets::FacetImage& facetImage,
                                   bool isPredictOnly) {
  const bool isImaginary = (imageIndex == 1);
  for (const std::shared_ptr<ImagingTableEntry>& facetEntry : facetGroup) {
    const size_t channelIndex = facetEntry->outputChannelIndex;
    facetImage.SetFacet(*facetEntry->facet, true);
    facetImage.CopyToFacet({fullImage.Data()});
    if (!isPredictOnly) {
      if (_settings.UseFacetCorrections()) {
        // Before deconvolution, the images are converted to 'flat noise' image
        // by multiplying them with the sqrt(average squared correction). This
        // is undone here to convert the model image back to flat gain.
        const double factor =
            _infoPerChannel[channelIndex]
                .averageFacetCorrection[facetEntry->facetIndex]
                .GetStokesIValue();
        facetImage *= 1.0 / std::sqrt(factor);
      }
    }
    imageCache.StoreFacet(facetImage, facetEntry->polarization, channelIndex,
                          facetEntry->facetIndex, isImaginary);
  }
}

void WSClean::initializeModelImages(const ImagingTableEntry& entry,
                                    PolarizationEnum polarization,
                                    size_t maxFacetGroupIndex) {
  _modelImages.SetWSCFitsWriter(
      createWSCFitsWriter(entry, polarization, false, true));

  if (_settings.continuedRun) {
    readExistingModelImages(entry, polarization, maxFacetGroupIndex);
  } else {
    // Set model to zero: already done if this is YX of XY/YX imaging combi
    if (!(polarization == Polarization::YX &&
          _settings.polarizations.count(Polarization::XY) != 0)) {
      Image modelImage(_settings.trimmedImageWidth,
                       _settings.trimmedImageHeight, 0.0f);
      _modelImages.Store(modelImage, polarization, entry.outputChannelIndex,
                         false);
      if (Polarization::IsComplex(polarization))
        _modelImages.Store(modelImage, polarization, entry.outputChannelIndex,
                           true);
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

    aocommon::FitsReader reader(prefix + PredictModelFileSuffix());
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
      if (_settings.mfWeighting)
        _image_weight_initializer->InitializeMf(_imagingTable,
                                                *_imageWeightCache);
    }

    WSCFitsWriter writer(reader);
    _modelImages.SetWSCFitsWriter(writer);

    Image buffer(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
    reader.Read(buffer.Data());
    for (size_t j = 0;
         j != _settings.trimmedImageWidth * _settings.trimmedImageHeight; ++j) {
      if (!std::isfinite(buffer[j]))
        throw std::runtime_error(
            "The input image contains non-finite values -- can't predict "
            "from an image with non-finite values");
    }
    _modelImages.Store(buffer, polarization, entry.outputChannelIndex, i == 1);
  }
}

bool WSClean::overrideImageSettings(const aocommon::FitsReader& reader) {
  bool resetGridder = false;
  if (_settings.trimmedImageWidth == 0 && _settings.trimmedImageHeight == 0) {
    _settings.trimmedImageWidth = reader.ImageWidth();
    _settings.trimmedImageHeight = reader.ImageHeight();
    _settings.RecalculateDerivedDimensions();
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

  resetModelColumns(groupTable.FacetGroups());
  _griddingTaskManager->Start(getMaxNrMSProviders() *
                              (groupTable.MaxFacetGroupIndex() + 1));
  _predictingWatch.Start();

  for (size_t groupIndex = 0; groupIndex != groupTable.IndependentGroupCount();
       ++groupIndex) {
    const ImagingTable independentGroup =
        groupTable.GetIndependentGroup(groupIndex);

    // Initialize the model images before entering the gridding loop. This is
    // necessary because in IDG mode, predicting Stokes I will require all
    // model images to have been initialized.
    ImagingTable::Groups facetGroups = independentGroup.FacetGroups();
    for (ImagingTable::Group& facetGroup : facetGroups) {
      // For facet-based prediction: facetGroup contains only a list of facets
      // from the same (full) image. The meta data for the full model image can
      // be inferred from the first entry in the facetGroup table
      _modelImages.Initialize(
          createWSCFitsWriter(*facetGroup.front(), false, true),
          _settings.polarizations.size(), _settings.channelsOut, _facetCount,
          _settings.prefixName + "-model");

      readExistingModelImages(*facetGroup.front(),
                              facetGroup.front()->polarization,
                              groupTable.MaxFacetGroupIndex());
      partitionModelIntoFacets({facetGroup}, true);
    }

    for (ImagingTable::Group& facetGroup : facetGroups) {
      if (!gridPolarizationsAtOnce || facetGroup.front()->polarization ==
                                          *_settings.polarizations.begin()) {
        Predict(facetGroup);
      }
    }  // facet groups of different polarizations
  }    // independent groups (channels)

  _griddingTaskManager->Finish();
  _predictingWatch.Pause();

  Logger::Info << "Inversion: " << _inversionWatch.ToString()
               << ", prediction: " << _predictingWatch.ToString()
               << ", cleaning: " << _deconvolutionWatch.ToString() << '\n';
}

void WSClean::resetModelColumns(const ImagingTable::Groups& facet_groups) {
  assert(!facet_groups.empty());
  if (facet_groups.front().size() > 1) {
    for (const ImagingTable::Group& facet_group : facet_groups) {
      resetModelColumns(*facet_group.front());
    }
  }
}

void WSClean::resetModelColumns(const ImagingTableEntry& entry) {
  std::vector<MsListItem> ms_list = _msHelper->InitializeMsList(entry);
  for (auto& ms : ms_list) {
    ms.ms_description->GetProvider()->ResetModelColumn();
  }
}

void WSClean::runFirstInversions(ImagingTable& groupTable,
                                 std::unique_ptr<PrimaryBeam>& primaryBeam,
                                 bool requestPolarizationsAtOnce,
                                 bool parallelizePolarizations) {
  std::vector<ImagingTable::Groups> facetGroupsQueue;
  if (requestPolarizationsAtOnce) {
    facetGroupsQueue.emplace_back(
        groupTable.FacetGroups([&](const ImagingTableEntry& entry) {
          return !entry.isDdPsf &&
                 (entry.polarization == *_settings.polarizations.begin());
        }));
  } else if (parallelizePolarizations) {
    facetGroupsQueue.emplace_back(groupTable.FacetGroups());
  } else {
    // Only use parallelism for entries with the same polarization.
    for (aocommon::PolarizationEnum polarization : _settings.polarizations) {
      facetGroupsQueue.emplace_back(
          groupTable.FacetGroups([&](const ImagingTableEntry& entry) {
            return !entry.isDdPsf && (entry.polarization == polarization);
          }));
    }
  }
  for (ImagingTable::Groups& facetGroups : facetGroupsQueue) {
    for (ImagingTable::Group& facetGroup : facetGroups) {
      runFirstInversionGroup(facetGroup, primaryBeam);
    }
    _griddingTaskManager->Finish();
  }

  if (requestPolarizationsAtOnce) {
    groupTable.AssignGridDataFromPolarization(*_settings.polarizations.begin());
  }
  if (!_settings.reuseDirty)
    stitchFacets(groupTable, _residualImages, _settings.isDirtySaved, false,
                 false);
}

void WSClean::runFirstInversionGroup(
    ImagingTable::Group& facetGroup,
    std::unique_ptr<PrimaryBeam>& primaryBeam) {
  const bool doMakePSF = _settings.deconvolutionIterationCount > 0 ||
                         _settings.makePSF || _settings.makePSFOnly;

  for (const std::shared_ptr<ImagingTableEntry>& entry : facetGroup) {
    const bool isLastPol =
        entry->polarization == *_settings.polarizations.rbegin();
    if (isLastPol) {
      ImageFilename imageName =
          ImageFilename(entry->outputChannelIndex, entry->outputIntervalIndex);
      if (_settings.applyPrimaryBeam || _settings.applyFacetBeam ||
          !_settings.facetSolutionFiles.empty()) {
        std::vector<MsListItem> msList = _msHelper->InitializeMsList(*entry);
        std::shared_ptr<ImageWeights> weights =
            _image_weight_initializer->Initialize(*entry, msList,
                                                  *_imageWeightCache);
        primaryBeam = std::make_unique<PrimaryBeam>(_settings);
        for (MsListItem& item : msList)
          primaryBeam->AddMS(std::move(item.ms_description));
        primaryBeam->SetPhaseCentre(_observationInfo.phaseCentreRA,
                                    _observationInfo.phaseCentreDec, _l_shift,
                                    _m_shift);
        // Only generate beam images for facetIndex == 0 in facet group
        if (entry->facetIndex == 0) {
          if (_settings.applyPrimaryBeam || _settings.applyFacetBeam) {
            primaryBeam->MakeOrReuse(imageName, *entry, std::move(weights),
                                     _settings.fieldIds[0]);
          } else {
            // This condition happens when no beam is constructed but facet
            // solutions have been applied. Make a unitary beam to be able to
            // apply the facet corrections later on.
            primaryBeam->MakeUnitary(*entry, imageName, _settings);
          }
        }
      }
    }
  }

  if (_settings.reuseDirty) {
    for (const std::shared_ptr<ImagingTableEntry>& entry : facetGroup) {
      loadExistingImage(*entry, false);
    }
  } else {
    ImageMain(facetGroup, true, !doMakePSF);
  }
}

void WSClean::runMajorIterations(ImagingTable& groupTable,
                                 std::unique_ptr<PrimaryBeam>& primaryBeam,
                                 bool requestPolarizationsAtOnce,
                                 bool parallelizePolarizations) {
  std::unique_ptr<radler::WorkTable> deconvolution_table =
      groupTable.CreateDeconvolutionTable(_settings.deconvolutionChannelCount,
                                          _psfImages, _modelImages,
                                          _residualImages);

  ImagingTable tableWithoutDdPsf(
      groupTable,
      [](const ImagingTableEntry& entry) { return !entry.isDdPsf; });

  ImagingTable::Groups facetGroups = tableWithoutDdPsf.FacetGroups(
      [](const ImagingTableEntry& entry) { return true; });

  _deconvolution.emplace(_settings.GetRadlerSettings(),
                         std::move(deconvolution_table),
                         minTheoreticalBeamSize(tableWithoutDdPsf));

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
        writeFirstResidualImages(tableWithoutDdPsf);
      const bool isFinished = !reachedMajorThreshold;
      if (isFinished) {
        writeModelImages(tableWithoutDdPsf);
      }

      partitionModelIntoFacets(facetGroups, false);
      if (_settings.deconvolutionMGain != 1.0) {
        resetModelColumns(facetGroups);
        _griddingTaskManager->Start(
            getMaxNrMSProviders() *
            (tableWithoutDdPsf.MaxFacetGroupIndex() + 1));

        if (requestPolarizationsAtOnce) {
          // Only request one polarization for each facet/channel.
          // The gridder will grid all polarizations.
          const aocommon::PolarizationEnum first_polarization =
              *_settings.polarizations.begin();
          _predictingWatch.Start();
          for (const ImagingTable::Group& facetGroup : facetGroups) {
            if (facetGroup.front()->polarization == first_polarization) {
              Predict(facetGroup);
            }
          }
          _griddingTaskManager->Finish();
          _predictingWatch.Pause();

          _inversionWatch.Start();
          for (ImagingTable::Group& facetGroup : facetGroups) {
            if (facetGroup.front()->polarization == first_polarization) {
              ImageMain(facetGroup, false, false);
            }
          }
          _griddingTaskManager->Finish();
          _inversionWatch.Pause();
        } else if (parallelizePolarizations) {
          _predictingWatch.Start();
          for (const ImagingTable::Group& facetGroup : facetGroups) {
            Predict(facetGroup);
          }
          _griddingTaskManager->Finish();
          _predictingWatch.Pause();

          _inversionWatch.Start();
          for (ImagingTable::Group& facetGroup : facetGroups) {
            ImageMain(facetGroup, false, false);
          }
          _griddingTaskManager->Finish();
          _inversionWatch.Pause();
        } else {  // only parallelize channels
          _predictingWatch.Start();
          for (aocommon::PolarizationEnum polarization :
               _settings.polarizations) {
            for (const ImagingTable::Group& facetGroup : facetGroups) {
              if (facetGroup.front()->polarization == polarization) {
                Predict(facetGroup);
              }
            }
            _griddingTaskManager->Finish();
          }
          _predictingWatch.Pause();

          _inversionWatch.Start();
          for (aocommon::PolarizationEnum polarization :
               _settings.polarizations) {
            for (ImagingTable::Group& facetGroup : facetGroups) {
              if (facetGroup.front()->polarization == polarization) {
                ImageMain(facetGroup, false, false);
              }
            }
            _griddingTaskManager->Finish();
          }
          _inversionWatch.Pause();
        }
        stitchFacets(tableWithoutDdPsf, _residualImages, false, false, false);
      }

      ++_majorIterationNr;
    } while (reachedMajorThreshold);

    --_majorIterationNr;
    Logger::Info << _majorIterationNr << " major iterations were performed.\n";
  }

  if (_settings.applyFacetBeam && _settings.deconvolutionIterationCount != 0) {
    // The model facet images have already been corrected for their average gain
    // correction, so the full image can just be re-stitch from the facet images
    // to make the "fpb" facet corrected model images.
    stitchFacets(tableWithoutDdPsf, _modelImages, false, false, true);
  }

  for (const ImagingTable::Group& facetGroup : facetGroups) {
    saveRestoredImagesForGroup(facetGroup, primaryBeam);
  }

  if (_settings.saveSourceList) {
    std::unique_ptr<radler::WorkTable> deconvolution_table =
        tableWithoutDdPsf.CreateDeconvolutionTable(
            _settings.deconvolutionChannelCount, _psfImages, _modelImages,
            _residualImages);
    ComponentListWriter componentListWriter(_settings,
                                            std::move(deconvolution_table));
    componentListWriter.SaveSourceList(
        *_deconvolution, _observationInfo.phaseCentreRA,
        _observationInfo.phaseCentreDec, _l_shift, _m_shift);
    if (usesBeam()) {
      componentListWriter.SavePbCorrectedSourceList(
          *_deconvolution, _observationInfo.phaseCentreRA,
          _observationInfo.phaseCentreDec, _l_shift, _m_shift);
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
  saveUVImage(image, entry, _infoPerChannel[entry.outputChannelIndex],
              isImaginary, prefix);
}

void WSClean::saveUVImage(const Image& image, const ImagingTableEntry& entry,
                          const OutputChannelInfo& channel_info,
                          bool isImaginary, const std::string& prefix) const {
  Image realUV(_settings.trimmedImageWidth, _settings.trimmedImageHeight),
      imagUV(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
  schaapcommon::math::Resampler fft(
      _settings.trimmedImageWidth, _settings.trimmedImageHeight,
      _settings.trimmedImageWidth, _settings.trimmedImageHeight, 1);
  fft.SingleFT(image.Data(), realUV.Data(), imagUV.Data());
  // Factors of 2 involved: because of SingleFT()
  // (also one from the fact that normF excludes a factor of two?)
  realUV *=
      channel_info.normalizationFactor /
      sqrt(0.5 * _settings.trimmedImageWidth * _settings.trimmedImageHeight);
  imagUV *=
      channel_info.normalizationFactor /
      sqrt(0.5 * _settings.trimmedImageWidth * _settings.trimmedImageHeight);
  WSCFitsWriter writer(createWSCFitsWriter(entry, isImaginary, false));
  writer.WriteUV(prefix + "-real.fits", realUV.Data());
  writer.WriteUV(prefix + "-imag.fits", imagUV.Data());
}

void WSClean::ApplyFacetCorrectionForSingleChannel(
    const ImagingTable& squared_group, CachedImageSet& image_cache) {
  for (size_t facet_entry_index = 0;
       facet_entry_index != squared_group.FacetCount(); ++facet_entry_index) {
    const ImagingTable polarized_entries =
        squared_group.GetFacet(facet_entry_index);
    // Load the image data of this set of polarizations into facet_images
    std::vector<Image> facet_images;
    for (const ImagingTableEntry& facet_entry : polarized_entries) {
      const schaapcommon::facets::BoundingBox& box =
          facet_entry.facet->GetTrimmedBoundingBox();
      Image& facet_image = facet_images.emplace_back(box.Width(), box.Height());
      image_cache.LoadFacet(facet_image.Data(), facet_entry.polarization,
                            facet_entry.outputChannelIndex,
                            facet_entry.facetIndex, facet_entry.facet, false);
    }

    const size_t channel_index = polarized_entries.begin()->outputChannelIndex;
    const size_t facet_index = polarized_entries.begin()->facetIndex;
    const double stokes_i_sqrt =
        std::sqrt(_infoPerChannel[channel_index]
                      .averageFacetCorrection[facet_index]
                      .GetStokesIValue());
    aocommon::HMC4x4 matrix = _infoPerChannel[channel_index]
                                  .averageFacetCorrection[facet_index]
                                  .GetMatrixValue();
    if (matrix.Invert()) {
      if (facet_images.size() == 4) {
        std::array<Image*, 4> images{&facet_images[0], &facet_images[1],
                                     &facet_images[2], &facet_images[3]};
        math::CorrectImagesForMuellerMatrix(matrix * stokes_i_sqrt, images);
      } else if (facet_images.size() == 2) {
        std::array<Image*, 2> images{&facet_images[0], &facet_images[1]};
        math::CorrectDualImagesForMuellerMatrix(matrix * stokes_i_sqrt, images);
      } else {
        throw std::runtime_error(
            "A polarized facet correction was requested that is not "
            "implemented");
      }
    } else {
      for (Image& image : facet_images)
        image = std::numeric_limits<float>::quiet_NaN();
    }

    // Save the image data of this set of polarizations
    for (size_t polarization_index = 0;
         polarization_index != facet_images.size(); ++polarization_index) {
      Image& facet_image = facet_images[polarization_index];
      const ImagingTableEntry& facet_entry =
          polarized_entries[polarization_index];
      image_cache.StoreFacet(facet_image, facet_entry.polarization,
                             facet_entry.outputChannelIndex,
                             facet_entry.facetIndex, facet_entry.facet, false);
    }
  }
}

void WSClean::stitchFacets(const ImagingTable& table,
                           CachedImageSet& image_cache, bool write_dirty,
                           bool is_psf, bool is_facet_pb_model) {
  if (_facetCount != 0) {
    Logger::Info << "Stitching facets onto full image...\n";
    // Allocate full image
    Image full_image(_settings.trimmedImageWidth, _settings.trimmedImageHeight);

    std::unique_ptr<Image> weight_image;
    // This loop iterates over the output channels
    for (size_t sq_group = 0; sq_group != table.SquaredGroupCount();
         ++sq_group) {
      const ImagingTable squared_group = table.GetSquaredGroup(sq_group);
      const bool apply_correction =
          !is_psf && !is_facet_pb_model && _settings.UseFacetCorrections();
      // If a matrix correction is needed, it is done before stitching. If a
      // scalar correction is needed, it's done "during" stitching to save time.
      const bool apply_matrix = _settings.polarizations.size() != 1;
      if (apply_correction && apply_matrix) {
        ApplyFacetCorrectionForSingleChannel(squared_group, image_cache);
      }

      // Initialize FacetImage with properties of stitched image.
      // There's only one image, because WSClean uses only 1 spectral term.
      schaapcommon::facets::FacetImage facet_image(
          _settings.trimmedImageWidth, _settings.trimmedImageHeight, 1);
      for (const ImagingTable::Group& facet_group :
           squared_group.FacetGroups()) {
        // The PSF is only once imaged for all polarizations
        if (!is_psf || facet_group.front()->polarization ==
                           *_settings.polarizations.begin()) {
          const size_t image_count = facet_group.front()->imageCount;
          for (size_t image_index = 0; image_index != image_count;
               ++image_index) {
            stitchSingleGroup(facet_group, image_index, image_cache,
                              write_dirty, is_psf, full_image, weight_image,
                              facet_image, table.MaxFacetGroupIndex(),
                              apply_correction && !apply_matrix,
                              is_facet_pb_model);
          }
        }
      }
    }
  }
}

void WSClean::stitchSingleGroup(const ImagingTable::Group& facetGroup,
                                size_t imageIndex, CachedImageSet& imageCache,
                                bool writeDirty, bool isPSF, Image& fullImage,
                                std::unique_ptr<Image>& weight_image,
                                schaapcommon::facets::FacetImage& facetImage,
                                size_t maxFacetGroupIndex, bool apply_scalar,
                                bool is_facet_pb_model) {
  const size_t feather_size =
      is_facet_pb_model ? 0 : _settings.GetFeatherSize();
  const bool isImaginary = (imageIndex == 1);
  fullImage = 0.0f;
  if (feather_size != 0) {
    if (weight_image)
      *weight_image = 0.0f;
    else
      weight_image = std::make_unique<Image>(
          _settings.trimmedImageWidth, _settings.trimmedImageHeight, 0.0f);
  }
  for (const std::shared_ptr<ImagingTableEntry>& facetEntry : facetGroup) {
    facetImage.SetFacet(*facetEntry->facet, true);
    imageCache.LoadFacet(facetImage.Data(0), facetEntry->polarization,
                         facetEntry->outputChannelIndex, facetEntry->facetIndex,
                         facetEntry->facet, isImaginary);

    // Apply beam/h5parm correction to the image
    if (apply_scalar) {
      const size_t channelIndex = facetEntry->outputChannelIndex;
      const double factor = _infoPerChannel[channelIndex]
                                .averageFacetCorrection[facetEntry->facetIndex]
                                .GetStokesIValue();
      facetImage *= 1.0 / std::sqrt(factor);
    }

    // Place facet image on full image
    if (feather_size != 0) {
      aocommon::Image mask = facetImage.MakeMask();
      const schaapcommon::facets::Facet& facet = *facetEntry->facet;
      const bool needs_padding =
          facet.GetConvolutionBox() != facet.GetTrimmedBoundingBox();
      size_t left = 0;
      size_t top = 0;
      if (needs_padding) {
        const auto clipped_subtract = [](size_t a, size_t b) {
          return a > b ? a - b : 0;
        };
        left = clipped_subtract(facet.GetTrimmedBoundingBox().Min().x,
                                facet.GetConvolutionBox().Min().x);
        top = clipped_subtract(facet.GetTrimmedBoundingBox().Min().y,
                               facet.GetConvolutionBox().Min().y);
        const size_t right =
            clipped_subtract(facet.GetConvolutionBox().Max().x,
                             facet.GetTrimmedBoundingBox().Max().x);
        const size_t bottom =
            clipped_subtract(facet.GetConvolutionBox().Max().y,
                             facet.GetTrimmedBoundingBox().Max().y);
        mask = mask.Pad(left, top, right, bottom);
      }
      tophat_convolution::Convolve(mask, feather_size);
      if (needs_padding) {
        mask = mask.TrimBox(left, top, facet.GetTrimmedBoundingBox().Width(),
                            facet.GetTrimmedBoundingBox().Height());
      }
      facetImage.AddWithMask(fullImage.Data(), weight_image->Data(), mask, 0);
    } else {
      facetImage.AddToImage({fullImage.Data()});
    }
  }
  if (weight_image) {
    for (size_t i = 0; i != fullImage.Size(); ++i) {
      if ((*weight_image)[i] != 0.0)
        fullImage[i] /= (*weight_image)[i];
      else
        fullImage[i] = 0.0;
    }
  }
  if (writeDirty) {
    initializeModelImages(*facetGroup.front(), facetGroup.front()->polarization,
                          maxFacetGroupIndex);
    _residualImages.SetWSCFitsWriter(
        createWSCFitsWriter(*facetGroup.front(), false, false));
    WSCFitsWriter writer(
        createWSCFitsWriter(*facetGroup.front(), isImaginary, false));
    Logger::Info << "Writing dirty image...\n";
    writer.WriteImage("dirty.fits", fullImage);
  }
  if (is_facet_pb_model) {
    WSCFitsWriter writer(
        createWSCFitsWriter(*facetGroup.front(), isImaginary, false));
    Logger::Info << "Writing facet-corrected primary-beam model image...\n";
    writer.WriteImage("model-fpb.fits", fullImage);
  }

  if (isPSF) {
    const ImagingTableEntry& entry = *facetGroup.front();
    processFullPSF(fullImage, entry);
  }

  if (!is_facet_pb_model) {
    const size_t channelIndex = facetGroup.front()->outputChannelIndex;
    const PolarizationEnum polarization = facetGroup.front()->polarization;
    imageCache.Store(fullImage, polarization, channelIndex, isImaginary);
  }
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
    const bool isFirstPol = (p == *_settings.polarizations.begin());
    templateEntry.polarization = p;
    if (p == Polarization::XY)
      templateEntry.imageCount = 2;
    else if (p == Polarization::YX)
      templateEntry.imageCount = 0;
    else
      templateEntry.imageCount = 1;

    if (_ddPsfCount && isFirstPol) {
      ImagingTableEntry ddPsfTemplateEntry(templateEntry);
      ddPsfTemplateEntry.isDdPsf = true;
      addFacetsToImagingTable(ddPsfTemplateEntry, _ddPsfCount);
    }
    addFacetsToImagingTable(templateEntry, _facetCount);

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

void WSClean::addFacetsToImagingTable(ImagingTableEntry& templateEntry,
                                      const size_t facet_count) {
  // Create a single entry (with facetIndex == 0) when facets are not used.
  const size_t facet_entry_count = std::max(facet_count, std::size_t(1));
  for (size_t f = 0; f != facet_entry_count; ++f) {
    auto entry = std::make_unique<ImagingTableEntry>(templateEntry);
    entry->facetIndex = f;
    entry->facet.reset();  // updateFacetsInImagingTable will set the facet.
    _imagingTable.AddEntry(std::move(entry));
  }
  ++templateEntry.facetGroupIndex;
}

void WSClean::updateFacetsInImagingTable(
    const std::vector<std::shared_ptr<schaapcommon::facets::Facet>>& facets,
    bool updateDdPsfs) {
  for (ImagingTableEntry& entry : _imagingTable) {
    if (entry.isDdPsf != updateDdPsfs) continue;
    assert(entry.facetIndex < facets.size());
    entry.facet = facets[entry.facetIndex];
    // Calculate phase center delta for entry
    entry.centreShiftX = entry.facet->GetUntrimmedBoundingBox().Centre().x -
                         _settings.trimmedImageWidth / 2;
    entry.centreShiftY = entry.facet->GetUntrimmedBoundingBox().Centre().y -
                         _settings.trimmedImageHeight / 2;
  }
}

WSCFitsWriter WSClean::createWSCFitsWriter(
    const ImagingTableEntry& entry, const OutputChannelInfo& channel_info,
    bool isImaginary, bool isModel) const {
  return WSCFitsWriter(entry, isImaginary, _settings, _deconvolution,
                       _observationInfo, _l_shift, _m_shift, _majorIterationNr,
                       _commandLine, channel_info, isModel, _lastStartTime);
}

WSCFitsWriter WSClean::createWSCFitsWriter(const ImagingTableEntry& entry,
                                           bool isImaginary,
                                           bool isModel) const {
  return WSCFitsWriter(entry, isImaginary, _settings, _deconvolution,
                       _observationInfo, _l_shift, _m_shift, _majorIterationNr,
                       _commandLine, _infoPerChannel[entry.outputChannelIndex],
                       isModel, _lastStartTime);
}

WSCFitsWriter WSClean::createWSCFitsWriter(const ImagingTableEntry& entry,
                                           PolarizationEnum polarization,
                                           bool isImaginary,
                                           bool isModel) const {
  return WSCFitsWriter(
      entry, polarization, isImaginary, _settings, _deconvolution,
      _observationInfo, _l_shift, _m_shift, _majorIterationNr, _commandLine,
      _infoPerChannel[entry.outputChannelIndex], isModel, _lastStartTime);
}

}  // namespace wsclean
