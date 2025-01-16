#include "settings.h"

#include <sstream>

#include <aocommon/logger.h>

#include <schaapcommon/h5parm/h5parm.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/TableRecord.h>

#include "../io/facetreader.h"

using aocommon::Logger;

namespace wsclean {

namespace {
/**
 * Determines the number of images in one dimension of an image grid.
 * @param image_size Size of the image, in some dimension.
 * @param max_grid_size Maximum size of a grid cell, in the same dimension.
 *        Zero indicates using a single grid cell.
 * @return The number of grid cells.
 */
size_t GetNumCells(const size_t image_size, const size_t max_grid_size) {
  size_t images = 1;
  if (max_grid_size > 0) {
    images = (image_size + max_grid_size - 1) / max_grid_size;
    // Since an image typically has an interesting object at its center,
    // ensure that the number of grid cells is odd. That object then resides in
    // one grid cell instead of two.
    images |= 1;
  }
  return images;
}
}  // namespace

void Settings::Validate() const {
  if (mode == ImagingMode) {
    if (trimmedImageWidth == 0 && trimmedImageHeight == 0)
      throw std::runtime_error("Image size has not been set.");

    if (trimmedImageWidth == 0 || trimmedImageHeight == 0)
      throw std::runtime_error(
          "Invalid image size given: one of the dimensions was zero.");

    if (pixelScaleX == 0.0 && pixelScaleY == 0.0)
      throw std::runtime_error("Pixel scale has not been set.");

    if (pixelScaleX == 0.0 || pixelScaleY == 0.0)
      throw std::runtime_error(
          "Invalid pixel scale given: one direction was set to zero.");
  } else if (mode == PredictMode) {
    if (joinedFrequencyDeconvolution)
      throw std::runtime_error(
          "Joined frequency deconvolution specified for prediction: "
          "prediction doesn't clean, parameter invalid");
    if (joinedPolarizationDeconvolution)
      throw std::runtime_error(
          "Joined polarization deconvolution specified for prediction: "
          "prediction doesn't clean, parameter invalid");
  }

  if (threadCount == 0)
    throw std::runtime_error("A thread count of zero (-j 0) is not valid");

  // antialiasingKernelSize should be odd
  if (antialiasingKernelSize % 2 == 0) {
    std::stringstream s;
    s << "Bad anti-aliasing kernel size given of " << antialiasingKernelSize
      << ". The kernel size has to be odd.";
    throw std::runtime_error(s.str());
  }

  if (visibilityReadMode != VisibilityReadMode::kFull) {
    if (!UseFacetCorrections()) {
      throw std::runtime_error(
          "-scalar-visibilities and -diagonal-visibilities must be combined "
          "with -apply-facet-solutions or -apply-facet-beam");
    }
    if (visibilityReadMode == VisibilityReadMode::kDiagonal &&
        (polarizations.size() != 1 ||
         *polarizations.begin() != aocommon::PolarizationEnum::StokesI)) {
      throw std::runtime_error(
          "-diagonal-visibilities can only be used when making Stokes I "
          "images");
    }
  }

  if (facetRegionFilename.empty() && ddPsfGridWidth == 1 &&
      ddPsfGridHeight == 1) {
    if (!facetSolutionFiles.empty())
      throw std::runtime_error(
          "A facet solution file can only be specified in conjunction with a "
          "facet regions file. Either remove -apply-facet-solutions from the "
          "command line, or specify a facet regions file with -facet-regions,"
          " or specify a grid of dd-psfs with -dd-psf-grid.");
    if (applyFacetBeam)
      throw std::runtime_error(
          "A facet beam can only applied if a facet regions file is specified, "
          "or if a dd-psf is used. "
          "Either remove -apply-facet-beam from the command line, or specify a "
          "regions file with -facet-regions, or specify a grid of dd-psfs with "
          "-dd-psf-grid.");
  } else {
    if (polarizations.size() > 1) {
      // -join-polarizations required in order to write the pb.fits images
      // in PrimaryBeam::CorrectImages
      if (applyFacetBeam && !joinedPolarizationDeconvolution &&
          deconvolutionIterationCount != 0) {
        throw std::runtime_error(
            "Can not apply the facet beam of multiple polarizations "
            "independently. Add -join-polarizations to the command line to "
            "apply the facet beam for multiple polarizations");
      }
      // This condition might become a bit more specific once xx,yy polarization
      // correction for h5 AND beam are implemented
      if (applyFacetBeam && !facetSolutionFiles.empty() &&
          !aocommon::Polarization::HasFullStokesPolarization(polarizations)) {
        throw std::runtime_error(
            "Applying h5parm solutions AND beam correction on multiple "
            "polarizations is not yet supported.");
      }
    }

    if (!facetSolutionFiles.empty()) {
      if (facetSolutionFiles.size() != 1 &&
          facetSolutionFiles.size() != filenames.size()) {
        throw std::runtime_error(
            "Incorrect number of facet solution files provided. The number of "
            "facet solution files should be either 1 or match the number of "
            "input measurement sets.");
      }

      const std::size_t nfacets = FacetReader::CountFacets(facetRegionFilename);
      for (const std::string& facetSolutionFile : facetSolutionFiles) {
        schaapcommon::h5parm::H5Parm h5parm =
            schaapcommon::h5parm::H5Parm(facetSolutionFile);
        const size_t nsources = h5parm.GetNumSources();
        if (nsources != nfacets) {
          const std::string message =
              "The number of source directions (" + std::to_string(nsources) +
              ") in h5 facet solution file " + facetSolutionFile +
              " does not match the number of facets (" +
              std::to_string(nfacets) + ") in the facet region file.";
          if (solutionDirectionsCheck)
            throw std::runtime_error(message);
          else
            Logger::Warn << message << '\n';
        }
      }
    }
  }

  if (facetRegionFilename.empty() && featherSize && *featherSize != 0) {
    throw std::runtime_error(
        "Parameter -feather-size was specified without enabling facetting.");
  }

  if (gridderType == GridderType::IDG) {
    const bool stokesIOnly =
        polarizations.size() == 1 &&
        *polarizations.begin() == aocommon::Polarization::StokesI;
    const bool allStokes =
        aocommon::Polarization::HasFullStokesPolarization(polarizations) &&
        polarizations.size() == 4;
    if (!allStokes && !stokesIOnly) {
      throw std::runtime_error(
          "When using IDG, it is only possible to either image Stokes I or to "
          "image all 4 Stokes polarizations: use -pol i or -pol iquv.");
    }
    if (polarizations.size() > 1 && !joinedPolarizationDeconvolution &&
        deconvolutionIterationCount != 0)
      throw std::runtime_error(
          "Deconvolving IDG images with multiple polarizations is only "
          "possible in joined "
          "polarization mode: use -join-polarizations or -link-polarizations.");
    if (trimmedImageWidth != trimmedImageHeight)
      throw std::runtime_error(
          "IDG can not yet make rectangular images -- this will be implemented "
          "at a later time.");
    if (parallelGridding != 1)
      throw std::runtime_error(
          "Parallel gridding can not be combined with IDG");
    if (applyPrimaryBeam)
      throw std::runtime_error(
          "IDG currently does not support -apply-primary-beam. Use a-term "
          "correction with -grid-with-beam instead.");
    if (applyFacetBeam)
      throw std::runtime_error(
          "IDG cannot apply facet based beam corrections. Remove facet related "
          "command line arguments and use -grid-with-beam "
          "instead.");
    if (!facetSolutionFiles.empty())
      throw std::runtime_error(
          "IDG cannot apply facet based direction dependent corrections. "
          "Remove -apply-facet-solution from the command line instruction.");

    if (baselineDependentAveragingInWavelengths != 0.0) {
      throw std::runtime_error(
          "IDG cannot be combined with (internally computed) "
          "baseline-dependent averaging. Please remove baseline-averaging "
          "option from your command.");
    }

    for (const auto& filename : filenames) {
      casacore::MeasurementSet ms(filename);
      const std::string& bda_factors = "BDA_FACTORS";
      const bool has_bda = ms.keywordSet().isDefined(bda_factors) &&
                           (ms.keywordSet().asTable(bda_factors).nrow() > 0);
      if (has_bda) {
        throw std::runtime_error(
            "IDG cannot be combined with the baseline-dependently averaged "
            "measurement set " +
            filename);
      }
    }
  }

  if (UseMpi()) {
    if (!masterDoesWork && nMpiNodes <= 1) {
      throw std::runtime_error(
          "Master was told not to work, but no other workers available.");
    }
    if (channelToNode.size() != channelsOut) {
      throw std::runtime_error(
          "Channel-to-node map should have exactly channels-out entries.");
    }
    for (size_t node : channelToNode) {
      if (!masterDoesWork && node == 0) {
        throw std::runtime_error(
            "Master was told not to work, but the channel-to-node map assigned "
            "a channel to node 0.");
      }
      if (node >= nMpiNodes) {
        throw std::runtime_error("Invalid node number in channel-to-node map.");
      }
    }
    const std::size_t kMessageLimit{std::numeric_limits<std::int32_t>::max()};
    if (maxMpiMessageSize > kMessageLimit) {
      throw std::runtime_error("MPI message size larger than " +
                               std::to_string(kMessageLimit) +
                               " is not supported.");
    }
  }

  if (gridWithBeam && gridderType != GridderType::IDG)
    throw std::runtime_error(
        "Can't grid with the beam without IDG: specify '-use-idg' to use IDG.");
  if (gridWithBeam && applyPrimaryBeam)
    throw std::runtime_error(
        "Can't simultaneously grid with the beam and apply the average beam: "
        "use either one.");
  if (gridWithBeam && !atermConfigFilename.empty())
    throw std::runtime_error(
        "Use of an aterm config file can't be combined with -grid-with-beam: "
        "add the beam to your aterm config and remove -grid-with-beam from the "
        "command line");
  if (gridderType != GridderType::IDG && !atermConfigFilename.empty())
    throw std::runtime_error(
        "Use of an aterm config file required IDG enabled: add -use-idg");

  if (baselineDependentAveragingInWavelengths != 0.0) {
    if (forceNoReorder)
      throw std::runtime_error(
          "Baseline dependent averaging can not be performed without "
          "reordering.");
    if (modelUpdateRequired)
      throw std::runtime_error(
          "Baseline dependent averaging can not update the model column (yet) "
          "-- you have to add -no-update-model-required.");
  }

  if (simulateNoise) {
    if (forceNoReorder)
      throw std::runtime_error(
          "Noise simulation can not be performed without reordering.");
  }

  if (channelsOut == 0)
    throw std::runtime_error(
        "You have specified 0 output channels -- at least one output channel "
        "is required.");

  if (joinedFrequencyDeconvolution && channelsOut == 1)
    throw std::runtime_error(
        "Joined frequency deconvolution was requested, but only one output "
        "channel is being requested. Did you forget -channels-out?");

  if (forceReorder && forceNoReorder)
    throw std::runtime_error(
        "Can not both force reordering and force not reordering!");

  if (deconvolutionChannelCount != 0 &&
      deconvolutionChannelCount != channelsOut &&
      spectralFittingMode ==
          schaapcommon::fitters::SpectralFittingMode::kNoFitting)
    throw std::runtime_error(
        "You have requested to deconvolve with a decreased number of channels "
        "(-deconvolution-channels), but you have not enabled spectral fitting. "
        "You should specify an interpolation function by enabling spectral "
        "fitting in order to interpolate the deconvolved channels back to the "
        "full number of channels. The most useful and common spectral fitting "
        "function is -fit-spectral-pol.");

  if (savePsfPb && !(applyPrimaryBeam || gridWithBeam))
    throw std::runtime_error(
        "You can not save the primary-beam corrected PSF without enabling "
        "primary beam correction: add -apply-primary-beam to your commandline "
        "or use IDG to apply the beam.");

  if (saveSourceList) {
    if (polarizations.size() != 1 ||
        (*polarizations.begin() != aocommon::Polarization::StokesI &&
         *polarizations.begin() != aocommon::Polarization::XX &&
         *polarizations.begin() != aocommon::Polarization::YY &&
         *polarizations.begin() != aocommon::Polarization::LL &&
         *polarizations.begin() != aocommon::Polarization::RR)) {
      throw std::runtime_error(
          "Saving a source list currently only works for Stokes I or pseudo "
          "Stokes I (XX, YY, LL or RR) imaging.");
    } else if (!IsSpectralFittingEnabled() && channelsOut > 1 &&
               joinedFrequencyDeconvolution) {
      throw std::runtime_error(
          "Saving a source list with multiple channels requires specifying a "
          "fitting method");
    }
  }

  if (saveSourceList && deconvolutionIterationCount == 0)
    throw std::runtime_error("A source list cannot be saved without cleaning.");

  if (!forcedSpectrumFilename.empty() &&
      spectralFittingMode !=
          schaapcommon::fitters::SpectralFittingMode::kLogPolynomial)
    throw std::runtime_error(
        "When using forced spectrum mode, it is required to fit logarithmic"
        "polynomials (i.e. spectral index + further terms). This "
        "implies you have to specify -fit-spectral-log-pol <N>, with N the"
        "number of terms.");

  if (parallelGridding != 1 && UseFacetCorrections() &&
      !schaapcommon::h5parm::H5Parm::IsThreadSafe()) {
    throw std::runtime_error(
        "Parallel gridding in combination with a facet beam or facet solutions,"
        " requires an HDF5 library that supports multi-threading.");
  }

  if (UseFacetCorrections() && polarizations.size() > 1 &&
      !joinedPolarizationDeconvolution) {
    throw std::runtime_error(
        "Applying beam or solutions per facet for multiple polarizations "
        "required joining or linking the polarizations.");
  }

  if (reuseDirty && (gridWithBeam || !atermConfigFilename.empty())) {
    throw std::runtime_error(
        "Reusing dirty image and beam/aterm corrections"
        " can not be combined, because the average beam is"
        " computed when the dirty image is made.");
  }

  checkPolarizations();
}

void Settings::checkPolarizations() const {
  const bool hasXY = polarizations.count(aocommon::Polarization::XY) != 0;
  const bool hasYX = polarizations.count(aocommon::Polarization::YX) != 0;
  if (joinedPolarizationDeconvolution) {
    if (polarizations.size() == 1)
      throw std::runtime_error(
          "Joined/linked polarization deconvolution requested, but only one "
          "polarization is being imaged. Specify multiple polarizations, or do "
          "not request to join the polarizations.");
  } else {
    if ((hasXY || hasYX) && deconvolutionIterationCount != 0)
      throw std::runtime_error(
          "You are imaging XY and/or YX polarizations and have enabled "
          "cleaning (niter!=0). This is not possible -- you have to specify "
          "'-join-polarizations' or disable cleaning.");
  }

  for (aocommon::PolarizationEnum p : linkedPolarizations) {
    if (polarizations.count(p) == 0) {
      std::ostringstream str;
      str << "Linked polarization cleaning was requested for polarization "
          << aocommon::Polarization::TypeToFullString(p)
          << ", but this polarization is not imaged.";
      throw std::runtime_error(str.str());
    }
  }

  if ((hasXY && !hasYX) || (!hasXY && hasYX))
    throw std::runtime_error(
        "You are imaging only one of the XY or YX polarizations. This is not "
        "possible -- you have to specify both XY and YX polarizations (the "
        "output of imaging both polarizations will be the XY and imaginary XY "
        "images).");
  if (IsSpectralFittingEnabled()) {
    if (!joinedFrequencyDeconvolution)
      throw std::runtime_error(
          "You have requested spectral fitting, but you are not joining "
          "channels. This is not possible: you probably want to turn channel "
          "joining on (add -join-channels).");
  }

  if (autoDeconvolutionThresholdSigma && autoMaskSigma) {
    if (*autoDeconvolutionThresholdSigma >= *autoMaskSigma)
      throw std::runtime_error(
          "The auto-masking threshold was smaller or equal to the "
          "auto-threshold. This does not make sense. Did you accidentally "
          "reverse the auto-mask and auto-threshold values?");
  }
  if (hasXY && hasYX && gridderType != GridderType::WStacking) {
    throw std::runtime_error(
        "Combined XY/YX imaging is not possible with gridders other than the "
        "w-stacking gridder. Either add '-gridder wstacking' to the command "
        "line, or image a different set of polarizations (e.g. iquv).");
  }
  if (UseFacetCorrections()) {
    const bool is_i =
        (polarizations == std::set{aocommon::PolarizationEnum::StokesI});
    const bool is_xx_yy =
        (polarizations == std::set{aocommon::PolarizationEnum::XX,
                                   aocommon::PolarizationEnum::YY});
    const bool is_iquv =
        aocommon::Polarization::HasFullStokesPolarization(polarizations);
    if (!(is_i || is_xx_yy || is_iquv)) {
      throw std::runtime_error(
          "Facet imaging with facet corrections is only possible for Stokes I, "
          "XX+YY or Stokes IQUV imaging");
    }
  }
}

void Settings::Propagate(bool verbose) {
  if (verbose) logImportantSettings();

  if (trimmedImageWidth % 2 != 0) {
    ++trimmedImageWidth;
    Logger::Warn << "Image width is not divisable by two: changing width to "
                 << trimmedImageWidth << '\n';
  }
  if (trimmedImageHeight % 2 != 0) {
    ++trimmedImageHeight;
    Logger::Warn << "Image height is not divisable by two: changing height to "
                 << trimmedImageHeight << '\n';
  }

  if (parallelDeconvolutionMaxThreads == 0) {
    parallelDeconvolutionMaxThreads = threadCount;
  }

  // When using IDG, polarizations should always be joined
  if (gridderType == GridderType::IDG) {
    joinedPolarizationDeconvolution = true;
  }

  if (mode == ImagingMode || mode == PredictMode) {
    RecalculateDerivedDimensions(verbose);
    doReorder = determineReorder();
    dataColumnName = determineDataColumn(verbose);
  }
}

void Settings::RecalculateDerivedDimensions(bool verbose) {
  paddedImageWidth =
      static_cast<size_t>(ceil(trimmedImageWidth * imagePadding));
  paddedImageHeight =
      static_cast<size_t>(ceil(trimmedImageHeight * imagePadding));
  // Make the width and height divisable by four.
  paddedImageWidth += (4 - (paddedImageWidth % 4)) % 4;
  paddedImageHeight += (4 - (paddedImageHeight % 4)) % 4;
  if (trimmedImageWidth != 0 && trimmedImageHeight != 0) {
    if (verbose)
      Logger::Debug << "Using image size of " << trimmedImageWidth << " x "
                    << trimmedImageHeight << ", padded to " << paddedImageWidth
                    << " x " << paddedImageHeight << ".\n";

    if (!makePSFOnly || parallelDeconvolutionMaxSize > 0) {
      parallelDeconvolutionGridWidth =
          GetNumCells(trimmedImageWidth, parallelDeconvolutionMaxSize);
      parallelDeconvolutionGridHeight =
          GetNumCells(trimmedImageHeight, parallelDeconvolutionMaxSize);
      if (ddPsfGridWidth > parallelDeconvolutionGridWidth ||
          ddPsfGridHeight > parallelDeconvolutionGridHeight) {
        Logger::Warn
            << "Warning: The DD PSF grid (" << ddPsfGridWidth << "x"
            << ddPsfGridHeight
            << ") has more cells than parallel deconvolution grid ("
            << parallelDeconvolutionGridWidth << "x"
            << parallelDeconvolutionGridHeight
            << ") in at least one dimension. Reducing the DD PSF grid to ";
        ddPsfGridWidth =
            std::min(ddPsfGridWidth, parallelDeconvolutionGridWidth);
        ddPsfGridHeight =
            std::min(ddPsfGridHeight, parallelDeconvolutionGridHeight);
        Logger::Warn << ddPsfGridWidth << "x" << ddPsfGridHeight << ".\n";
      }
    }
  }
}

radler::Settings Settings::GetRadlerSettings() const {
  radler::Settings radler_settings;

  radler_settings.trimmed_image_width = trimmedImageWidth;
  radler_settings.trimmed_image_height = trimmedImageHeight;
  radler_settings.channels_out = channelsOut;
  radler_settings.pixel_scale.x = pixelScaleX;
  radler_settings.pixel_scale.y = pixelScaleY;
  radler_settings.thread_count = threadCount;
  radler_settings.prefix_name = prefixName;
  radler_settings.linked_polarizations = linkedPolarizations;
  radler_settings.parallel.grid_width = parallelDeconvolutionGridWidth;
  radler_settings.parallel.grid_height = parallelDeconvolutionGridHeight;
  radler_settings.parallel.max_threads = parallelDeconvolutionMaxThreads;
  radler_settings.auto_threshold_sigma = autoDeconvolutionThresholdSigma;
  radler_settings.absolute_threshold =
      absoluteDeconvolutionThreshold.value_or(0.0);
  radler_settings.auto_mask_sigma = autoMaskSigma;
  radler_settings.absolute_auto_mask_threshold = absoluteAutoMaskThreshold;
  radler_settings.minor_loop_gain = deconvolutionGain;
  radler_settings.major_loop_gain = deconvolutionMGain;
  radler_settings.local_rms.method = localRMSMethod;
  radler_settings.local_rms.strength = localRMSStrength;
  radler_settings.local_rms.window = localRMSWindow;
  radler_settings.local_rms.image = localRMSImage;
  radler_settings.save_source_list = saveSourceList;
  radler_settings.minor_iteration_count = deconvolutionIterationCount;
  radler_settings.major_iteration_count = majorIterationCount;
  radler_settings.allow_negative_components = allowNegativeComponents;
  radler_settings.stop_on_negative_components = stopOnNegativeComponents;
  radler_settings.squared_joins = squaredJoins;
  radler_settings.spectral_correction_frequency = spectralCorrectionFrequency;
  radler_settings.spectral_correction = spectralCorrection;
  radler_settings.border_ratio = deconvolutionBorderRatio;
  radler_settings.fits_mask = fitsDeconvolutionMask;
  radler_settings.casa_mask = casaDeconvolutionMask;
  if (horizonMask) {
    radler_settings.horizon_mask_distance = horizonMaskDistance;
  }
  if (!forcedSpectrumFilename.empty())
    radler_settings.spectral_fitting.mode =
        schaapcommon::fitters::SpectralFittingMode::kForcedTerms;
  else
    radler_settings.spectral_fitting.mode = spectralFittingMode;
  radler_settings.spectral_fitting.terms = spectralFittingTerms;
  radler_settings.spectral_fitting.forced_filename = forcedSpectrumFilename;
  radler_settings.algorithm_type = algorithmType;

  switch (algorithmType) {
    case radler::AlgorithmType::kAdaptiveScalePixel:
    case radler::AlgorithmType::kMultiscale:
      radler_settings.multiscale.fast_sub_minor_loop =
          multiscaleFastSubMinorLoop;
      radler_settings.multiscale.sub_minor_loop_gain = multiscaleGain;
      radler_settings.multiscale.scale_bias = multiscaleDeconvolutionScaleBias;
      radler_settings.multiscale.max_scales = multiscaleMaxScales;
      radler_settings.multiscale.convolution_padding =
          multiscaleConvolutionPadding;
      radler_settings.multiscale.scale_list.assign(multiscaleScaleList.begin(),
                                                   multiscaleScaleList.end());
      radler_settings.multiscale.shape = multiscaleShapeFunction;
      break;
    case radler::AlgorithmType::kIuwt:
      // IUWT has no algorithm-specific settings
      break;
    case radler::AlgorithmType::kMoreSane:
      radler_settings.more_sane.location = moreSaneLocation;
      radler_settings.more_sane.arguments = moreSaneArgs;
      radler_settings.more_sane.sigma_levels.assign(moreSaneSigmaLevels.begin(),
                                                    moreSaneSigmaLevels.end());
      break;
    case radler::AlgorithmType::kPython:
      radler_settings.python.filename = pythonDeconvolutionFilename;
      break;
    case radler::AlgorithmType::kGenericClean:
      radler_settings.generic.use_sub_minor_optimization =
          useSubMinorOptimization;
      break;
  }

  return radler_settings;
}

aocommon::PolarizationEnum Settings::GetProviderPolarization(
    aocommon::PolarizationEnum entry_polarization) const {
  const bool xx_and_yy =
      polarizations ==
      std::set{aocommon::PolarizationEnum::XX, aocommon::PolarizationEnum::YY};
  const bool stokes_i =
      polarizations.size() == 1 &&
      (*polarizations.begin()) == aocommon::Polarization::StokesI;
  if (gridderType == GridderType::IDG) {
    if (stokes_i) {
      if ((ddPsfGridWidth > 1 || ddPsfGridHeight > 1) && gridWithBeam) {
        return aocommon::Polarization::StokesI;
      } else {
        return aocommon::Polarization::DiagonalInstrumental;
      }
    } else {
      return aocommon::Polarization::Instrumental;
    }
  } else if (visibilityReadMode == VisibilityReadMode::kDiagonal) {
    return aocommon::Polarization::DiagonalInstrumental;
  } else if (visibilityReadMode == VisibilityReadMode::kScalar) {
    return entry_polarization;
  } else if ((xx_and_yy || stokes_i) && UseFacetCorrections()) {
    return aocommon::Polarization::Instrumental;
  } else {
    bool requires_instrumental = false;
    for (aocommon::PolarizationEnum p : polarizations) {
      if (aocommon::Polarization::IsStokes(p) &&
          p != aocommon::PolarizationEnum::StokesI) {
        requires_instrumental = UseFacetCorrections();
        break;
      }
    }
    if (requires_instrumental) {
      return aocommon::Polarization::Instrumental;
    } else {
      return entry_polarization;
    }
  }
}

bool Settings::determineReorder() const {
  const bool prefer_reordering =
      (channelsOut != 1) || (polarizations.size() >= 4) ||
      (deconvolutionMGain != 1.0) ||
      (baselineDependentAveragingInWavelengths != 0.0) || simulateNoise ||
      !facetRegionFilename.empty();
  return (prefer_reordering || forceReorder) && !forceNoReorder;
}

std::string Settings::determineDataColumn(bool verbose) const {
  // If no column specified, determine column to use
  if (mode == PredictMode) return "DATA";
  std::string col = dataColumnName;
  if (col.empty()) {
    casacore::MeasurementSet ms(filenames.front());
    bool hasCorrected = ms.tableDesc().isColumn("CORRECTED_DATA");
    if (hasCorrected) {
      if (verbose)
        Logger::Info
            << "First measurement set has corrected data: tasks will be "
               "applied on the corrected data column.\n";
      col = "CORRECTED_DATA";
    } else {
      if (verbose)
        Logger::Info
            << "No corrected data in first measurement set: tasks will "
               "be applied on the data column.\n";
      col = "DATA";
    }
  }
  return col;
}

void Settings::logImportantSettings() const {
  Logger::Debug << "Number of threads selected: " << threadCount << '\n';
}

}  // namespace wsclean
