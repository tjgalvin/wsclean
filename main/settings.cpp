#include "settings.h"

#include "../io/facetreader.h"

#include <aocommon/logger.h>

#include <schaapcommon/h5parm/h5parm.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/TableRecord.h>

#include <sstream>

using aocommon::Logger;

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

  if (facetRegionFilename.empty()) {
    if (!facetSolutionFiles.empty())
      throw std::runtime_error(
          "A facet solution file can only be specified in conjunction with a "
          "facet regions file. Either remove -apply-facet-solutions from the "
          "command line, or specify a facet regions file with -facet-regions.");
    if (applyFacetBeam)
      throw std::runtime_error(
          "A facet beam can only applied if a facet regions file is specified. "
          "Either remove -apply-facet-beam from the command line, or specify a "
          "regions file with -facet-regions.");
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
      if (applyFacetBeam && !facetSolutionFiles.empty()) {
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

      const std::size_t nfacets =
          FacetReader::ReadFacets(facetRegionFilename).size();
      for (const std::string& facetSolutionFile : facetSolutionFiles) {
        schaapcommon::h5parm::H5Parm h5parm =
            schaapcommon::h5parm::H5Parm(facetSolutionFile);
        const size_t nsources = h5parm.GetNumSources();
        if (nsources != nfacets) {
          throw std::runtime_error(
              "Number of source directions in one of the h5 facet solution "
              "files does not match the number of facets in the facet "
              "definition file.");
        }
      }
    }
  }

  if (useIDG) {
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

  if (gridWithBeam && !useIDG)
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
  if (!useIDG && !atermConfigFilename.empty())
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
      spectralFittingMode == SpectralFittingMode::NoFitting)
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

  if (saveSourceList &&
      (polarizations.size() != 1 ||
       (*polarizations.begin()) != aocommon::Polarization::StokesI))
    throw std::runtime_error(
        "Saving a source list currently only works for Stokes I imaging.");

  if (saveSourceList && deconvolutionIterationCount == 0)
    throw std::runtime_error("A source list cannot be saved without cleaning.");

  if (!forcedSpectrumFilename.empty() &&
      (spectralFittingMode != SpectralFittingMode::LogPolynomial ||
       spectralFittingTerms != 2))
    throw std::runtime_error(
        "When using forced spectrum mode, currently it is required to fit "
        "logarithmic polynomials (i.e. spectral index + further terms). This "
        "implies you have to specify -fit-spectral-log-pol 2.");

  if (parallelGridding != 1 &&
      (applyFacetBeam || !facetSolutionFiles.empty()) &&
      !schaapcommon::h5parm::H5Parm::IsThreadSafe()) {
    throw std::runtime_error(
        "Parallel gridding in combination with a facet beam or facet solutions,"
        " requires an HDF5 library that supports multi-threading.");
  }

  checkPolarizations();
}

void Settings::checkPolarizations() const {
  bool hasXY = polarizations.count(aocommon::Polarization::XY) != 0;
  bool hasYX = polarizations.count(aocommon::Polarization::YX) != 0;
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
    if (joinedPolarizationDeconvolution)
      throw std::runtime_error(
          "You have requested spectral fitting, but you are joining multiple "
          "polarizations. This is not supported. You probably want to turn off "
          "the joining of polarizations (leave out -join-polarizations).");
    if (!joinedFrequencyDeconvolution)
      throw std::runtime_error(
          "You have requested spectral fitting, but you are not joining "
          "channels. This is not possible: you probably want to turn channel "
          "joining on (add -join-channels).");
  }

  if (autoDeconvolutionThreshold && autoMask) {
    if (autoDeconvolutionThresholdSigma >= autoMaskSigma)
      throw std::runtime_error(
          "The auto-masking threshold was smaller or equal to the "
          "auto-threshold. This does not make sense. Did you accidentally "
          "reverse the auto-mask and auto-threshold values?");
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

  // When using IDG with aterms, a PSF must be made, because the beam
  // image is created during the PSF imaging stage.
  if (useIDG && (!atermConfigFilename.empty() || gridWithBeam)) {
    makePSF = true;
  }

  // When using IDG, polarizations should always be joined
  if (useIDG) {
    joinedPolarizationDeconvolution = true;
  }

  if (mode == ImagingMode || mode == PredictMode) {
    RecalculatePaddedDimensions(verbose);
    doReorder = determineReorder();
    dataColumnName = determineDataColumn(verbose);
  }
}

void Settings::RecalculatePaddedDimensions(bool verbose) {
  paddedImageWidth = (size_t)ceil(trimmedImageWidth * imagePadding);
  paddedImageHeight = (size_t)ceil(trimmedImageHeight * imagePadding);
  // Make the width and height divisable by four.
  paddedImageWidth += (4 - (paddedImageWidth % 4)) % 4;
  paddedImageHeight += (4 - (paddedImageHeight % 4)) % 4;
  if (trimmedImageWidth != 0 && trimmedImageHeight != 0) {
    if (verbose)
      Logger::Debug << "Using image size of " << trimmedImageWidth << " x "
                    << trimmedImageHeight << ", padded to " << paddedImageWidth
                    << " x " << paddedImageHeight << ".\n";
  }
}

DeconvolutionSettings Settings::GetDeconvolutionSettings() const {
  DeconvolutionSettings deconvolutionSettings;

  deconvolutionSettings.trimmedImageWidth = trimmedImageWidth;
  deconvolutionSettings.trimmedImageHeight = trimmedImageHeight;
  deconvolutionSettings.channelsOut = channelsOut;
  deconvolutionSettings.pixelScaleX = pixelScaleX;
  deconvolutionSettings.pixelScaleY = pixelScaleY;
  deconvolutionSettings.threadCount = threadCount;
  deconvolutionSettings.prefixName = prefixName;

  deconvolutionSettings.linkedPolarizations = linkedPolarizations;
  deconvolutionSettings.parallelDeconvolutionMaxSize =
      parallelDeconvolutionMaxSize;
  deconvolutionSettings.parallelDeconvolutionMaxThreads =
      parallelDeconvolutionMaxThreads;

  deconvolutionSettings.deconvolutionThreshold = deconvolutionThreshold;
  deconvolutionSettings.deconvolutionGain = deconvolutionGain;
  deconvolutionSettings.deconvolutionMGain = deconvolutionMGain;
  deconvolutionSettings.autoDeconvolutionThreshold = autoDeconvolutionThreshold;
  deconvolutionSettings.autoMask = autoMask;
  deconvolutionSettings.autoDeconvolutionThresholdSigma =
      autoDeconvolutionThresholdSigma;
  deconvolutionSettings.autoMaskSigma = autoMaskSigma;
  deconvolutionSettings.localRMS = localRMS;
  deconvolutionSettings.localRMSWindow = localRMSWindow;
  deconvolutionSettings.localRMSMethod = localRMSMethod;
  deconvolutionSettings.saveSourceList = saveSourceList;
  deconvolutionSettings.deconvolutionIterationCount =
      deconvolutionIterationCount;
  deconvolutionSettings.majorIterationCount = majorIterationCount;
  deconvolutionSettings.allowNegativeComponents = allowNegativeComponents;
  deconvolutionSettings.stopOnNegativeComponents = stopOnNegativeComponents;
  deconvolutionSettings.useMultiscale = useMultiscale;
  deconvolutionSettings.useSubMinorOptimization = useSubMinorOptimization;
  deconvolutionSettings.squaredJoins = squaredJoins;
  deconvolutionSettings.spectralCorrectionFrequency =
      spectralCorrectionFrequency;
  deconvolutionSettings.spectralCorrection = spectralCorrection;
  deconvolutionSettings.multiscaleFastSubMinorLoop = multiscaleFastSubMinorLoop;
  deconvolutionSettings.multiscaleGain = multiscaleGain;
  deconvolutionSettings.multiscaleDeconvolutionScaleBias =
      multiscaleDeconvolutionScaleBias;
  deconvolutionSettings.multiscaleMaxScales = multiscaleMaxScales;
  deconvolutionSettings.multiscaleConvolutionPadding =
      multiscaleConvolutionPadding;
  deconvolutionSettings.multiscaleScaleList.assign(multiscaleScaleList.begin(),
                                                   multiscaleScaleList.end());
  deconvolutionSettings.multiscaleShapeFunction = multiscaleShapeFunction;
  deconvolutionSettings.deconvolutionBorderRatio = deconvolutionBorderRatio;
  deconvolutionSettings.fitsDeconvolutionMask = fitsDeconvolutionMask;
  deconvolutionSettings.casaDeconvolutionMask = casaDeconvolutionMask;
  deconvolutionSettings.horizonMask = horizonMask;
  deconvolutionSettings.horizonMaskDistance = horizonMaskDistance;
  deconvolutionSettings.localRMSImage = localRMSImage;
  deconvolutionSettings.pythonDeconvolutionFilename =
      pythonDeconvolutionFilename;
  deconvolutionSettings.useMoreSaneDeconvolution = useMoreSaneDeconvolution;
  deconvolutionSettings.useIUWTDeconvolution = useIUWTDeconvolution;
  deconvolutionSettings.iuwtSNRTest = iuwtSNRTest;
  deconvolutionSettings.moreSaneLocation = moreSaneLocation;
  deconvolutionSettings.moreSaneArgs = moreSaneArgs;
  deconvolutionSettings.moreSaneSigmaLevels.assign(moreSaneSigmaLevels.begin(),
                                                   moreSaneSigmaLevels.end());
  deconvolutionSettings.spectralFittingMode = spectralFittingMode;
  deconvolutionSettings.spectralFittingTerms = spectralFittingTerms;
  deconvolutionSettings.forcedSpectrumFilename = forcedSpectrumFilename;
  deconvolutionSettings.deconvolutionChannelCount = deconvolutionChannelCount;

  return deconvolutionSettings;
}

bool Settings::determineReorder() const {
  return ((channelsOut != 1) || (polarizations.size() >= 4) ||
          (deconvolutionMGain != 1.0) ||
          (baselineDependentAveragingInWavelengths != 0.0) || simulateNoise ||
          forceReorder) &&
         !forceNoReorder;
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
