#include "settings.h"

#include "../io/logger.h"

#include <sstream>

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
          "Using IDG with multiple polarizations is only possible in joined "
          "polarization mode: use -join-polarizations or -link-polarizations.");
    if (trimmedImageWidth != trimmedImageHeight)
      throw std::runtime_error(
          "IDG can not yet make rectangular images -- this will be implemented "
          "at a later time.");
    if (parallelGridding != 1)
      throw std::runtime_error(
          "Parallel gridding can not be combined with IDG");
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
  if (useDifferentialLofarBeam && !(gridWithBeam || applyPrimaryBeam))
    throw std::runtime_error(
        "Differential beam correction was requested, but no beam correction is "
        "applied. Use either IDG and grid with the beam, or apply the average "
        "beam.");

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
      spectralFittingMode == NoSpectralFitting)
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
