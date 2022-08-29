#include "commandline.h"
#include "wsclean.h"

#include <wscversion.h>

#include "../structures/numberlist.h"

#include <aocommon/fits/fitswriter.h>
#include <aocommon/logger.h>
#include <aocommon/radeccoord.h>
#include <aocommon/units/angle.h>
#include <aocommon/units/fluxdensity.h>

#include <schaapcommon/fitters/spectralfitter.h>
#include <schaapcommon/h5parm/jonesparameters.h>

#include <boost/algorithm/string.hpp>

#include <iostream>
#include <optional>
#include <string>
#include <sstream>

using aocommon::Logger;
using aocommon::units::Angle;
using aocommon::units::FluxDensity;

namespace {
void IncArgi(int& argi, int argc) {
  ++argi;
  if (argi >= argc)
    throw std::runtime_error("Unexpected end of command line arguments");
}

void Deprecated(bool isSlave, const std::string& param,
                const std::string& replacement) {
  if (!isSlave)
    Logger::Warn << "!!! WARNING: Parameter \'-" << param
                 << "\' is deprecated and will be removed in a future version "
                    "of WSClean.\n"
                 << "!!!          Use parameter \'-" << replacement
                 << "\' instead.\n";
}

void PrintHeader() {
  Logger::Info << "\n"
                  "WSClean version " WSCLEAN_VERSION_STR
                  " (" WSCLEAN_VERSION_DATE
                  ")\n"
                  "This software package is released under the GPL version 3.\n"
                  "Author: André Offringa (offringa@gmail.com).\n\n";
#ifndef NDEBUG
  Logger::Info
      << "\n"
         "WARNING: Symbol NDEBUG was not defined; this WSClean version was\n"
         "compiled as a DEBUG version. This can seriously affect "
         "performance!\n\n";
#endif
}

void PrintHelp() {
  std::cout << R"(Syntax: wsclean [options] <input-ms> [<2nd-ms> [..]]
Will create cleaned images of the input ms(es).
If multiple mses are specified, they need to be phase-rotated to the same point on the sky.

Options can be:

  ** GENERAL OPTIONS **
-version
   Print WSClean's version and exit.
-j <threads>
   Specify number of computing threads to use, i.e., number of cpu cores that will be used.
   Default: use all cpu cores.
-parallel-gridding <n>
   Will execute multiple gridders simultaneously. This can make things faster in certain cases,
   but will increase memory usage.
-parallel-reordering <n>
   Process the reordering with multipliple threads.
-no-work-on-master
   In MPI runs, do not use the master for gridding. This may be useful if the
   resources such as memory of the master are limited.
-mem <percentage>
   Limit memory usage to the given fraction of the total system memory. This is an approximate value.
   Default: 100.
-abs-mem <memory limit>
   Like -mem, but this specifies a fixed amount of memory in gigabytes.
-verbose (or -v)
   Increase verbosity of output.
-log-time
   Add date and time to each line in the output.
-quiet
   Do not output anything but errors.
-reorder
-no-reorder
   Force or disable reordering of Measurement Set. This can be faster when the measurement set needs to
   be iterated several times, such as with many major iterations or in channel imaging mode.
   Default: only reorder when in channel imaging mode.
-temp-dir <directory>
   Set the temporary directory used when reordering files. Default: same directory as input measurement set.
-update-model-required (default), and
-no-update-model-required
   These two options specify whether the model data column is required to
   contain valid model data after imaging. It can save time to not update
   the model data column.
-no-dirty
   Do not save the dirty image.
-save-first-residual
   Save the residual after the first iteration.
-save-weights
   Save the gridded weights in the a fits file named <image-prefix>-weights.fits.
-save-uv
   Save the gridded uv plane, i.e., the FFT of the residual image. The UV plane is complex, hence
   two images will be output: <prefix>-uv-real.fits and <prefix>-uv-imag.fits.
-reuse-psf <prefix>
   Load the psf(s) from the given prefix and skip the inversion for the psf image.
-reuse-dirty <prefix>
   Load the dirty from the given prefix and skip the inversion for the dirty image.
-apply-primary-beam
   Calculate and apply the primary beam and save images for the Jones components, with weighting identical to the
   weighting as used by the imager. Only available for instruments
   supported by EveryBeam.
-reuse-primary-beam
   If a primary beam image exists on disk, reuse those images.
-use-differential-lofar-beam
   Assume the visibilities have already been beam-corrected for the reference direction.
   By default, WSClean will use the information in the measurement set to determine
   if the differential beam should be applied for obtaining proper flux levels.
-primary-beam-limit <limit>
   Level at which to trim the beam when performing image-based beam
   correction,. Default: 0.005.
-scalar-beam
   In the case of Stokes I imaging, this will take the average of
   1/XX and 1/YY instead of the inverted Mueller matrix.
-mwa-path <path>
   Set path where to find the MWA beam file(s).
-save-psf-pb
   When applying beam correction, also save the primary-beam corrected PSF image.
-pb-grid-size <npixel>
   Specify the grid size in number of pixels at which to evaluate the primary beam.
   Typically, the primary beam is calculated at a coarse resolution grid
   and interpolated, to reduce the time spent in evaluating the beam.
   This parameter controls the resolution of the grid at which to evaluate
   the primary beam. Default: 32.
-dd-psf-grid <width> <height>
   This parameter enables direction-dependent psfs.
   Select the grid size (number of cells in both directions).
   Default: 1 1 (no direction-dependent psfs).
-beam-model
   Specify the beam model, only relevant for SKA and LOFAR. Available models are Hamaker, Lobes, OskarDipole, OskarSphericalWave.
   Input is case insensitive. Default is Hamaker for LOFAR and
   OskarSphericalWave for SKA.
-beam-mode
   [DEBUGGING ONLY] Manually specify the beam mode. Only relevant for simulated SKA measurement sets.
   Available modes are array_factor, element and full.
   Input is case insensitive. Default is full.
-beam-normalisation-mode
    [DEBUGGING ONLY] Manually specify the normalisation of the beam. Only relevant for simulated SKA measurement sets.
    Available modes are none, preapplied, full, and amplitude. Default is preapplied.
-dry-run
   Parses the command line and quits afterwards. No imaging is done.

  ** WEIGHTING OPTIONS **
-weight <weightmode>
   Weightmode can be: natural, uniform, briggs. Default: uniform. When using Briggs' weighting,
   add the robustness parameter, like: "-weight briggs 0.5".
-super-weight <factor>
   Increase the weight gridding box size, similar to Casa's superuniform weighting scheme. Default: 1.0
   The factor can be rational and can be less than one for subpixel weighting.
-mf-weighting
   In spectral mode, calculate the weights as if the image was made using MF. This makes sure that the sum of
   channel images equals the MF weights. Otherwise, the channel image will become a bit more naturally weighted.
   This is only relevant for weighting modes that require gridding (i.e., Uniform, Briggs').
   Default: off, unless -join-channels is specified.
-no-mf-weighting
   Opposite of -ms-weighting; can be used to turn off MF weighting in -join-channels mode.
-weighting-rank-filter <level>
   Filter the weights and set high weights to the local mean. The level parameter specifies
   the filter level; any value larger than level*localmean will be set to level*localmean.
-weighting-rank-filter-size <size>
   Set size of weighting rank filter. Default: 16.
-taper-gaussian <beamsize>
   Taper the weights with a Gaussian function. This will reduce the contribution of long baselines.
   The beamsize is by default in asec, but a unit can be specified ("2amin").
-taper-tukey <lambda>
   Taper the outer weights with a Tukey transition. Lambda specifies the size of the transition; use in
   combination with -maxuv-l.
-taper-inner-tukey <lambda>
   Taper the weights with a Tukey transition. Lambda specifies the size of the transition; use in
   combination with -minuv-l.
-taper-edge <lambda>
   Taper the weights with a rectangle, to keep a space of lambda between the edge and gridded visibilities.
-taper-edge-tukey <lambda>
   Taper the edge weights with a Tukey window. Lambda is the size of the Tukey transition. When -taper-edge
   is also specified, the Tukey transition starts inside the inner rectangle.
-use-weights-as-taper
   Will not use visibility weights when determining the imaging weights.
   This has the effect that e.g. uniform weighting can be modified by increasing
   the visibility weight of certain baselines. Without this option, uniform imaging
   weights absorb the visibility weight to make the weighting truly uniform.
-store-imaging-weights
   Will store the imaging weights in a column named 'IMAGING_WEIGHT_SPECTRUM'.

  ** INVERSION OPTIONS **
-name <image-prefix>
   Use image-prefix as prefix for output files. Default is 'wsclean'.
-size <width> <height>
   Set the output image size in number of pixels (without padding).
-padding <factor>
   Pad images by the given factor during inversion to avoid aliasing. Default: 1.2 (=20%).
-scale <pixel-scale>
   Scale of a pixel. Default unit is degrees, but can be specificied, e.g. -scale 20asec. Default: 0.01deg.
-predict
   Only perform a single prediction for an existing image. Doesn't do any imaging or cleaning.
   The input images should have the same name as the model output images would have in normal imaging mode.
-continue
   Will continue an earlier WSClean run. Earlier model images will be read and model visibilities will be
   subtracted to create the first dirty residual. CS should have been used in the earlier run, and model data
   should have been written to the measurement set for this to work. Default: off.
-subtract-model
   Subtract the model from the data column in the first iteration. This can be used to reimage
   an already cleaned image, e.g. at a different resolution.
-gridder <type>
   Set gridder type: direct-ft, idg, wgridder, tuned-wgridder, or wstacking.
-channels-out <count>
   Splits the bandwidth and makes count nr. of images. Default: 1.
-shift <ra> <dec>
   Shift the phase centre to the given location. The shift is along
   the tangential plane.
-gap-channel-division
   In case of irregular frequency spacing, this option can be used to not try and split channels
   to make the output channel bandwidth similar, but instead to split largest gaps first.
-channel-division-frequencies <list>
   Split the bandwidth at the specified frequencies (in Hz) before the normal bandwidth
   division is performed. This can e.g. be useful for imaging multiple bands with irregular
   number of channels.
-nwlayers <nwlayers>
   Number of w-layers to use. Default: minimum suggested #w-layers for first MS.
-nwlayers-factor <factor>
   Use automatic calculation of the number of w-layers, but multiple that number by
   the given factor. This can e.g. be useful for increasing w-accuracy.
-nwlayers-for-size <width> <height>
   Use the minimum suggested w-layers for an image of the given size. Can e.g. be used to increase
   accuracy when predicting small part of full image.
-no-small-inversion and -small-inversion
   Perform inversion at the Nyquist resolution and upscale the image to the requested image size afterwards.
   This speeds up inversion considerably, but makes aliasing slightly worse. This effect is
   in most cases <1%. Default: on.
-grid-mode <"nn", "kb" or "rect">
   Kernel and mode used for gridding: kb = Kaiser-Bessel (default with 7 pixels), nn = nearest
   neighbour (no kernel), more options: rect, kb-no-sinc, gaus, bn. Default: kb.
-kernel-size <size>
   Gridding antialiasing kernel size. Default: 7.
-oversampling <factor>
   Oversampling factor used during gridding. Default: 63.
-make-psf
   Always make the psf, even when no cleaning is performed.
-make-psf-only
   Only make the psf, no images are made.
-visibility-weighting-mode [normal/squared/unit]
   Specify visibility weighting modi. Affects how the weights (normally) stored in
   WEIGHT_SPECTRUM column are applied. Useful for estimating e.g. EoR power spectra errors.
-baseline-averaging <size-in-wavelengths>
   Enable baseline-dependent averaging. The specified size is in number of wavelengths (i.e., uvw-units). One way
   to calculate this is with <baseline in nr. of lambdas> * 2pi * <acceptable integration in s> / (24*60*60).
-simulate-noise <stddev-in-jy>
   Will replace every visibility by a Gaussian distributed value with given standard deviation before imaging.
-simulate-baseline-noise <filename>
   Like -simulate-noise, but the stddevs are provided per baseline, in a text file
   with antenna1 and antenna2 indices and the stddev per line, separated by spaces, e.g. "0 1 3.14".
-idg-mode [cpu/gpu/hybrid]
   Sets the IDG mode. Default: cpu. Hybrid is recommended when a GPU is available.
-wgridder-accuracy <value>
   Set the w-gridding accuracy. Default: 1e-4
   Useful range: 1e-2 to 1e-6

  ** A-TERM GRIDDING **
-aterm-config <filename>
   Specify a parameter set describing how a-terms should be applied. Please refer to the documentation for
   details of the configuration file format. Applying a-terms is only possible when IDG is enabled.
-grid-with-beam
   Apply a-terms to correct for the primary beam. This is only possible when IDG is enabled.
-beam-aterm-update <seconds>
   Set the ATerm update time in seconds. The default is every 300 seconds.
   It also sets the interval over which to calculate the primary beam when using
   -apply-primary-beam when not gridding with the beam.
-aterm-kernel-size <double>
   Kernel size reserved for aterms by IDG.
-apply-facet-solutions <path-to-file> <name1[,name2]>
   Apply solutions from the provided (h5) file per facet when gridding facet based images.
   Provided file is assumed to be in H5Parm format.
   Filename is followed by a comma separated list of strings specifying which sol tabs from the provided H5Parm file are used.
-apply-facet-beam
   Apply beam gains to facet center when gridding facet based images or direction dependent psfs
-facet-beam-update <seconds>
   Set the facet beam update time in seconds. The default is every 120 seconds.
-save-aterms
   Output a fits file for every aterm update, containing the applied image for every station.

  ** DATA SELECTION OPTIONS **
-pol <list>
   Default: 'I'. Possible values: XX, XY, YX, YY, I, Q, U, V, RR, RL, LR or LL (case insensitive).
   It is allowed but not necessary to separate with commas, e.g.: 'xx,xy,yx,yy'.   Two or four polarizations can be joinedly cleaned (see '-joinpolarizations'), but
   this is not the default. I, Q, U and V polarizations will be directly calculated from
   the visibilities, which might require correction to get to real IQUV values. The
   'xy' polarization will output both a real and an imaginary image, which allows calculating
   true Stokes polarizations for those telescopes.
-interval <start-index> <end-index>
   Only image the given time interval. Indices specify the timesteps, end index is exclusive.
   Default: image all time steps.
-intervals-out <count>
   Number of intervals to image inside the selected global interval. Default: 1
-even-timesteps
   Only select even timesteps. Can be used together with -odd-timesteps to determine noise values.
-odd-timesteps
   Only select odd timesteps.
-channel-range <start-channel> <end-channel>
   Only image the given channel range. Indices specify channel indices, end index is exclusive.
   Default: image all channels.
-field <list>
   Image the given field id(s). A comma-separated list of field ids can be provided. When multiple
   fields are given, all fields should have the same phase centre. Specifying '-field all' will image
   all fields in the measurement set. Default: first field (id 0).
-spws <list>
   Selects only the spws given in the list. list should be a comma-separated list of integers. Default: all spws.
-data-column <columnname>
   Default: CORRECTED_DATA if it exists, otherwise DATA will be used.
-maxuvw-m <meters>
-minuvw-m <meters>
   Set the min/max baseline distance in meters.
-maxuv-l <lambda>
-minuv-l <lambda>
   Set the min/max uv distance in lambda.
-maxw <percentage>
   Do not grid visibilities with a w-value higher than the given percentage of the max w, to save speed.
   Default: grid everything

  ** DECONVOLUTION OPTIONS **
-niter <niter>
   Maximum number of clean iterations to perform. Default: 0 (=no cleaning)
-nmiter <nmiter>
   Maximum number of major clean (inversion/prediction) iterations. Default: 20.   A value of 0 means no limit.
-threshold <threshold>
   Stopping clean thresholding in Jy. Default: 0.0
-auto-threshold <sigma>
   Estimate noise level using a robust estimator and stop at sigma x stddev.
-auto-mask <sigma>
   Construct a mask from found components and when a threshold of sigma is reached, continue
   cleaning with the mask down to the normal threshold.
-local-rms
   Instead of using a single RMS for auto thresholding/masking, use a spatially varying
   RMS image.
-local-rms-window
   Size of window for creating the RMS background map, in number of PSFs. Default: 25 psfs.
-local-rms-method
   Either 'rms' (default, uses sliding window RMS) or 'rms-with-min' (use max(window rms, 0.3 x window min)).
-gain <gain>
   Cleaning gain: Ratio of peak that will be subtracted in each iteration. Default: 0.1
-mgain <gain>
   Cleaning gain for major iterations: Ratio of peak that will be subtracted in each major
   iteration. To use major iterations, 0.85 is a good value. Default: 1.0
-join-polarizations
   Perform deconvolution by searching for peaks in the sum of squares of the polarizations,
   but subtract components from the individual images. Only possible when imaging two or four Stokes
   or linear parameters. Default: off.
-link-polarizations <pollist>
   Links all polarizations to be cleaned from the given list: components are found in the
   given list, but cleaned from all polarizations.
-facet-regions <facets.reg>
   Split the image into facets using the facet regions defined in  the facets.reg file. Default: off.
-join-channels
   Perform deconvolution by searching for peaks in the MF image,
but subtract components from individual channels.
   This will turn on mf-weighting by default. Default: off.
-spectral-correction <reffreq> <term list>
   Enable correction of the given spectral function inside deconvolution.
   This can e.g. avoid downweighting higher frequencies because of
   reduced flux density. 1st term is total flux, 2nd is si, 3rd curvature, etc.
   Example: -spectral-correction 150e6 83.084,-0.699,-0.110
-no-fast-subminor
   Do not use the subminor loop optimization during (non-multiscale) cleaning. Default: use the optimization.
-multiscale
   Clean on different scales. This is a new algorithm. Default: off.
   This parameter invokes the optimized multiscale algorithm published by Offringa & Smirnov (2017).
-multiscale-scale-bias
   Parameter to prevent cleaning small scales in the large-scale iterations. A lower
   bias will give more focus to larger scales. Default: 0.6
-multiscale-max-scales <n>
   Set the maximum number of scales that WSClean should use in multiscale cleaning.
   Only relevant when -multiscale-scales is not set. Default: unlimited.
-multiscale-scales <comma-separated list of sizes in pixels>
   Sets a list of scales to use in multi-scale cleaning. If unset, WSClean will select the delta
   (zero) scale, scales starting at four times the synthesized PSF, and increase by a factor of
   two until the maximum scale is reached or the maximum number of scales is reached.
   Example: -multiscale-scales 0,5,12.5
-multiscale-shape <shape>
   Sets the shape function used during multi-scale clean. Either 'tapered-quadratic' (default) or 'gaussian'.
-multiscale-gain <gain>
   Size of step made in the subminor loop of multi-scale. Default currently 0.2, but shows sign of instability.
   A value of 0.1 might be more stable.
-multiscale-convolution-padding <padding>
   Size of zero-padding for convolutions during the multi-scale cleaning. Default: 1.1
-no-multiscale-fast-subminor
   Disable the 'fast subminor loop' optimization, that will only search a part of the
   image during the multi-scale subminor loop. The optimization is on by default.
-python-deconvolution <filename>
   Run a custom deconvolution algorithm written in Python. See manual
   for the interface.
-iuwt
   Use the IUWT deconvolution algorithm.
-iuwt-snr-test / -no-iuwt-snr-test
   Stop (/do not stop) IUWT when the SNR decreases. This might help limitting divergence, but can
   occasionally also stop the algorithm too early. Default: no SNR test.
-moresane-ext <location>
   Use the MoreSane deconvolution algorithm, installed at the specified location.
-moresane-arg <arguments>
   Pass the specified arguments to moresane. Note that multiple parameters have to be
   enclosed in quotes.
-moresane-sl <sl1,sl2,...>
   MoreSane --sigmalevel setting for each major loop iteration. Useful to start at high
   levels and go down with subsequent loops, e.g. 20,10,5
-save-source-list
   Saves the found clean components as a BBS/DP3 text sky model. This parameter
   enables Gaussian shapes during multi-scale cleaning (-multiscale-shape gaussian).
-clean-border <percentage>
   Set the border size in which no cleaning is performed, in percentage of the width/height of the image.
   With an image size of 1000 and clean border of 1%, each border is 10 pixels. Default: 0%
-fits-mask <mask>
   Use the specified fits-file as mask during cleaning.
-casa-mask <mask>
   Use the specified CASA mask as mask during cleaning.
-horizon-mask <distance>
   Use a mask that avoids cleaning emission beyond the horizon. Distance is an angle (e.g. "5deg")
   that (when positive) decreases the size of the mask to stay further away from the horizon.
-no-negative
   Do not allow negative components during cleaning. Not the default.
-negative
   Default on: opposite of -nonegative.
-stop-negative
   Stop on negative components. Not the default.
-fit-spectral-pol <nterms>
   Fit a polynomial over frequency to each clean component. This has only effect
   when the channels are joined with -join-channels.
-fit-spectral-log-pol <nterms>
   Like fit-spectral-pol, but fits a logarithmic polynomial over frequency instead.
-force-spectrum <fitsfile>
   Uses the fits file to force spectral indices (or other/more terms)   during the deconvolution.
-deconvolution-channels <nchannels>
   Decrease the number of channels as specified by -channels-out to the given number for
   deconvolution. Only possible in combination with one of the -fit-spectral options.
   Proper residuals/restored images will only be returned when mgain < 1.
-squared-channel-joining
   Use with -join-channels to perform peak finding in the sum of squared values over
   channels, instead of the normal sum. This is useful for imaging QU polarizations
   with non-zero rotation measures, for which the normal sum is insensitive.
-parallel-deconvolution <maxsize>
   Deconvolve subimages in parallel. Subimages will be at most of the given size.
-deconvolution-threads <n>
   Number of threads to use during deconvolution. On machines with a large nr of cores, this may be used to decrease the memory usage.
   If not specified, the number of threads during deconvolution is controlled with the -j option.

  ** RESTORATION OPTIONS **
-restore <input residual> <input model> <output image>
   Restore the model image onto the residual image and save it in output image. By
   default, the beam parameters are read from the residual image. If this parameter
   is given, wsclean will do the restoring and then exit: no cleaning is performed.
-restore-list <input residual> <input list> <output image>
   Restore a source list onto the residual image and save it in output image. Except
   for the model input format, this parameter behaves equal to -restore.
-beam-size <arcsec>
   Set a circular beam size (FWHM) in arcsec for restoring the clean components. This is
   the same as -beam-shape <size> <size> 0.
-beam-shape <maj in arcsec> <min in arcsec> <position angle in deg>
   Set the FWHM beam shape for restoring the clean components. Defaults units for maj and min are arcsec, and
   degrees for PA. Can be overriden, e.g. '-beam-shape 1amin 1amin 3deg'. Default: shape of PSF.
-fit-beam
   Determine beam shape by fitting the PSF (default if PSF is made).
-no-fit-beam
   Do not determine beam shape from the PSF.
-beam-fitting-size <factor>
   Use a fitting box the size of <factor> times the theoretical beam size for fitting a Gaussian to the PSF.
-theoretic-beam
   Write the beam in output fits files as calculated from the longest projected baseline.
   This method results in slightly less accurate beam size/integrated fluxes, but provides a beam size
   without making the PSF for quick imaging. Default: off.
-circular-beam
   Force the beam to be circular: bmin will be set to bmaj.
-elliptical-beam
   Allow the beam to be elliptical. Default.

For detailed help, check the WSClean website: https://wsclean.readthedocs.io/ .
)";
}

std::vector<std::string> ParseStringList(const char* param) {
  std::vector<std::string> list;
  boost::split(list, param, [](char c) { return c == ','; });
  return list;
}

size_t ParseSizeT(const char* param, const char* name) {
  char* endptr;
  errno = 0;
  const long v = strtol(param, &endptr, 0);
  if (*endptr != 0 || endptr == param || errno != 0) {
    std::ostringstream msg;
    msg << "Could not parse value '" << param << "' for parameter -" << name
        << " to an integer";
    throw std::runtime_error(msg.str());
  }
  if (v < 0) {
    std::ostringstream msg;
    msg << "Invalid value (" << v << ") for parameter -" << name;
    throw std::runtime_error(msg.str());
  }
  return v;
}

double ParseDouble(const char* param, const char* name) {
  char* endptr;
  const double v = std::strtod(param, &endptr);
  if (*endptr != 0 || endptr == param || !std::isfinite(v)) {
    std::ostringstream msg;
    msg << "Could not parse value '" << param << "' for parameter -" << name
        << " to a (double-precision) floating point value";
    throw std::runtime_error(msg.str());
  }
  return v;
}

double ParseDouble(const char* param, double lowerLimit, const char* name,
                   bool inclusive = true) {
  const double v = ParseDouble(param, name);
  if (v < lowerLimit || (v <= lowerLimit && !inclusive)) {
    std::ostringstream msg;
    msg << "Parameter value for -" << name << " was " << v << " but ";
    if (inclusive)
      msg << "is not allowed to be smaller than " << lowerLimit;
    else
      msg << "has to be larger than " << lowerLimit;
    throw std::runtime_error(msg.str());
  }
  return v;
}

}  // namespace

bool CommandLine::ParseWithoutValidation(WSClean& wsclean, int argc,
                                         const char* argv[], bool isSlave) {
  if (argc < 2) {
    if (!isSlave) {
      PrintHeader();
      PrintHelp();
    }
    return false;
  }

  Settings& settings = wsclean.GetSettings();
  int argi = 1;
  bool mfWeighting = false, noMFWeighting = false, dryRun = false;
  std::optional<double> atermKernelSize;
  Logger::SetVerbosity(Logger::kNormalVerbosity);
  while (argi < argc && argv[argi][0] == '-') {
    const std::string param =
        argv[argi][1] == '-' ? (&argv[argi][2]) : (&argv[argi][1]);
    if (param == "version") {
      if (!isSlave) {
        PrintHeader();
#ifdef HAVE_EVERYBEAM
        Logger::Info << "EveryBeam is available.\n";
#endif
#ifdef HAVE_IDG
        Logger::Info << "IDG is available.\n";
#endif
        Logger::Info << "WGridder is available.\n";
      }
      return false;
    } else if (param == "help") {
      if (!isSlave) {
        PrintHeader();
        PrintHelp();
      }
      return false;
    } else if (param == "quiet") {
      Logger::SetVerbosity(Logger::kQuietVerbosity);
    } else if (param == "v" || param == "verbose") {
      Logger::SetVerbosity(Logger::kVerboseVerbosity);
    } else if (param == "log-time") {
      Logger::SetLogTime(true);
    } else if (param == "temp-dir") {
      IncArgi(argi, argc);
      settings.temporaryDirectory = argv[argi];
    } else if (param == "save-weights") {
      settings.isWeightImageSaved = true;
    } else if (param == "save-uv") {
      settings.isUVImageSaved = true;
    } else if (param == "reuse-psf") {
      IncArgi(argi, argc);
      settings.reusePsf = true;
      settings.reusePsfPrefix = argv[argi];
    } else if (param == "reuse-dirty") {
      IncArgi(argi, argc);
      settings.reuseDirty = true;
      settings.reuseDirtyPrefix = argv[argi];
    } else if (param == "predict") {
      settings.mode = Settings::PredictMode;
    } else if (param == "continue") {
      settings.continuedRun = true;
      // Always make a PSF -- otherwise no beam size is available for
      // restoring the existing model.
      settings.makePSF = true;
    } else if (param == "subtract-model") {
      settings.subtractModel = true;
    } else if (param == "size") {
      size_t width = ParseSizeT(argv[argi + 1], "size"),
             height = ParseSizeT(argv[argi + 2], "size");
      settings.trimmedImageWidth = width;
      settings.trimmedImageHeight = height;
      argi += 2;
    } else if (param == "padding") {
      IncArgi(argi, argc);
      settings.imagePadding = ParseDouble(argv[argi], 1.0, "padding");
    } else if (param == "scale") {
      IncArgi(argi, argc);
      settings.pixelScaleX =
          Angle::Parse(argv[argi], "scale parameter", Angle::kDegrees);
      settings.pixelScaleY = settings.pixelScaleX;
    } else if (param == "nwlayers") {
      IncArgi(argi, argc);
      settings.nWLayers = ParseSizeT(argv[argi], "nwlayers");
    } else if (param == "nwlayers-factor") {
      IncArgi(argi, argc);
      settings.nWLayersFactor =
          ParseDouble(argv[argi], 0.0, "nwlayers-factor", false);
    } else if (param == "nwlayers-for-size") {
      settings.widthForNWCalculation =
          ParseSizeT(argv[argi + 1], "nwlayers-for-size");
      settings.heightForNWCalculation =
          ParseSizeT(argv[argi + 2], "nwlayers-for-size");
      argi += 2;
    } else if (param == "gain") {
      IncArgi(argi, argc);
      settings.deconvolutionGain = ParseDouble(argv[argi], 0.0, "gain", false);
    } else if (param == "mgain") {
      IncArgi(argi, argc);
      settings.deconvolutionMGain = ParseDouble(argv[argi], 0.0, "mgain");
    } else if (param == "niter") {
      IncArgi(argi, argc);
      settings.deconvolutionIterationCount = ParseSizeT(argv[argi], "niter");
    } else if (param == "nmiter") {
      IncArgi(argi, argc);
      settings.majorIterationCount = ParseSizeT(argv[argi], "nmiter");
    } else if (param == "threshold") {
      IncArgi(argi, argc);
      settings.deconvolutionThreshold = FluxDensity::Parse(
          argv[argi], "threshold parameter", FluxDensity::kJansky);
    } else if (param == "auto-threshold") {
      IncArgi(argi, argc);
      settings.autoDeconvolutionThreshold = true;
      settings.autoDeconvolutionThresholdSigma =
          ParseDouble(argv[argi], 0.0, "auto-threshold");
    } else if (param == "auto-mask") {
      IncArgi(argi, argc);
      settings.autoMask = true;
      settings.autoMaskSigma = ParseDouble(argv[argi], 0.0, "auto-mask");
    } else if (param == "local-rms") {
      settings.localRMSMethod = radler::LocalRmsMethod::kRmsWindow;
    } else if (param == "local-rms-window") {
      IncArgi(argi, argc);
      settings.localRMSMethod = radler::LocalRmsMethod::kRmsWindow;
      settings.localRMSWindow =
          ParseDouble(argv[argi], 0.0, "local-rms-window", false);
    } else if (param == "local-rms-image") {
      IncArgi(argi, argc);
      settings.localRMSMethod = radler::LocalRmsMethod::kRmsWindow;
      settings.localRMSImage = argv[argi];
    } else if (param == "local-rms-method") {
      IncArgi(argi, argc);
      std::string method = argv[argi];
      if (method == "rms")
        settings.localRMSMethod = radler::LocalRmsMethod::kRmsWindow;
      else if (method == "rms-with-min")
        settings.localRMSMethod = radler::LocalRmsMethod::kRmsAndMinimumWindow;
      else
        throw std::runtime_error("Unknown RMS background method specified");
    } else if (param == "data-column") {
      IncArgi(argi, argc);
      settings.dataColumnName = argv[argi];
    } else if (param == "pol") {
      IncArgi(argi, argc);
      settings.polarizations = aocommon::Polarization::ParseList(argv[argi]);
    } else if (param == "beam-model") {
      IncArgi(argi, argc);
      std::string beamModel = argv[argi];
      boost::to_upper(beamModel);
      if (beamModel == "HAMAKER" || beamModel == "LOBES" ||
          beamModel == "OSKARDIPOLE" || beamModel == "OSKARSPHERICALWAVE") {
        settings.beamModel = beamModel;
      } else {
        throw std::runtime_error(
            "Invalid beam-model: should be either Hamaker, Lobes, OskarDipole "
            "or OskarSphericalWave (case insensitive)");
      }
    } else if (param == "beam-mode") {
      IncArgi(argi, argc);
      std::string beamMode = argv[argi];
      boost::to_upper(beamMode);
      if (beamMode == "ARRAY_FACTOR" || beamMode == "ELEMENT" ||
          beamMode == "FULL") {
        settings.beamMode = beamMode;
      } else {
        throw std::runtime_error(
            "Invalid beam-mode: should be either array_factor, element or full "
            "(case insensitive)");
      }
    } else if (param == "beam-normalisation-mode") {
      IncArgi(argi, argc);
      std::string beamNormalisationMode = argv[argi];
      boost::to_upper(beamNormalisationMode);
      if (beamNormalisationMode == "NONE" ||
          beamNormalisationMode == "PREAPPLIED" ||
          beamNormalisationMode == "AMPLITUDE" ||
          beamNormalisationMode == "FULL") {
        settings.beamNormalisationMode = beamNormalisationMode;
      } else {
        throw std::runtime_error(
            "Invalid beam-normalisation-mode: should be either none, "
            "preapplied, amplitude or full "
            "(case insensitive)");
      }
    } else if (param == "apply-primary-beam") {
      settings.applyPrimaryBeam = true;
    } else if (param == "reuse-primary-beam") {
      settings.reusePrimaryBeam = true;
    } else if (param == "use-differential-lofar-beam") {
      // pre_applied_or_full is the beam normalisation mode
      // that implements the behaviour of
      // the old use_differential_beam option of EveryBeam
      settings.beamNormalisationMode = "preapplied_or_full";
    } else if (param == "primary-beam-limit") {
      IncArgi(argi, argc);
      settings.primaryBeamLimit =
          ParseDouble(argv[argi], 0.0, "primary-beam-limit");
    } else if (param == "scalar-beam") {
      settings.useScalarPrimaryBeam = true;
    } else if (param == "mwa-path") {
      IncArgi(argi, argc);
      settings.mwaPath = argv[argi];
    } else if (param == "dry-run") {
      dryRun = true;
    } else if (param == "save-psf-pb") {
      settings.savePsfPb = true;
    } else if (param == "pb-grid-size") {
      IncArgi(argi, argc);
      settings.primaryBeamGridSize = ParseSizeT(argv[argi], "pb-grid-size");
    } else if (param == "dd-psf-grid") {
      settings.ddPsfGridWidth = ParseSizeT(argv[argi + 1], "dd-psf-grid");
      settings.ddPsfGridHeight = ParseSizeT(argv[argi + 2], "dd-psf-grid");
      argi += 2;
    } else if (param == "negative") {
      settings.allowNegativeComponents = true;
    } else if (param == "no-negative") {
      settings.allowNegativeComponents = false;
    } else if (param == "stop-negative") {
      settings.stopOnNegativeComponents = true;
    } else if (param == "python-deconvolution") {
      IncArgi(argi, argc);
      settings.algorithmType = radler::AlgorithmType::kPython;
      settings.pythonDeconvolutionFilename = argv[argi];
      settings.deconvolutionIterationCount =
          std::max(size_t{1}, settings.deconvolutionIterationCount);
    } else if (param == "iuwt") {
      settings.algorithmType = radler::AlgorithmType::kIuwt;
      // Currently (WSClean 1.9, 2015-08-19) IUWT deconvolution
      // seems not to work when allowing negative components. The algorithm
      // becomes unstable. Hence, turn negative components off.
      settings.allowNegativeComponents = false;
    } else if (param == "iuwt-snr-test") {
      settings.iuwtSNRTest = true;
    } else if (param == "no-iuwt-snr-test") {
      settings.iuwtSNRTest = false;
    } else if (param == "moresane-ext") {
      IncArgi(argi, argc);
      settings.algorithmType = radler::AlgorithmType::kMoreSane;
      settings.moreSaneLocation = argv[argi];
    } else if (param == "moresane-arg") {
      IncArgi(argi, argc);
      settings.moreSaneArgs = argv[argi];
    } else if (param == "moresane-sl") {
      IncArgi(argi, argc);
      settings.moreSaneSigmaLevels = NumberList::ParseDoubleList(argv[argi]);
    } else if (param == "make-psf") {
      settings.makePSF = true;
    } else if (param == "make-psf-only") {
      settings.makePSFOnly = true;
    } else if (param == "name") {
      IncArgi(argi, argc);
      settings.prefixName = argv[argi];
    } else if (param == "grid-mode") {
      IncArgi(argi, argc);
      std::string gridModeStr = argv[argi];
      boost::to_lower(gridModeStr);
      if (gridModeStr == "kb" || gridModeStr == "kaiserbessel" ||
          gridModeStr == "kaiser-bessel")
        settings.gridMode = GriddingKernelMode::KaiserBessel;
      else if (gridModeStr == "bn")
        settings.gridMode = GriddingKernelMode::BlackmanNuttall;
      else if (gridModeStr == "bh")
        settings.gridMode = GriddingKernelMode::BlackmanHarris;
      else if (gridModeStr == "gaus")
        settings.gridMode = GriddingKernelMode::Gaussian;
      else if (gridModeStr == "rect")
        settings.gridMode = GriddingKernelMode::Rectangular;
      else if (gridModeStr == "kb-no-sinc")
        settings.gridMode = GriddingKernelMode::KaiserBesselWithoutSinc;
      else if (gridModeStr == "nn" || gridModeStr == "nearestneighbour")
        settings.gridMode = GriddingKernelMode::NearestNeighbour;
      else
        throw std::runtime_error(
            "Invalid gridding mode: should be either kb (Kaiser-Bessel), nn "
            "(NearestNeighbour), bn, bh, gaus, kb-no-sinc or rect");
    } else if (param == "small-inversion") {
      settings.smallInversion = true;
    } else if (param == "no-small-inversion") {
      settings.smallInversion = false;
    } else if (param == "interval") {
      settings.startTimestep = ParseSizeT(argv[argi + 1], "interval");
      settings.endTimestep = ParseSizeT(argv[argi + 2], "interval");
      argi += 2;
    } else if (param == "intervals-out") {
      IncArgi(argi, argc);
      settings.intervalsOut = atoi(argv[argi]);
    } else if (param == "even-timesteps") {
      settings.evenOddTimesteps = MSSelection::EvenTimesteps;
    } else if (param == "odd-timesteps") {
      settings.evenOddTimesteps = MSSelection::OddTimesteps;
    } else if (param == "channel-range") {
      settings.startChannel = ParseSizeT(argv[argi + 1], "channel-range");
      settings.endChannel = ParseSizeT(argv[argi + 2], "channel-range");
      argi += 2;
    } else if (param == "shift") {
      settings.hasShift = true;
      settings.shiftRA = aocommon::RaDecCoord::ParseRA(argv[argi + 1]);
      settings.shiftDec = aocommon::RaDecCoord::ParseDec(argv[argi + 2]);
      argi += 2;
    } else if (param == "channels-out") {
      IncArgi(argi, argc);
      settings.channelsOut = ParseSizeT(argv[argi], "channels-out");
    } else if (param == "gap-channel-division") {
      settings.divideChannelsByGaps = true;
    } else if (param == "channel-division-frequencies") {
      IncArgi(argi, argc);
      settings.divideChannelFrequencies =
          NumberList::ParseDoubleList(argv[argi]);
    } else if (param == "facet-regions") {
      IncArgi(argi, argc);
      settings.facetRegionFilename = argv[argi];
    } else if (param == "join-polarizations") {
      settings.joinedPolarizationDeconvolution = true;
    } else if (param == "link-polarizations") {
      IncArgi(argi, argc);
      settings.joinedPolarizationDeconvolution = true;
      settings.linkedPolarizations =
          aocommon::Polarization::ParseList(argv[argi]);
    } else if (param == "join-channels") {
      settings.joinedFrequencyDeconvolution = true;
    } else if (param == "mf-weighting" || param == "mfs-weighting") {
      mfWeighting = true;
      // mfs was renamed to mf in wsclean 2.7
      if (param != "mf-weighting") Deprecated(isSlave, param, "mf-weighting");
    } else if (param == "no-mf-weighting" || param == "no-mfs-weighting") {
      noMFWeighting = true;
      // mfs was renamed to mf in wsclean 2.7
      if (param != "no-mf-weighting")
        Deprecated(isSlave, param, "no-mf-weighting");
    } else if (param == "spectral-correction") {
      settings.spectralCorrectionFrequency =
          ParseDouble(argv[argi + 1], 0.0, "spectral-correction", false);
      aocommon::UVector<double> list =
          NumberList::ParseDoubleList(argv[argi + 2]);
      settings.spectralCorrection.assign(list.begin(), list.end());
      argi += 2;
    } else if (param == "taper-gaussian") {
      IncArgi(argi, argc);
      double taperBeamSize =
          Angle::Parse(argv[argi], "Gaussian taper", Angle::kArcseconds);
      settings.gaussianTaperBeamSize = taperBeamSize;
    } else if (param == "taper-edge") {
      IncArgi(argi, argc);
      settings.edgeTaperInLambda = ParseDouble(argv[argi], 0.0, "taper-edge");
    } else if (param == "taper-edge-tukey") {
      IncArgi(argi, argc);
      settings.edgeTukeyTaperInLambda =
          ParseDouble(argv[argi], 0.0, "taper-edge-tukey");
    } else if (param == "taper-tukey") {
      IncArgi(argi, argc);
      settings.tukeyTaperInLambda = ParseDouble(argv[argi], 0.0, "taper-tukey");
    } else if (param == "taper-inner-tukey") {
      IncArgi(argi, argc);
      settings.tukeyInnerTaperInLambda =
          ParseDouble(argv[argi], 0.0, "taper-inner-tukey");
    } else if (param == "use-weights-as-taper") {
      settings.useWeightsAsTaper = true;
    } else if (param == "store-imaging-weights") {
      settings.writeImagingWeightSpectrumColumn = true;
    } else if (param == "no-fast-subminor") {
      settings.useSubMinorOptimization = false;
    } else if (param == "multiscale") {
      settings.algorithmType = radler::AlgorithmType::kMultiscale;
    } else if (param == "multiscale-gain") {
      IncArgi(argi, argc);
      settings.multiscaleGain =
          ParseDouble(argv[argi], 0.0, "multiscale-gain", false);
    } else if (param == "multiscale-scale-bias") {
      IncArgi(argi, argc);
      settings.multiscaleDeconvolutionScaleBias =
          ParseDouble(argv[argi], 0.0, "multiscale-scale-bias", false);
    } else if (param == "multiscale-max-scales") {
      IncArgi(argi, argc);
      settings.multiscaleMaxScales =
          ParseSizeT(argv[argi], "multiscale-max-scales");
    } else if (param == "multiscale-scales") {
      IncArgi(argi, argc);
      settings.multiscaleScaleList = NumberList::ParseDoubleList(argv[argi]);
    } else if (param == "multiscale-shape") {
      IncArgi(argi, argc);
      std::string shape = argv[argi];
      if (shape == "tapered-quadratic")
        settings.multiscaleShapeFunction =
            radler::MultiscaleShape::kTaperedQuadraticShape;
      else if (shape == "gaussian")
        settings.multiscaleShapeFunction =
            radler::MultiscaleShape::kGaussianShape;
      else
        throw std::runtime_error("Unknown multiscale shape function given");
    } else if (param == "multiscale-convolution-padding") {
      IncArgi(argi, argc);
      settings.multiscaleConvolutionPadding =
          ParseDouble(argv[argi], 1.0, "multiscale-convolution-padding");
    } else if (param == "no-multiscale-fast-subminor") {
      settings.multiscaleFastSubMinorLoop = false;
    } else if (param == "weighting-rank-filter") {
      IncArgi(argi, argc);
      settings.rankFilterLevel =
          ParseDouble(argv[argi], 0.0, "weighting-rank-filter");
    } else if (param == "weighting-rank-filter-size") {
      IncArgi(argi, argc);
      settings.rankFilterSize =
          ParseSizeT(argv[argi], "weighting-rank-filter-size");
    } else if (param == "save-source-list") {
      settings.saveSourceList = true;
      settings.multiscaleShapeFunction =
          radler::MultiscaleShape::kGaussianShape;
    } else if (param == "clean-border") {
      IncArgi(argi, argc);
      settings.deconvolutionBorderRatio =
          ParseDouble(argv[argi], 0.0, "clean-border") * 0.01;
    } else if (param == "fits-mask") {
      IncArgi(argi, argc);
      settings.fitsDeconvolutionMask = argv[argi];
    } else if (param == "casa-mask") {
      IncArgi(argi, argc);
      settings.casaDeconvolutionMask = argv[argi];
    } else if (param == "horizon-mask") {
      IncArgi(argi, argc);
      settings.horizonMask = true;
      settings.horizonMaskDistance =
          Angle::Parse(argv[argi], "horizon mask distance", Angle::kDegrees);
    } else if (param == "fit-spectral-pol") {
      IncArgi(argi, argc);
      settings.spectralFittingMode =
          schaapcommon::fitters::SpectralFittingMode::kPolynomial;
      settings.spectralFittingTerms =
          ParseSizeT(argv[argi], "fit-spectral-pol");
    } else if (param == "fit-spectral-log-pol") {
      IncArgi(argi, argc);
      settings.spectralFittingMode =
          schaapcommon::fitters::SpectralFittingMode::kLogPolynomial;
      settings.spectralFittingTerms =
          ParseSizeT(argv[argi], "fit-spectral-log-pol");
    } else if (param == "force-spectrum") {
      IncArgi(argi, argc);
      settings.forcedSpectrumFilename = argv[argi];
    } else if (param == "deconvolution-channels") {
      IncArgi(argi, argc);
      settings.deconvolutionChannelCount =
          ParseSizeT(argv[argi], "deconvolution-channels");
    } else if (param == "squared-channel-joining") {
      settings.squaredJoins = true;
    } else if (param == "parallel-deconvolution") {
      IncArgi(argi, argc);
      settings.parallelDeconvolutionMaxSize =
          ParseSizeT(argv[argi], "parallel-deconvolution");
    } else if (param == "deconvolution-threads") {
      IncArgi(argi, argc);
      settings.parallelDeconvolutionMaxThreads =
          ParseSizeT(argv[argi], "deconvolution-threads");
    } else if (param == "field") {
      IncArgi(argi, argc);
      if (argv[argi] == std::string("all"))
        settings.fieldIds.assign(1, MSSelection::ALL_FIELDS);
      else {
        aocommon::UVector<int> list = NumberList::ParseIntList(argv[argi]);
        settings.fieldIds.assign(list.begin(), list.end());
      }
    } else if (param == "spws") {
      IncArgi(argi, argc);
      aocommon::UVector<int> list = NumberList::ParseIntList(argv[argi]);
      settings.spectralWindows.insert(list.begin(), list.end());
    } else if (param == "weight") {
      IncArgi(argi, argc);
      std::string weightArg = argv[argi];
      if (weightArg == "natural")
        settings.weightMode = WeightMode(WeightMode::NaturalWeighted);
      else if (weightArg == "uniform")
        settings.weightMode = WeightMode(WeightMode::UniformWeighted);
      else if (weightArg == "briggs") {
        IncArgi(argi, argc);
        settings.weightMode =
            WeightMode::Briggs(ParseDouble(argv[argi], "weight briggs"));
      } else
        throw std::runtime_error("Unknown weighting mode specified");
    } else if (param == "super-weight") {
      IncArgi(argi, argc);
      settings.weightMode.SetSuperWeight(
          ParseDouble(argv[argi], 0.0, "super-weight"));
    } else if (param == "restore" || param == "restore-list") {
      if (param == "restore")
        settings.mode = Settings::RestoreMode;
      else
        settings.mode = Settings::RestoreListMode;
      settings.restoreInput = argv[argi + 1];
      settings.restoreModel = argv[argi + 2];
      settings.restoreOutput = argv[argi + 3];
      argi += 3;
    } else if (param == "beam-size") {
      IncArgi(argi, argc);
      double beam = Angle::Parse(argv[argi], "beam size", Angle::kArcseconds);
      settings.manualBeamMajorSize = beam;
      settings.manualBeamMinorSize = beam;
      settings.manualBeamPA = 0.0;
    } else if (param == "beam-shape") {
      double beamMaj = Angle::Parse(argv[argi + 1], "beam shape, major axis",
                                    Angle::kArcseconds);
      double beamMin = Angle::Parse(argv[argi + 2], "beam shape, minor axis",
                                    Angle::kArcseconds);
      double beamPA = Angle::Parse(argv[argi + 3], "beam shape, position angle",
                                   Angle::kDegrees);
      argi += 3;
      settings.manualBeamMajorSize = beamMaj;
      settings.manualBeamMinorSize = beamMin;
      settings.manualBeamPA = beamPA;
    } else if (param == "fit-beam") {
      settings.fittedBeam = true;
    } else if (param == "no-fit-beam") {
      settings.fittedBeam = false;
    } else if (param == "beam-fitting-size") {
      IncArgi(argi, argc);
      settings.beamFittingBoxSize =
          ParseDouble(argv[argi], 0.0, "beam-fitting-size", false);
    } else if (param == "theoretic-beam") {
      settings.theoreticBeam = true;
      settings.fittedBeam = false;
    } else if (param == "circular-beam") {
      settings.circularBeam = true;
    } else if (param == "elliptical-beam") {
      settings.circularBeam = false;
    } else if (param == "kernel-size") {
      IncArgi(argi, argc);
      settings.antialiasingKernelSize = ParseSizeT(argv[argi], "kernel-size");
    } else if (param == "oversampling") {
      IncArgi(argi, argc);
      settings.overSamplingFactor = ParseSizeT(argv[argi], "oversampling");
    } else if (param == "reorder") {
      settings.forceReorder = true;
      settings.forceNoReorder = false;
    } else if (param == "no-reorder") {
      settings.forceNoReorder = true;
      settings.forceReorder = false;
    } else if (param == "update-model-required") {
      settings.modelUpdateRequired = true;
    } else if (param == "no-update-model-required") {
      settings.modelUpdateRequired = false;
    } else if (param == "j") {
      IncArgi(argi, argc);
      settings.threadCount = ParseSizeT(argv[argi], "j");
    } else if (param == "parallel-reordering") {
      IncArgi(argi, argc);
      settings.parallelReordering =
          ParseSizeT(argv[argi], "parallel-reordering");
    } else if (param == "parallel-gridding") {
      IncArgi(argi, argc);
      settings.parallelGridding = ParseSizeT(argv[argi], "parallel-gridding");
    } else if (param == "no-work-on-master") {
      settings.masterDoesWork = false;
    } else if (param == "mem") {
      IncArgi(argi, argc);
      settings.memFraction = ParseDouble(argv[argi], 0.0, "mem", false) / 100.0;
    } else if (param == "abs-mem") {
      IncArgi(argi, argc);
      settings.absMemLimit = ParseDouble(argv[argi], 0.0, "abs-mem", false);
    } else if (param == "maxuvw-m") {
      IncArgi(argi, argc);
      settings.maxUVWInMeters = ParseDouble(argv[argi], 0.0, "maxuvw-m", false);
    } else if (param == "minuvw-m") {
      IncArgi(argi, argc);
      settings.minUVWInMeters = ParseDouble(argv[argi], 0.0, "minuvw-m");
    } else if (param == "maxuv-l") {
      IncArgi(argi, argc);
      settings.maxUVInLambda = ParseDouble(argv[argi], 0.0, "maxuv-l", false);
    } else if (param == "minuv-l") {
      IncArgi(argi, argc);
      settings.minUVInLambda = ParseDouble(argv[argi], 0.0, "minuv-l");
    } else if (param == "maxw") {
      // This was to test the optimization suggested in Tasse et al., 2013,
      // Appendix C.
      IncArgi(argi, argc);
      settings.wLimit = ParseDouble(argv[argi], 0.0, "maxw");
    } else if (param == "baseline-averaging") {
      IncArgi(argi, argc);
      settings.baselineDependentAveragingInWavelengths =
          ParseDouble(argv[argi], 0.0, "baseline-averaging", false);
    } else if (param == "simulate-noise") {
      IncArgi(argi, argc);
      settings.simulateNoise = true;
      settings.simulatedNoiseStdDev =
          ParseDouble(argv[argi], 0.0, "simulate-noise");
    } else if (param == "simulate-baseline-noise") {
      IncArgi(argi, argc);
      settings.simulateNoise = true;
      settings.simulatedBaselineNoiseFilename = argv[argi];
    } else if (param == "aterm-config") {
      IncArgi(argi, argc);
      settings.atermConfigFilename = argv[argi];
    } else if (param == "grid-with-beam") {
      settings.gridWithBeam = true;
    } else if (param == "beam-aterm-update") {
      IncArgi(argi, argc);
      double val = ParseDouble(argv[argi], 0.0, "beam-aterm-update");
      settings.beamAtermUpdateTime = val;
      settings.primaryBeamUpdateTime = std::max<size_t>(val, 1.0);
    } else if (param == "aterm-kernel-size") {
      IncArgi(argi, argc);
      atermKernelSize = ParseDouble(argv[argi], 0.0, "aterm-kernel-size");
    } else if (param == "apply-facet-solutions") {
      IncArgi(argi, argc);
      settings.facetSolutionFiles = ParseStringList(argv[argi]);
      IncArgi(argi, argc);
      settings.facetSolutionTables = ParseStringList(argv[argi]);
      if (settings.facetSolutionTables.size() > 2) {
        throw std::runtime_error(
            "List of solution tables (soltabs) should contain at most two "
            "entries.");
      }
    } else if (param == "apply-facet-beam") {
      settings.applyFacetBeam = true;
    } else if (param == "facet-beam-update") {
      IncArgi(argi, argc);
      settings.facetBeamUpdateTime =
          ParseDouble(argv[argi], 0.0, "facet-beam-update");
    } else if (param == "save-aterms") {
      settings.saveATerms = true;
    } else if (param == "visibility-weighting-mode") {
      IncArgi(argi, argc);
      std::string modeStr = argv[argi];
      boost::to_lower(modeStr);
      if (modeStr == "normal")
        settings.visibilityWeightingMode =
            VisibilityWeightingMode::NormalVisibilityWeighting;
      else if (modeStr == "squared")
        settings.visibilityWeightingMode =
            VisibilityWeightingMode::SquaredVisibilityWeighting;
      else if (modeStr == "unit")
        settings.visibilityWeightingMode =
            VisibilityWeightingMode::UnitVisibilityWeighting;
      else
        throw std::runtime_error("Unknown weighting mode: " + modeStr);
    } else if (param == "direct-ft") {
      Deprecated(isSlave, param, "gridder");
      settings.gridderType = GridderType::DirectFT;
      settings.imagePadding = 1.0;
      settings.smallInversion = false;
    } else if (param == "direct-ft-precision") {
      IncArgi(argi, argc);
      std::string precStr = argv[argi];
      if (precStr == "float")
        settings.directFTPrecision = DirectFTPrecision::Float;
      else if (precStr == "double")
        settings.directFTPrecision = DirectFTPrecision::Double;
      else if (precStr == "ldouble")
        settings.directFTPrecision = DirectFTPrecision::LongDouble;
      else
        throw std::runtime_error(
            "Invalid direct ft precision specified. Allowed options: float, "
            "double and ldouble.");
    } else if (param == "gridder") {
      IncArgi(argi, argc);
      const std::string gridder_str = argv[argi];
      if (gridder_str == "idg") {
#if !defined(HAVE_IDG)
        throw std::runtime_error(
            "WSClean was not compiled with IDG: to use it, install IDG and "
            "recompile WSClean");
#endif
        settings.gridderType = GridderType::IDG;
        settings.smallInversion = false;
      } else if (gridder_str == "wgridder") {
        settings.gridderType = GridderType::WGridder;
      } else if (gridder_str == "tuned-wgridder") {
        settings.gridderType = GridderType::TunedWGridder;
      } else if (gridder_str == "wstacking") {
        settings.gridderType = GridderType::WStacking;
      } else if (gridder_str == "direct-ft") {
        settings.gridderType = GridderType::DirectFT;
      } else
        throw std::runtime_error("Invalid gridder requested: '" + gridder_str +
                                 "'");
    } else if (param == "use-idg") {
      Deprecated(isSlave, param, "gridder");
#if !defined(HAVE_IDG)
      throw std::runtime_error(
          "WSClean was not compiled with IDG: to use it, install IDG and "
          "recompile WSClean");
#endif
      settings.gridderType = GridderType::IDG;
      settings.smallInversion = false;
    } else if (param == "idg-mode") {
      IncArgi(argi, argc);
      std::string mode =
          boost::algorithm::to_lower_copy(std::string(argv[argi]));
      if (mode == "cpu")
        settings.idgMode = Settings::IDG_CPU;
      else if (mode == "gpu")
        settings.idgMode = Settings::IDG_GPU;
      else if (mode == "hybrid")
        settings.idgMode = Settings::IDG_HYBRID;
      else
        throw std::runtime_error("Unknown IDG mode: " + mode);
    } else if (param == "use-wgridder") {
      Deprecated(isSlave, param, "gridder");
      settings.gridderType = GridderType::WGridder;
    } else if (param == "wgridder-accuracy") {
      IncArgi(argi, argc);
      settings.wgridderAccuracy =
          ParseDouble(argv[argi], 0.0, "wgridder-accuracy", false);
    } else if (param == "no-dirty") {
      settings.isDirtySaved = false;
    } else if (param == "save-first-residual") {
      settings.isFirstResidualSaved = true;
    } else {
      throw std::runtime_error("Unknown parameter: " + param);
    }

    ++argi;
  }

  if (argi == argc && settings.mode != Settings::RestoreMode &&
      settings.mode != Settings::RestoreListMode)
    throw std::runtime_error("No input measurement sets given.");

  // Done parsing.

  // We print the header only now, because the logger has now been set up
  // and possibly set to quiet.
  if (!isSlave) PrintHeader();

  const size_t defaultAtermSize = settings.atermConfigFilename.empty() ? 5 : 16;
  settings.atermKernelSize = atermKernelSize.value_or(defaultAtermSize);

  settings.mfWeighting =
      (settings.joinedFrequencyDeconvolution && !noMFWeighting) || mfWeighting;

  // Joined polarizations is implemented by linking all polarizations
  if (settings.joinedPolarizationDeconvolution &&
      settings.linkedPolarizations.empty()) {
    settings.linkedPolarizations = settings.polarizations;
  }

  for (int i = argi; i != argc; ++i) settings.filenames.push_back(argv[i]);

  std::ostringstream commandLineStr;
  commandLineStr << "wsclean";
  for (int i = 1; i != argc; ++i) commandLineStr << ' ' << argv[i];
  wsclean.SetCommandLine(commandLineStr.str());

  return !dryRun;
}

void CommandLine::Validate(WSClean& wsclean) {
  wsclean.GetSettings().Validate();
  wsclean.GetSettings().Propagate();
}

void CommandLine::Run(class WSClean& wsclean) {
  const Settings& settings = wsclean.GetSettings();
  switch (settings.mode) {
    case Settings::RestoreMode:
      WSCFitsWriter::Restore(settings);
      break;
    case Settings::RestoreListMode:
      WSCFitsWriter::RestoreList(settings);
      break;
    case Settings::PredictMode:
      wsclean.RunPredict();
      break;
    case Settings::ImagingMode:
      wsclean.RunClean();
      break;
  }
}
