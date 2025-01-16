#ifndef WSCLEAN_SETTINGS_H
#define WSCLEAN_SETTINGS_H

#include <cassert>
#include <cstdint>
#include <limits>
#include <optional>

#include "../gridding/visibilityweightingmode.h"
#include "../gridding/wstackinggridder.h"

#include "../structures/msselection.h"
#include "../structures/weightmode.h"

#include <aocommon/system.h>

#include <schaapcommon/fitters/spectralfitter.h>
#include <schaapcommon/reordering/storagemanagertype.h>

#include <radler/settings.h>

namespace wsclean {

enum class DirectFTPrecision { Float, Double, LongDouble };
enum class GridderType { WStacking, WGridder, TunedWGridder, DirectFT, IDG };
enum class VisibilityReadMode { kScalar, kDiagonal, kFull };

/**
 * This class describes all settings for a single WSClean run.
 * @sa WSClean
 */
class Settings {
 public:
  Settings() = default;

  void Validate() const;

  void Propagate(bool verbose = true);

  void RecalculateDerivedDimensions(bool verbose = true);

  std::vector<std::string> filenames;
  enum Mode {
    ImagingMode,
    PredictMode,
    RestoreMode,
    RestoreListMode
  } mode = ImagingMode;
  GridderType gridderType = GridderType::WGridder;
  size_t paddedImageWidth = 0, paddedImageHeight = 0;
  size_t trimmedImageWidth = 0, trimmedImageHeight = 0;
  bool hasShift = false;
  double shiftRA = 0.0, shiftDec = 0.0;
  double imagePadding = 1.2;
  size_t widthForNWCalculation = 0, heightForNWCalculation = 0;
  size_t channelsOut = 1, intervalsOut = 1;
  enum schaapcommon::reordering::MSSelection::EvenOddSelection
      evenOddTimesteps = schaapcommon::reordering::MSSelection::kAllTimesteps;
  bool divideChannelsByGaps = false;
  aocommon::UVector<double> divideChannelFrequencies;
  double pixelScaleX = 0.0, pixelScaleY = 0.0;
  std::string restoreModel, restoreInput, restoreOutput;
  double manualBeamMajorSize = 0.0, manualBeamMinorSize = 0.0;
  double manualBeamPA = 0.0;
  bool fittedBeam = true, theoreticBeam = false, circularBeam = false;
  double beamFittingBoxSize = 10.0;
  bool continuedRun = false;
  double memFraction = 1.0, absMemLimit = 0.0;
  double minUVWInMeters = 0.0, maxUVWInMeters = 0.0, minUVInLambda = 0.0;
  double maxUVInLambda = 0.0;
  double wLimit = 0.0;
  double rankFilterLevel = 3.0;
  size_t rankFilterSize = 16;
  double gaussianTaperBeamSize = 0.0, tukeyTaperInLambda = 0.0;
  double tukeyInnerTaperInLambda = 0.0;
  double edgeTaperInLambda = 0.0, edgeTukeyTaperInLambda = 0.0;
  bool useWeightsAsTaper = false;
  size_t nWLayers = 0;
  double nWLayersFactor = 1.0;
  size_t antialiasingKernelSize = 7, overSamplingFactor = 1023;
  size_t threadCount = aocommon::system::ProcessorCount();
  size_t parallelReordering = 4, parallelGridding = 1;
  size_t nMpiNodes = 0;  // 0 if MPI is disabled.
  size_t maxMpiMessageSize =
      std::numeric_limits<std::int32_t>::max();  // Only used if MPI is enabled.
  bool masterDoesWork = true;
  /// Maps output channel indices to node indices, when using MPI.
  std::vector<size_t> channelToNode;
  std::vector<size_t> fieldIds{0};
  size_t startTimestep = 0, endTimestep = 0;
  size_t startChannel = 0, endChannel = 0;
  std::string dataColumnName;
  std::string modelColumnName = "MODEL_DATA";
  schaapcommon::reordering::StorageManagerType modelStorageManager =
      schaapcommon::reordering::StorageManagerType::Default;
  std::set<aocommon::PolarizationEnum> polarizations{
      aocommon::Polarization::StokesI};
  std::string facetRegionFilename;
  std::optional<size_t> featherSize;
  std::set<size_t> spectralWindows;
  WeightMode weightMode{WeightClass::Uniform};
  std::string prefixName = "wsclean";
  bool joinedPolarizationDeconvolution = false;
  bool joinedFrequencyDeconvolution = false;
  bool minGridResolution = true;
  bool makePSF = false;
  bool makePSFOnly = false;
  bool isWeightImageSaved = false;
  bool isUVImageSaved = false;
  bool isDirtySaved = true;
  bool isFirstResidualSaved = false;
  bool reusePsf = false, reuseDirty = false;
  std::string reusePsfPrefix, reuseDirtyPrefix;
  bool writeImagingWeightSpectrumColumn = false;
  std::string temporaryDirectory;
  bool forceReorder = false;
  bool forceNoReorder = false;
  bool doReorder = true;
  bool reuseReorder = false;
  bool saveReorder = false;
  bool subtractModel = false, modelUpdateRequired = true, mfWeighting = false;
  size_t fullResOffset = 0, fullResWidth = 0, fullResPad = 0;
  std::string beamModel, beamMode = "default";
  std::string beamNormalisationMode = "preapplied";
  bool applyPrimaryBeam = false, reusePrimaryBeam = false;
  bool savePsfPb = false;
  double primaryBeamLimit = 0.005;
  std::string mwaPath;
  size_t primaryBeamGridSize = 32, primaryBeamUpdateTime = 1800;
  size_t ddPsfGridHeight = 1, ddPsfGridWidth = 1;
  DirectFTPrecision directFTPrecision = DirectFTPrecision::Double;
  double wgridderAccuracy = 1e-4;
  std::string atermConfigFilename;
  double atermKernelSize = 5.0;
  bool gridWithBeam = false;
  double beamAtermUpdateTime = 300.0;  // in seconds.
  std::vector<std::string> facetSolutionFiles;
  std::vector<std::string> facetSolutionTables;
  bool solutionDirectionsCheck = true;
  VisibilityReadMode visibilityReadMode = VisibilityReadMode::kFull;
  bool applyFacetBeam = false;
  double facetBeamUpdateTime = 120.0;  // in seconds.
  bool saveATerms = false;
  enum IDGMode {
    IDG_DEFAULT,
    IDG_GPU,
    IDG_CPU,
    IDG_HYBRID
  } idgMode = IDG_DEFAULT;
  enum GriddingKernelMode gridMode = GriddingKernelMode::KaiserBessel;
  enum VisibilityWeightingMode visibilityWeightingMode =
      VisibilityWeightingMode::NormalVisibilityWeighting;
  double baselineDependentAveragingInWavelengths = 0.0;
  bool simulateNoise = false;
  double simulatedNoiseStdDev = 0.0;
  std::string simulatedBaselineNoiseFilename;
  bool compound_tasks = false;
  bool shared_facet_reads = false;

  /** @{
   * These settings all relate to the deconvolution.
   */
  std::set<aocommon::PolarizationEnum> linkedPolarizations;
  size_t parallelDeconvolutionMaxSize = 0;
  size_t parallelDeconvolutionGridWidth = 0;   // derived setting
  size_t parallelDeconvolutionGridHeight = 0;  // derived setting
  size_t parallelDeconvolutionMaxThreads = 0;
  std::optional<double> autoDeconvolutionThresholdSigma;
  std::optional<double> absoluteDeconvolutionThreshold;
  std::optional<double> autoMaskSigma;
  std::optional<double> absoluteAutoMaskThreshold;
  double deconvolutionGain = 0.1, deconvolutionMGain = 1.0;
  double localRMSWindow = 25.0;
  double localRMSStrength = 1.0;
  radler::LocalRmsMethod localRMSMethod = radler::LocalRmsMethod::kNone;
  bool saveSourceList = false;
  size_t deconvolutionIterationCount = 0;
  size_t majorIterationCount = 12;
  bool allowNegativeComponents = true, stopOnNegativeComponents = false;
  bool useSubMinorOptimization = true, squaredJoins = false;
  double spectralCorrectionFrequency = 0.0;
  std::vector<float> spectralCorrection;
  bool multiscaleFastSubMinorLoop = true;
  double multiscaleGain = 0.2, multiscaleDeconvolutionScaleBias = 0.6;
  size_t multiscaleMaxScales = 0;
  double multiscaleConvolutionPadding = 1.1;
  aocommon::UVector<double> multiscaleScaleList;
  radler::MultiscaleShape multiscaleShapeFunction =
      radler::MultiscaleShape::kTaperedQuadraticShape;
  double deconvolutionBorderRatio = 0.0;
  std::string fitsDeconvolutionMask, casaDeconvolutionMask;
  bool horizonMask = false;
  double horizonMaskDistance = 0.0;
  std::string localRMSImage;
  std::string pythonDeconvolutionFilename;
  bool iuwtSNRTest = false;
  std::string moreSaneLocation, moreSaneArgs;
  aocommon::UVector<double> moreSaneSigmaLevels;
  schaapcommon::fitters::SpectralFittingMode spectralFittingMode =
      schaapcommon::fitters::SpectralFittingMode::kNoFitting;
  size_t spectralFittingTerms = 0;
  std::string forcedSpectrumFilename;
  /**
   * The number of channels used during deconvolution. This can be used to
   * image with more channels than deconvolution. Before deconvolution,
   * channels are averaged, and after deconvolution they are interpolated.
   * It is 0 when all channels should be used.
   */
  size_t deconvolutionChannelCount = 0;

  /**
   * Type of deconvolution algorithm.
   */
  radler::AlgorithmType algorithmType = radler::AlgorithmType::kGenericClean;
  /**
   * @}
   */

  /**
   * @brief Extract the settings that are relevant to the deconvolution.
   * Currently, it duplicates the existing settings into a DeconvolutionSettings
   * object.
   */
  radler::Settings GetRadlerSettings() const;

  schaapcommon::reordering::MSSelection GetMSSelection() const {
    schaapcommon::reordering::MSSelection selection;
    selection.SetInterval(startTimestep, endTimestep);
    selection.SetFieldIds(fieldIds);
    selection.SetMinUVWInM(minUVWInMeters);
    selection.SetMaxUVWInM(maxUVWInMeters);
    selection.SetEvenOrOddTimesteps(evenOddTimesteps);
    return selection;
  }

  bool IsSpectralFittingEnabled() const {
    return spectralFittingMode !=
           schaapcommon::fitters::SpectralFittingMode::kNoFitting;
  }

  size_t GetFeatherSize() const {
    if (featherSize) {
      return *featherSize;
    } else {
      // Return the default: 1% of sqrt(width * height)
      return std::ceil(std::sqrt(trimmedImageWidth * trimmedImageHeight) *
                       0.01);
    }
  }

  bool UseMpi() const { return nMpiNodes > 0; }

  bool IsBandSelected(size_t band_index) const {
    // An empty selection means that all bands are selected
    return spectralWindows.empty() ||
           spectralWindows.find(band_index) != spectralWindows.end();
  }

  /**
   * Determines if the gridder uses diagonal instrumental or full instrumental
   * polarizations.
   */
  aocommon::PolarizationEnum GetProviderPolarization(
      aocommon::PolarizationEnum entry_polarization) const;

  /**
   * True if either a h5parm solution file was specified or the beam is applied
   * while faceting.
   */
  bool UseFacetCorrections() const {
    return applyFacetBeam || !facetSolutionFiles.empty();
  }

 private:
  void checkPolarizations() const;
  bool determineReorder() const;
  std::string determineDataColumn(bool verbose) const;
  void logImportantSettings() const;
};

}  // namespace wsclean

#endif
