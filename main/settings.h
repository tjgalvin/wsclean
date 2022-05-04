#ifndef WSCLEAN_SETTINGS_H
#define WSCLEAN_SETTINGS_H

#include "../gridding/wstackinggridder.h"

#include "../gridding/visibilityweightingmode.h"
#include "../structures/weightmode.h"

#include "../structures/msselection.h"

#include <aocommon/system.h>

#include <schaapcommon/fitters/spectralfitter.h>

#include <radler/settings.h>

enum class DirectFTPrecision { Float, Double, LongDouble };

/**
 * This class describes all settings for a single WSClean run.
 * @sa WSClean
 */
class Settings {
 public:
  Settings();

  void Validate() const;

  void Propagate(bool verbose = true);

  void RecalculatePaddedDimensions(bool verbose = true);

  std::vector<std::string> filenames;
  enum Mode { ImagingMode, PredictMode, RestoreMode, RestoreListMode } mode;
  size_t paddedImageWidth, paddedImageHeight;
  size_t trimmedImageWidth, trimmedImageHeight;
  bool hasShift;
  double shiftRA, shiftDec;
  double imagePadding;
  size_t widthForNWCalculation, heightForNWCalculation;
  size_t channelsOut, intervalsOut;
  enum MSSelection::EvenOddSelection evenOddTimesteps;
  bool divideChannelsByGaps;
  aocommon::UVector<double> divideChannelFrequencies;
  double pixelScaleX, pixelScaleY;
  std::string restoreModel, restoreInput, restoreOutput;
  double manualBeamMajorSize, manualBeamMinorSize, manualBeamPA;
  bool fittedBeam, theoreticBeam, circularBeam;
  double beamFittingBoxSize;
  bool continuedRun;
  double memFraction, absMemLimit;
  double minUVWInMeters, maxUVWInMeters, minUVInLambda, maxUVInLambda, wLimit,
      rankFilterLevel;
  size_t rankFilterSize;
  double gaussianTaperBeamSize, tukeyTaperInLambda, tukeyInnerTaperInLambda,
      edgeTaperInLambda, edgeTukeyTaperInLambda;
  bool useWeightsAsTaper;
  size_t nWLayers;
  double nWLayersFactor;
  size_t antialiasingKernelSize, overSamplingFactor, threadCount,
      parallelReordering, parallelGridding;
  bool useMPI, masterDoesWork;
  std::vector<size_t> fieldIds;
  size_t startTimestep, endTimestep;
  size_t startChannel, endChannel;
  size_t predictionChannels;
  std::string dataColumnName;
  std::set<aocommon::PolarizationEnum> polarizations;
  std::string facetRegionFilename;
  std::set<size_t> spectralWindows;
  WeightMode weightMode;
  std::string prefixName;
  bool joinedPolarizationDeconvolution;
  bool joinedFrequencyDeconvolution;
  bool smallInversion, makePSF, makePSFOnly, isWeightImageSaved, isUVImageSaved,
      isDirtySaved, isFirstResidualSaved;
  bool reusePsf, reuseDirty;
  std::string reusePsfPrefix, reuseDirtyPrefix;
  bool writeImagingWeightSpectrumColumn;
  std::string temporaryDirectory;
  bool forceReorder, forceNoReorder, doReorder;
  bool subtractModel, modelUpdateRequired, mfWeighting;
  size_t fullResOffset, fullResWidth, fullResPad;
  std::string beamModel;
  std::string beamMode;
  std::string beamNormalisationMode;
  bool applyPrimaryBeam, reusePrimaryBeam, savePsfPb;
  double primaryBeamLimit;
  std::string mwaPath;
  size_t primaryBeamGridSize, primaryBeamUpdateTime;
  bool directFT;
  DirectFTPrecision directFTPrecision;
  bool useIDG, useWGridder;
  double wgridderAccuracy;
  std::string atermConfigFilename;
  double atermKernelSize;
  bool gridWithBeam;
  double beamAtermUpdateTime;  // in seconds.
  std::vector<std::string> facetSolutionFiles;
  std::vector<std::string> facetSolutionTables;
  bool applyFacetBeam;
  double facetBeamUpdateTime;  // in seconds.
  bool saveATerms;
  enum IDGMode { IDG_DEFAULT, IDG_GPU, IDG_CPU, IDG_HYBRID } idgMode;
  enum GridMode gridMode;
  enum VisibilityWeightingMode visibilityWeightingMode;
  double baselineDependentAveragingInWavelengths;
  bool simulateNoise;
  double simulatedNoiseStdDev;
  std::string simulatedBaselineNoiseFilename;

  /** @{
   * These settings all relate to the deconvolution.
   */
  std::set<aocommon::PolarizationEnum> linkedPolarizations;
  size_t parallelDeconvolutionMaxSize;
  size_t parallelDeconvolutionMaxThreads;
  double deconvolutionThreshold;
  double deconvolutionGain;
  double deconvolutionMGain;
  bool autoDeconvolutionThreshold, autoMask;
  double autoDeconvolutionThresholdSigma, autoMaskSigma;
  double localRMSWindow;
  radler::LocalRmsMethod localRMSMethod;
  bool saveSourceList;
  size_t deconvolutionIterationCount, majorIterationCount;
  bool allowNegativeComponents, stopOnNegativeComponents;
  bool useSubMinorOptimization, squaredJoins;
  double spectralCorrectionFrequency;
  std::vector<float> spectralCorrection;
  bool multiscaleFastSubMinorLoop;
  double multiscaleGain, multiscaleDeconvolutionScaleBias;
  size_t multiscaleMaxScales;
  double multiscaleConvolutionPadding;
  aocommon::UVector<double> multiscaleScaleList;
  radler::MultiscaleShape multiscaleShapeFunction;

  double deconvolutionBorderRatio;
  std::string fitsDeconvolutionMask, casaDeconvolutionMask;
  bool horizonMask;
  double horizonMaskDistance;
  std::string localRMSImage;
  std::string pythonDeconvolutionFilename;
  bool iuwtSNRTest;
  std::string moreSaneLocation, moreSaneArgs;
  aocommon::UVector<double> moreSaneSigmaLevels;
  schaapcommon::fitters::SpectralFittingMode spectralFittingMode;
  size_t spectralFittingTerms;
  std::string forcedSpectrumFilename;
  /**
   * The number of channels used during deconvolution. This can be used to
   * image with more channels than deconvolution. Before deconvolution,
   * channels are averaged, and after deconvolution they are interpolated.
   * It is 0 when all channels should be used.
   */
  size_t deconvolutionChannelCount;

  /**
   * Type of deconvolution algorithm.
   */
  radler::AlgorithmType algorithmType;
  /**
   * @}
   */

  /**
   * @brief Extract the settings that are relevant to the deconvolution.
   * Currently, it duplicates the existing settings into a DeconvolutionSettings
   * object.
   */
  radler::Settings GetRadlerSettings() const;

  MSSelection GetMSSelection() const {
    MSSelection selection;
    selection.SetInterval(startTimestep, endTimestep);
    selection.SetFieldIds(fieldIds);
    selection.SetMinUVWInM(minUVWInMeters);
    selection.SetMaxUVWInM(maxUVWInMeters);
    selection.SetEvenOrOddTimesteps(evenOddTimesteps);
    return selection;
  }

  bool IsSpectralFittingEnabled() const {
    return spectralFittingMode !=
           schaapcommon::fitters::SpectralFittingMode::NoFitting;
  }

 private:
  void checkPolarizations() const;
  bool determineReorder() const;
  std::string determineDataColumn(bool verbose) const;
  void logImportantSettings() const;
};

inline Settings::Settings()
    : filenames(),
      mode(ImagingMode),
      paddedImageWidth(0),
      paddedImageHeight(0),
      trimmedImageWidth(0),
      trimmedImageHeight(0),
      hasShift(false),
      shiftRA(0.0),
      shiftDec(0.0),
      imagePadding(1.2),
      widthForNWCalculation(0),
      heightForNWCalculation(0),
      channelsOut(1),
      intervalsOut(1),
      evenOddTimesteps(MSSelection::AllTimesteps),
      divideChannelsByGaps(false),
      divideChannelFrequencies(),
      pixelScaleX(0.0),
      pixelScaleY(0.0),
      restoreModel(),
      restoreInput(),
      restoreOutput(),
      manualBeamMajorSize(0.0),
      manualBeamMinorSize(0.0),
      manualBeamPA(0.0),
      fittedBeam(true),
      theoreticBeam(false),
      circularBeam(false),
      beamFittingBoxSize(10.0),
      continuedRun(false),
      memFraction(1.0),
      absMemLimit(0.0),
      minUVWInMeters(0.0),
      maxUVWInMeters(0.0),
      minUVInLambda(0.0),
      maxUVInLambda(0.0),
      wLimit(0.0),
      rankFilterLevel(3.0),
      rankFilterSize(16),
      gaussianTaperBeamSize(0.0),
      tukeyTaperInLambda(0.0),
      tukeyInnerTaperInLambda(0.0),
      edgeTaperInLambda(0.0),
      edgeTukeyTaperInLambda(0.0),
      useWeightsAsTaper(false),
      nWLayers(0),
      nWLayersFactor(1.0),
      antialiasingKernelSize(7),
      overSamplingFactor(1023),
      threadCount(aocommon::system::ProcessorCount()),
      parallelReordering(1),
      parallelGridding(1),
      useMPI(false),
      masterDoesWork(true),
      fieldIds{0},
      startTimestep(0),
      endTimestep(0),
      startChannel(0),
      endChannel(0),
      predictionChannels(0),
      dataColumnName(),
      polarizations({aocommon::Polarization::StokesI}),
      facetRegionFilename(),
      weightMode(WeightMode::UniformWeighted),
      prefixName("wsclean"),
      joinedPolarizationDeconvolution(false),
      joinedFrequencyDeconvolution(false),
      smallInversion(true),
      makePSF(false),
      makePSFOnly(false),
      isWeightImageSaved(false),
      isUVImageSaved(false),
      isDirtySaved(true),
      isFirstResidualSaved(false),
      reusePsf(false),
      reuseDirty(false),
      reusePsfPrefix(),
      reuseDirtyPrefix(),
      writeImagingWeightSpectrumColumn(false),
      temporaryDirectory(),
      forceReorder(false),
      forceNoReorder(false),
      doReorder(true),
      subtractModel(false),
      modelUpdateRequired(true),
      mfWeighting(false),
      fullResOffset(0),
      fullResWidth(0),
      fullResPad(0),
      beamModel(),
      beamMode("default"),
      beamNormalisationMode("preapplied"),
      applyPrimaryBeam(false),
      reusePrimaryBeam(false),
      savePsfPb(false),
      primaryBeamLimit(0.005),
      primaryBeamGridSize(32),
      primaryBeamUpdateTime(1800),
      directFT(false),
      directFTPrecision(DirectFTPrecision::Double),
      useIDG(false),
      useWGridder(false),
      wgridderAccuracy(1e-4),
      atermConfigFilename(),
      atermKernelSize(5.0),
      gridWithBeam(false),
      beamAtermUpdateTime(300.0),
      facetSolutionFiles(),
      facetSolutionTables(),
      applyFacetBeam(false),
      facetBeamUpdateTime(120.0),
      saveATerms(false),
      idgMode(IDG_DEFAULT),
      gridMode(GridMode::KaiserBesselKernel),
      visibilityWeightingMode(
          VisibilityWeightingMode::NormalVisibilityWeighting),
      baselineDependentAveragingInWavelengths(0.0),
      simulateNoise(false),
      simulatedNoiseStdDev(0.0),
      simulatedBaselineNoiseFilename(),
      // Deconvolution default settings:
      linkedPolarizations(),
      parallelDeconvolutionMaxSize(0),
      parallelDeconvolutionMaxThreads(0),
      deconvolutionThreshold(0.0),
      deconvolutionGain(0.1),
      deconvolutionMGain(1.0),
      autoDeconvolutionThreshold(false),
      autoMask(false),
      autoDeconvolutionThresholdSigma(0.0),
      autoMaskSigma(0.0),
      localRMSWindow(25.0),
      localRMSMethod(radler::LocalRmsMethod::kNone),
      saveSourceList(false),
      deconvolutionIterationCount(0),
      majorIterationCount(20),
      allowNegativeComponents(true),
      stopOnNegativeComponents(false),
      useSubMinorOptimization(true),
      squaredJoins(false),
      spectralCorrectionFrequency(0.0),
      spectralCorrection(),
      multiscaleFastSubMinorLoop(true),
      multiscaleGain(0.2),
      multiscaleDeconvolutionScaleBias(0.6),
      multiscaleMaxScales(0),
      multiscaleConvolutionPadding(1.1),
      multiscaleScaleList(),
      multiscaleShapeFunction(radler::MultiscaleShape::kTaperedQuadraticShape),
      deconvolutionBorderRatio(0.0),
      fitsDeconvolutionMask(),
      casaDeconvolutionMask(),
      horizonMask(false),
      horizonMaskDistance(0.0),
      pythonDeconvolutionFilename(),
      iuwtSNRTest(false),
      moreSaneLocation(),
      moreSaneArgs(),
      spectralFittingMode(
          schaapcommon::fitters::SpectralFittingMode::NoFitting),
      spectralFittingTerms(0),
      forcedSpectrumFilename(),
      deconvolutionChannelCount(0),
      algorithmType(radler::AlgorithmType::kGenericClean) {}

#endif
