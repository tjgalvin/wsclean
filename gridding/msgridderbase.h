#ifndef MS_GRIDDER_BASE_H
#define MS_GRIDDER_BASE_H

#include "gridmode.h"

#include <aocommon/banddata.h>
#include <aocommon/image.h>
#include <aocommon/polarization.h>
#include <aocommon/imagecoordinates.h>

#include "../structures/observationinfo.h"
#include "../structures/msselection.h"
#include "../structures/weightmode.h"

#include "visibilitymodifier.h"
#include "visibilityweightingmode.h"

#include "../main/settings.h"

#include "../scheduling/metadatacache.h"
#include "../scheduling/griddingtaskmanager.h"

#include <aocommon/uvector.h>

#include <mutex>
#include <memory>

class MSReader;
namespace schaapcommon {
namespace h5parm {
class H5Parm;
class SolTab;
}  // namespace h5parm
}  // namespace schaapcommon

enum class PsfMode {
  kNone,    // Not a psf, grid the visibilities in the MS
  kSingle,  // Grid generated visibilities for a point source at the centre of
            // the main image
  kDirectionDependent  // Grid generated visibilities for a point source at the
                       // centre of the current facet
};

class MSGridderBase {
 public:
  MSGridderBase(const Settings& settings);
  virtual ~MSGridderBase();

  size_t ImageWidth() const { return _imageWidth; }
  size_t ImageHeight() const { return _imageHeight; }
  double ImagePadding() const { return _imagePadding; }
  double PixelSizeX() const { return _pixelSizeX; }
  double PixelSizeY() const { return _pixelSizeY; }
  size_t ActualWGridSize() const { return _actualWGridSize; }

  void ClearMeasurementSetList() {
    _measurementSets.clear();
    _selections.clear();
  }
  class MSProvider& MeasurementSet(size_t index) const {
    return *_measurementSets[index];
  }
  const MSSelection& Selection(size_t index) const {
    return _selections[index];
  }
  size_t MeasurementSetCount() const { return _measurementSets.size(); }
  void AddMeasurementSet(class MSProvider* msProvider,
                         const MSSelection& selection) {
    _measurementSets.push_back(msProvider);
    _selections.push_back(selection);
  }

  const std::string& DataColumnName() const { return _dataColumnName; }
  bool IsFacet() const { return _isFacet; }
  PsfMode GetPsfMode() const { return _psfMode; }
  bool DoSubtractModel() const { return _doSubtractModel; }
  bool SmallInversion() const { return _smallInversion; }
  aocommon::PolarizationEnum Polarization() const { return _polarization; }
  WeightMode Weighting() const { return _weighting; }
  const ImageWeights* GetImageWeights() const {
    return _precalculatedWeightInfo;
  }
  bool IsComplex() const { return _isComplex; }

  VisibilityWeightingMode GetVisibilityWeightingMode() const {
    return _visibilityWeightingMode;
  }
  bool StoreImagingWeights() const { return _storeImagingWeights; }

  void SetFacetIndex(size_t facetIndex) { _facetIndex = facetIndex; }
  void SetFacetGroupIndex(size_t index) { _facetGroupIndex = index; }
  /**
   * @brief In case of facet-based imaging, the model data in the @param
   * MSProvider is reset to zeros in every major cycle, and predicted data
   * should be add-assigned to the model data (_isFacet = true) rather
   * than overwriting it. For "standard" imaging, the model data should
   * be overwritten (_isFacet = false).
   */
  void SetIsFacet(bool isFacet) { _isFacet = isFacet; }
  void SetImageWidth(size_t imageWidth) { _imageWidth = imageWidth; }
  void SetImageHeight(size_t imageHeight) { _imageHeight = imageHeight; }
  void SetActualWGridSize(size_t actualWGridSize) {
    _actualWGridSize = actualWGridSize;
  }
  void SetPsfMode(PsfMode psf_mode) { _psfMode = psf_mode; }
  void SetImagePadding(double imagePadding) { _imagePadding = imagePadding; }
  void SetPolarization(aocommon::PolarizationEnum polarization) {
    _polarization = polarization;
  }
  void SetIsComplex(bool isComplex) { _isComplex = isComplex; }
  void SetDoSubtractModel(bool doSubtractModel) {
    _doSubtractModel = doSubtractModel;
  }

  void SetWriterLockManager(WriterLockManager* writerLockManager) {
    _writerLockManager = writerLockManager;
  }

  void SetImageWeights(const class ImageWeights* weights) {
    _precalculatedWeightInfo = weights;
  }

  /**
   * If this is the first gridder iteration, the gridder may output more
   * information.
   */
  bool IsFirstIteration() const { return _isFirstIteration; }
  void SetIsFirstIteration(bool isFirstIteration) {
    _isFirstIteration = isFirstIteration;
  }

  void SetStoreImagingWeights(bool storeImagingWeights) {
    _storeImagingWeights = storeImagingWeights;
  }

  virtual void Invert() = 0;

  virtual void Predict(std::vector<aocommon::Image>&& images) = 0;

  virtual std::vector<aocommon::Image> ResultImages() = 0;

  void SetPhaseCentreRA(const double phaseCentreRA) {
    _phaseCentreRA = phaseCentreRA;
  }
  void SetPhaseCentreDec(const double phaseCentreDec) {
    _phaseCentreDec = phaseCentreDec;
  }
  double PhaseCentreRA() const { return _phaseCentreRA; }
  double PhaseCentreDec() const { return _phaseCentreDec; }
  void SetLShift(const double l_shift) { _l_shift = l_shift; }
  void SetMShift(const double m_shift) { _m_shift = m_shift; }

  void SetMainImageDL(const double main_image_dl) {
    _mainImageDL = main_image_dl;
  }
  void SetMainImageDM(const double main_image_dm) {
    _mainImageDM = main_image_dm;
  }

  void SetFacetDirection(double ra, double dec) {
    _visibilityModifier.SetFacetDirection(ra, dec);
  }

  double FacetDirectionRA() const {
    return _visibilityModifier.FacetDirectionRA();
  }
  double FacetDirectionDec() const {
    return _visibilityModifier.FacetDirectionDec();
  }
  double LShift() const { return _l_shift; }
  double MShift() const { return _m_shift; }
  double MainImageDL() const { return _mainImageDL; }
  double MainImageDM() const { return _mainImageDM; }

  /**
   * Deallocate any data that is no longer necessary, but all methods
   * will still return results from the imaging, with the exception of
   * ImageReal/ImageResult().
   */
  virtual void FreeImagingData() {}

  GriddingKernelMode GetGridMode() const { return _gridMode; }
  void SetGridMode(GriddingKernelMode gridMode) { _gridMode = gridMode; }

  size_t TrimWidth() const { return _trimWidth; }
  size_t TrimHeight() const { return _trimHeight; }
  bool HasTrimSize() const { return _trimWidth != 0 || _trimHeight != 0; }
  void SetTrimSize(size_t trimWidth, size_t trimHeight) {
    _trimWidth = trimWidth;
    _trimHeight = trimHeight;
  }

  double StartTime() const { return _startTime; }
  bool HasDenormalPhaseCentre() const {
    return _l_shift != 0.0 || _m_shift != 0.0;
  }
  double ImageWeight() const { return _totalWeight; }
  double NormalizationFactor() const { return _totalWeight; }
  double BeamSize() const { return _theoreticalBeamSize; }

  /**
   * This is the sum of the weights as given by the measurement set, before the
   * image weighting is applied.
   */
  double VisibilityWeightSum() const { return _visibilityWeightSum; }
  /**
   * The number of visibilities that were gridded.
   */
  size_t GriddedVisibilityCount() const { return _griddedVisibilityCount; }
  /**
   * The maximum weight, after having applied the imaging weighting.
   */
  double MaxGriddedWeight() const { return _maxGriddedWeight; }
  /**
   * The effective number of visibilities, taking into account imaging weighting
   * and visibility weighting. This number is relative to the "best" visibility:
   * if one visibility with a weight of 10 and 5 visibilities with
   * a weight of 4 were gridded, the effective number of visibilities is
   * (10 + 5 x 4) / 10 = 3
   */
  double EffectiveGriddedVisibilityCount() const {
    return totalWeight() / MaxGriddedWeight();
  }

  void SetMetaDataCache(std::unique_ptr<MetaDataCache> cache) {
    _metaDataCache = std::move(cache);
  }
  std::unique_ptr<MetaDataCache> AcquireMetaDataCache() {
    return std::move(_metaDataCache);
  }

 protected:
  size_t ActualInversionWidth() const { return _actualInversionWidth; }
  size_t ActualInversionHeight() const { return _actualInversionHeight; }
  double ActualPixelSizeX() const { return _actualPixelSizeX; }
  double ActualPixelSizeY() const { return _actualPixelSizeY; }

  struct MSData {
   public:
    MSData();
    MSData(const MSData& source) = delete;
    ~MSData() = default;
    MSData& operator=(const MSData& source) = delete;

    MSProvider* msProvider;
    size_t msIndex;
    size_t dataDescId;
    aocommon::BandData bandData;
    size_t startChannel, endChannel;
    size_t matchingRows, totalRowsProcessed;
    double minW, maxW, maxWWithFlags, maxBaselineUVW, maxBaselineInM;
    size_t rowStart, rowEnd;
    double integrationTime;
    std::vector<std::string> antennaNames;

    aocommon::BandData SelectedBand() const {
      return aocommon::BandData(bandData, startChannel, endChannel);
    }
  };

  struct InversionRow {
    double uvw[3];
    size_t rowId;
    std::complex<float>* data;
  };

  /**
   * Initializes MS related data members, i.e. the @param _telescope and the
   * @param _pointResponse data in case a beam is applied on the facets and
   * EveryBeam is available and the @param _predictReader data member in case
   * @param isPredict is true.
   */
  void StartMeasurementSet(const MSGridderBase::MSData& msData, bool isPredict);

  /**
   * Read the visibilities from the msprovider, and apply weights and flags.
   *
   * This function applies both the selected method of visibility weighting
   * (i.e. the weights that are normally stored in the WEIGHT_SPECTRUM column)
   * and the imaging weight (coming from uniform or Briggs weighting, etc). To
   * read the data, this function requires scratch weight and model buffers to
   * store intermediate values in. Even if the caller does not need these
   * values, they still need to provide an already allocated buffer. This is to
   * avoid having to allocate memory within this method.
   * @tparam PolarizationCount Normally set to one when imaging a single
   * polarization, but set to 2 or 4 for IDG as it images multiple polarizations
   * at once.
   * @param msProvider The measurement set provider
   * @param rowData The resulting weighted data
   * @param curBand The spectral band currently being imaged
   * @param weightBuffer An allocated buffer to store intermediate weights in.
   * After returning from the call, these values will hold the full applied
   * weight (i.e. visibility weight * imaging weight).
   * @param modelBuffer An allocated buffer to store intermediate model data in.
   * @param isSelected Per visibility whether that visibility will be gridded in
   * this pass. When the visibility is not gridded, its weight will not be added
   * to the relevant sums (visibility count, weight sum, etc.).
   * @param gain_mode Selects which entry or entries in the gain matrix
   * (provided by EveryBeam and/or an h5 solution) file to use for correcting
   * the visibilities.
   */
  template <size_t PolarizationCount>
  void readAndWeightVisibilities(MSReader& msReader,
                                 const std::vector<std::string>& antennaNames,
                                 InversionRow& rowData,
                                 const aocommon::BandData& curBand,
                                 float* weightBuffer,
                                 std::complex<float>* modelBuffer,
                                 const bool* isSelected, GainMode gain_mode);

  /**
   * @brief Write (modelled) visibilities to MS, provides an interface to
   * MSProvider::WriteModel(). Method can be templated on the number of
   * polarizations (1, 2 or 4). The gain_mode can be used to
   * select an entry or entries from the gain matrix that should be used for the
   * correction.
   */
  template <size_t PolarizationCount>
  void writeVisibilities(MSProvider& msProvider,
                         const std::vector<std::string>& antennaNames,
                         const aocommon::BandData& curBand,
                         std::complex<float>* buffer, GainMode gain_mode);

  double _maxW, _minW;

  virtual size_t getSuggestedWGridSize() const = 0;

  void resetVisibilityCounters() {
    _griddedVisibilityCount = 0;
    _totalWeight = 0.0;
    _maxGriddedWeight = 0.0;
    _visibilityWeightSum = 0.0;
  }

  double totalWeight() const { return _totalWeight; }

  void initializeMSDataVector(std::vector<MSData>& msDataVector);

  std::unique_ptr<MetaDataCache> _metaDataCache;

  template <size_t PolarizationCount>
  static void rotateVisibilities(const aocommon::BandData& bandData,
                                 double shiftFactor,
                                 std::complex<float>* dataIter);

  const Settings& _settings;

 private:
  static std::vector<std::string> getAntennaNames(
      const casacore::MSAntenna& msAntenna);

  void resetMetaData() { _hasFrequencies = false; }

  void calculateMSLimits(const aocommon::BandData& selectedBand,
                         double startTime) {
    if (_hasFrequencies) {
      _freqLow = std::min(_freqLow, selectedBand.LowestFrequency());
      _freqHigh = std::max(_freqHigh, selectedBand.HighestFrequency());
      _bandStart = std::min(_bandStart, selectedBand.BandStart());
      _bandEnd = std::max(_bandEnd, selectedBand.BandEnd());
      _startTime = std::min(_startTime, startTime);
    } else {
      _freqLow = selectedBand.LowestFrequency();
      _freqHigh = selectedBand.HighestFrequency();
      _bandStart = selectedBand.BandStart();
      _bandEnd = selectedBand.BandEnd();
      _startTime = startTime;
      _hasFrequencies = true;
    }
  }

  template <size_t NPolInMSProvider>
  void calculateWLimits(MSGridderBase::MSData& msData);

  void initializeMeasurementSet(MSGridderBase::MSData& msData,
                                MetaDataCache::Entry& cacheEntry,
                                bool isCacheInitialized);

  void calculateOverallMetaData(const MSData* msDataVector);
  bool hasWGridSize() const { return _wGridSize != 0; }
  void initializeBandData(const casacore::MeasurementSet& ms,
                          MSGridderBase::MSData& msData);
  void initializePointResponse(const MSGridderBase::MSData& msData);
  void initializePredictReader(MSProvider& msProvider);

  template <size_t PolarizationCount, GainMode GainEntry>
  void readAndWeightVisibilities(MSReader& msReader,
                                 const std::vector<std::string>& antennaNames,
                                 InversionRow& rowData,
                                 const aocommon::BandData& curBand,
                                 float* weightBuffer,
                                 std::complex<float>* modelBuffer,
                                 const bool* isSelected);

  template <size_t PolarizationCount, GainMode GainEntry>
  void writeVisibilities(MSProvider& msProvider,
                         const std::vector<std::string>& antennaNames,
                         const aocommon::BandData& curBand,
                         std::complex<float>* buffer);

  /**
   * @brief Applies both the conjugated h5 parm
   * solutions to the visibilities and computes the weight corresponding to the
   * combined effect.
   *
   * @param apply_forward If true, also apply the forward (non-conjugated) gain.
   *                      Used for generating a direction dependent psf, where
   *                      both the (forward) gain needs to be applied for the
   *                      predict/degridding step
   *                      and the conjugate gain for the gridding step
   */
  template <size_t PolarizationCount, GainMode GainEntry>
  void ApplyConjugatedH5Parm(MSReader& msReader,
                             const std::vector<std::string>& antennaNames,
                             InversionRow& rowData,
                             const aocommon::BandData& curBand,
                             const float* weightBuffer,
                             bool apply_forward = false);

#ifdef HAVE_EVERYBEAM
  /**
   * @brief Applies the conjugated facet beam to the visibilities and computes
   * the weight corresponding to the combined effect.
   *
   * @param apply_forward If true, also apply the forward (non-conjugated) gain.
   *                      Used for generating a direction dependent psf, where
   *                      both the (forward) gain needs to be applied for the
   *                      predict/degridding step
   *                      and the conjugate gain for the gridding step
   */

  template <size_t PolarizationCount, GainMode GainEntry>
  void ApplyConjugatedFacetBeam(MSReader& msReader, InversionRow& rowData,
                                const aocommon::BandData& curBand,
                                const float* weightBuffer,
                                bool apply_forward = false);

  /**
   * @brief Applies both the conjugated facet beam and the conjugated h5 parm
   * solutions to the visibilities and computes the weight corresponding to the
   * combined effect.
   *
   * @param apply_forward If true, also apply the forward (non-conjugated) gain.
   *                      Used for generating a direction dependent psf, where
   *                      both the (forward) gain needs to be applied for the
   *                      predict/degridding step
   *                      and the conjugate gain for the gridding step
   */
  template <size_t PolarizationCount, GainMode GainEntry>
  void ApplyConjugatedFacetDdEffects(
      MSReader& msReader, const std::vector<std::string>& antennaNames,
      InversionRow& rowData, const aocommon::BandData& curBand,
      const float* weightBuffer, bool apply_forward = false);
#endif  // HAVE_EVERYBEAM

  size_t _actualInversionWidth;
  size_t _actualInversionHeight;
  double _actualPixelSizeX;
  double _actualPixelSizeY;
  double _phaseCentreRA, _phaseCentreDec, _l_shift, _m_shift;
  double _mainImageDL, _mainImageDM;
  size_t _facetIndex;
  /// @p _facetGroupIndex and @p _msIndex in conjunction with the @p
  /// MeasurementSetCount() determine the index in the _writerGroupLocks vector,
  /// having size FacetGroupCount() * MeasurementSetCount(). These variable are
  /// only relevant for prediction.
  size_t _facetGroupIndex;
  size_t _msIndex;
  /// @see SetIsFacet()
  bool _isFacet;
  double _imagePadding;
  size_t _imageWidth, _imageHeight;
  size_t _trimWidth, _trimHeight;
  double _pixelSizeX, _pixelSizeY;
  size_t _wGridSize, _actualWGridSize;
  std::vector<MSProvider*> _measurementSets;
  std::string _dataColumnName;
  PsfMode _psfMode;
  bool _doSubtractModel, _smallInversion;
  double _wLimit;
  const class ImageWeights* _precalculatedWeightInfo;
  aocommon::PolarizationEnum _polarization;
  bool _isComplex;
  WeightMode _weighting;
  bool _isFirstIteration;
  std::vector<MSSelection> _selections;
  VisibilityWeightingMode _visibilityWeightingMode;
  GriddingKernelMode _gridMode;
  bool _storeImagingWeights;
  double _theoreticalBeamSize;

  bool _hasFrequencies;
  double _freqHigh, _freqLow;
  double _bandStart, _bandEnd;
  double _startTime;

  size_t _griddedVisibilityCount;
  double _totalWeight;
  double _maxGriddedWeight;
  double _visibilityWeightSum;

  aocommon::UVector<float> _scratchImageWeights;

  std::unique_ptr<MSReader> _predictReader;
  WriterLockManager* _writerLockManager;

  VisibilityModifier _visibilityModifier;
};

template <size_t PolarizationCount>
inline void MSGridderBase::readAndWeightVisibilities(
    MSReader& msReader, const std::vector<std::string>& antennaNames,
    InversionRow& rowData, const aocommon::BandData& curBand,
    float* weightBuffer, std::complex<float>* modelBuffer,
    const bool* isSelected, GainMode gain_mode) {
  switch (gain_mode) {
    case GainMode::kXX:
      readAndWeightVisibilities<PolarizationCount, GainMode::kXX>(
          msReader, antennaNames, rowData, curBand, weightBuffer, modelBuffer,
          isSelected);
      break;
    case GainMode::kYY:
      readAndWeightVisibilities<PolarizationCount, GainMode::kYY>(
          msReader, antennaNames, rowData, curBand, weightBuffer, modelBuffer,
          isSelected);
      break;
    case GainMode::kDiagonal:
      readAndWeightVisibilities<PolarizationCount, GainMode::kDiagonal>(
          msReader, antennaNames, rowData, curBand, weightBuffer, modelBuffer,
          isSelected);
      break;
    case GainMode::kFull:
      if constexpr (PolarizationCount == 2 || PolarizationCount == 4)
        readAndWeightVisibilities<PolarizationCount, GainMode::kFull>(
            msReader, antennaNames, rowData, curBand, weightBuffer, modelBuffer,
            isSelected);
      else
        throw std::runtime_error(
            "Invalid combination of visibility polarizations and gain mode");
      break;
  }
}

template <size_t PolarizationCount>
inline void MSGridderBase::writeVisibilities(
    MSProvider& msProvider, const std::vector<std::string>& antennaNames,
    const aocommon::BandData& curBand, std::complex<float>* buffer,
    GainMode gain_mode) {
  switch (gain_mode) {
    case GainMode::kXX:
      writeVisibilities<PolarizationCount, GainMode::kXX>(
          msProvider, antennaNames, curBand, buffer);
      break;
    case GainMode::kYY:
      writeVisibilities<PolarizationCount, GainMode::kYY>(
          msProvider, antennaNames, curBand, buffer);
      break;
    case GainMode::kDiagonal:
      writeVisibilities<PolarizationCount, GainMode::kDiagonal>(
          msProvider, antennaNames, curBand, buffer);
      break;
    case GainMode::kFull:
      if constexpr (PolarizationCount == 2 || PolarizationCount == 4)
        writeVisibilities<PolarizationCount, GainMode::kFull>(
            msProvider, antennaNames, curBand, buffer);
      else
        throw std::runtime_error(
            "Invalid combination of visibility polarizations and gain mode");
      break;
  }
}

#endif
