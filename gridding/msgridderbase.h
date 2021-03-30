#ifndef MS_GRIDDER_BASE_H
#define MS_GRIDDER_BASE_H

#include "measurementsetgridder.h"

#include "../main/settings.h"

#include "../scheduling/metadatacache.h"

#include "../structures/multibanddata.h"

#ifdef HAVE_EVERYBEAM
#include <EveryBeam/pointresponse/pointresponse.h>
#endif

#include <aocommon/uvector.h>

#include <mutex>
#include <memory>

class MSGridderBase : public MeasurementSetGridder {
 public:
  MSGridderBase(const Settings& settings);
  ~MSGridderBase();

  virtual double StartTime() const final override { return _startTime; }
  virtual bool HasDenormalPhaseCentre() const final override {
    return _phaseCentreDL != 0.0 || _phaseCentreDM != 0.0;
  }
  virtual double ImageWeight() const final override { return _totalWeight; }
  virtual double NormalizationFactor() const final override {
    return _totalWeight;
  }
  virtual double BeamSize() const final override {
    return _theoreticalBeamSize;
  }

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
  int64_t getAvailableMemory(double memFraction, double absMemLimit);

  struct MSData {
   public:
    MSData();
    ~MSData();
    class MSProvider* msProvider;
    size_t msIndex;
    MultiBandData bandData;
    size_t startChannel, endChannel;
    size_t matchingRows, totalRowsProcessed;
    double minW, maxW, maxWWithFlags, maxBaselineUVW, maxBaselineInM;
    size_t rowStart, rowEnd;
    double integrationTime;

    MultiBandData SelectedBand() const {
      return MultiBandData(bandData, startChannel, endChannel);
    }

   private:
    MSData(const MSData& source);

    void operator=(const MSData& source);
  };

  struct InversionRow {
    double uvw[3];
    size_t dataDescId, rowId;
    std::complex<float>* data;
  };

  void resetMetaData() { _hasFrequencies = false; }

  void calculateMSLimits(const MultiBandData& selectedBand, double startTime) {
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
   * polarization, but set to 4 for IDG as it images all polarizations at once.
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
   */
  template <size_t PolarizationCount>
  void readAndWeightVisibilities(MSProvider& msProvider, InversionRow& rowData,
                                 const BandData& curBand, float* weightBuffer,
                                 std::complex<float>* modelBuffer,
                                 const bool* isSelected);

  /**
   * @brief Write (modelled) visibilities to MS, provides an interface to
   * MSProvider::WriteModel()
   */
  void writeVisibilities(MSProvider& msProvider, size_t rowId,
                         const std::complex<float>* buffer) const;

  double _maxW, _minW;
  double _theoreticalBeamSize;
  size_t _actualInversionWidth, _actualInversionHeight;
  double _actualPixelSizeX, _actualPixelSizeY;

  virtual size_t getSuggestedWGridSize() const = 0;

  void resetVisibilityCounters() {
    _griddedVisibilityCount = 0;
    _totalWeight = 0.0;
    _maxGriddedWeight = 0.0;
    _visibilityWeightSum = 0.0;
  }

  double totalWeight() const { return _totalWeight; }

  void initializeMSDataVector(std::vector<MSData>& msDataVector);

  std::unique_ptr<struct MetaDataCache> _metaDataCache;

  template <size_t PolarizationCount>
  static void rotateVisibilities(const BandData& bandData, double shiftFactor,
                                 std::complex<float>* dataIter);

  const Settings& _settings;

 private:
  void initializeBandData(casacore::MeasurementSet& ms,
                          MSGridderBase::MSData& msData);

  bool _hasFrequencies;
  double _freqHigh, _freqLow;
  double _bandStart, _bandEnd;
  double _startTime;

  size_t _griddedVisibilityCount;
  double _totalWeight;
  double _maxGriddedWeight;
  double _visibilityWeightSum;

  aocommon::UVector<float> _scratchWeights;

#ifdef HAVE_EVERYBEAM
  // _telescope attribute needed to keep the telecope in _point_response alive
  std::unique_ptr<everybeam::telescope::Telescope> _telescope;
  std::unique_ptr<everybeam::pointresponse::PointResponse> _pointResponse;
  aocommon::UVector<std::complex<float>> _cachedResponse;
#endif
};

#endif
