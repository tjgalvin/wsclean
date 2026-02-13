#ifndef WSCLEAN_GRIDDING_MS_GRIDDER_H_
#define WSCLEAN_GRIDDING_MS_GRIDDER_H_

#include "gridmode.h"

#include "../structures/weightmode.h"

#include "msgridderdata.h"
#include "msprovidercollection.h"

#include "../main/settings.h"

#include <aocommon/uvector.h>

namespace wsclean {

/**
 * MsGridder implements an interface that gridders can implement in order to do
 * gridding of one or multiple measurement sets. In addition to the interface,
 * MsGridder also stores various data for gridding and methods for calculating
 * and working with the data.
 *
 * Data/Methods relating to an individual measurement set are provided by @ref
 * MsGridderData while those that are more general and the gridding interface
 * itself reside on MsGridder.
 *
 * Everything required to prepare data for gridding can be accessed via @ref
 * MsGridderData while @ref MsGridder is needed for the actual gridding
 * operation. E.g., when doing shared reads, @ref MsGridderManager instantiates
 * a single @ref MsGridderData, prepares a single shared data buffer and then
 * passes that to the MsGridders for the gridding to be performed.
 */
class MsGridder : public MsGridderData {
 public:
  MsGridder(const Settings& settings,
            MsProviderCollection& ms_provider_collection);
  ~MsGridder() override;

  size_t ImageWidth() const { return image_width_; }
  size_t ImageHeight() const { return image_height_; }
  double ImagePadding() const { return image_padding_; }

  size_t ActualWGridSize() const { return actual_w_grid_size_; }

  const std::string& DataColumnName() const { return data_column_name_; }
  bool SmallInversion() const { return small_inversion_; }
  WeightMode Weighting() const { return weighting_; }
  bool IsComplex() const { return is_complex_; }

  void SetFacetIndex(size_t facet_index) { facet_index_ = facet_index; }

  void SetImageWidth(size_t image_width) { image_width_ = image_width; }
  void SetImageHeight(size_t image_height) { image_height_ = image_height; }
  void SetActualWGridSize(size_t actual_w_grid_size) {
    actual_w_grid_size_ = actual_w_grid_size;
  }
  void SetImagePadding(double image_padding) { image_padding_ = image_padding; }
  void SetIsComplex(bool is_complex) { is_complex_ = is_complex; }

  /**
   * When processing the first gridder task, the gridder may output more
   * information.
   */
  bool IsFirstTask() const { return is_first_task_; }
  void SetIsFirstTask(bool is_first_task) { is_first_task_ = is_first_task; }

  /**
   * To handle inversion a gridder must implement the following
   * functions.
   * If a gridder only makes one pass per MS:
   *   @ref StartInversion()
   *   @ref GridMeasurementSet()
   *   @ref  FinishInversion()
   * If multiple passes then additionally:
   *   @ref GetNInversionPasses()
   *   @ref StartInversionPass()
   *   @ref FinishInversionPass()
   * For inversion with shared memory and reads:
   *   @ref CalculateConstantMemory()
   *   @ref CalculateMaxRowsInMemory()
   *   @ref GridSharedMeasurementSetChunk()
   */
  virtual void StartInversion() = 0;
  virtual size_t GetNInversionPasses() const { return 1; }
  virtual void StartInversionPass(size_t pass_index){};
  /** @return The number of visibility rows processed */
  virtual size_t GridMeasurementSet(
      const MsProviderCollection::MsData& ms_data) = 0;
  virtual void FinishInversionPass(size_t pass_index){};
  virtual void FinishInversion() = 0;
  /** @return Constant memory overhead for the gridder in bytes. i.e. Not
   * including per row memory usage. */
  virtual size_t CalculateConstantMemory() const { return 0; }
  /** Calculate the maximum amount of rows the gridder expects to be able to fit
   * in memory given a specific available memory size and constant memory
   * overhead.
   * @param additional_per_row_consumption External consumption per row that
   * also needs to be taken into account; e.g. for shared reads the gridder
   * manager may need to cache antenna pairs and correction time offsets.
   * @param num_polarizations_stored The number of visibility polarizations per
   * channel that will be kept in memory. Often this will be 1 as we collapse
   * the visibilities before gridding, but in some cases; e.g. shared reads we
   * need to keep all the polatizations in memory
   * @param additional_per_row_consumption Is taken as a double in order to
   * maintain precision, because in some cases its value can represent shared
   * memory across multiple facets (shared_per_row_consumption/num_facets) and
   * therefore might be fractional.
   * @return The number of rows that will fit in memory.
   */
  virtual size_t CalculateMaxRowsInMemory(
      int64_t available_memory, size_t constant_memory,
      double additional_per_row_consumption, size_t per_row_uvw_consumption,
      size_t channel_count, size_t num_polarizations_stored) const {
    return 0;
  }
  /** Takes a chunk of visibilities pre-populated by the caller, that does
   * not have facet solution pre-applied but for which the correction sums
   * have already been computed, and performs gridding on the chunk.
   * Corrections will be applied "on the fly" during the gridding.
   * In addition to the visibilities, various additional data is required in
   * order to compute the corrections as needed.
   *
   * It is expected that the average correction has already been calculated by
   * summing the corrections via @ref ApplyCorrections<kSum>() and that this
   * gridder is already populated with the resulting sums.
   *
   * @param n_polarizations The number of polarizations per visibility in @ref
   * visibilities
   * @param antennas Pointer to n_rows `std::pair<size_t, size_t>` containing
   * the antenna pair for each row of visibilities.
   * @param uvws Pointer to n_rows*3 doubles containing UVW in meters.
   *        U(row) := uvws[3*row  ]
   *        V(row) := uvws[3*row+1]
   *        W(row) := uvws[3*row+2]
   * @param frequencies Pointer to n_channels doubles containing channel
   * frequencies
   * @param visibilities Pointer to n_rows*n_channels*n_polarizations
   * `complex<float>` containing weighted but uncorrected and not yet collapsed
   * visibilities: visibility(row, chan) := vis[row*n_chan + chan]
   * @param time_offsets Pointer to n_rows `size_t` containing the time offset
   * as calculated by @ref CacheParmResponse() for the corresponding visibility
   * row when applying @ref ApplyCorrections<ModifierBehaviour::kSum>() on it
   * For further explanation see @ref VisibilityCallbackBuffer::time_offsets_
   */
  virtual void GridSharedMeasurementSetChunk(
      bool apply_corrections, size_t n_polarizations, size_t n_rows,
      const double* uvws, const double* frequencies,
      const aocommon::BandData& selected_band, size_t data_desc_id,
      const std::pair<size_t, size_t>* antennas,
      const std::complex<float>* visibilities, const size_t* time_offsets,
      size_t n_antennas, const std::vector<std::complex<float>>& parm_response,
      const BeamResponseCacheChunk& beam_response) {
    throw std::runtime_error("Gridder does not yet support shared reading");
  }

  /**
   * To handle prediction a gridder must implement the following
   * functions.
   * If a gridder only makes one pass per MS:
   *   @ref StartPredict()
   *   @ref PredictMeasurementSet()
   *   @ref FinishPredict()
   * If multiple passes then additionally:
   *   @ref GetNPredictPasses()
   *   @ref StartPredictPass()
   *   @ref FinishPredictPass()
   * For prediction with shared memory and writes:
   *   @ref PredictChunk()
   */
  virtual void StartPredict(std::vector<aocommon::Image>&& images) = 0;
  virtual size_t GetNPredictPasses() const { return 1; }
  virtual void StartPredictPass(size_t pass_index){};
  /** @return The number of visibility rows processed */
  virtual size_t PredictMeasurementSet(
      const MsProviderCollection::MsData& ms_data) = 0;
  virtual void PredictChunk(size_t n_rows, size_t n_channels,
                            const double* frequencies, const double* uvws,
                            std::complex<float>* visibilities) const {
    throw std::runtime_error(
        "Selected gridder does not yet support shared predict. Try a different "
        "gridder, e.g. wgridder.");
  }
  virtual void FinishPredictPass(size_t pass_index){};
  virtual void FinishPredict() = 0;

  virtual std::vector<aocommon::Image> ResultImages() = 0;

  void SetPhaseCentreRA(const double phase_centre_ra) {
    phase_centre_ra_ = phase_centre_ra;
  }
  void SetPhaseCentreDec(const double phase_centre_dec) {
    phase_centre_dec_ = phase_centre_dec;
  }
  double PhaseCentreRA() const { return phase_centre_ra_; }
  double PhaseCentreDec() const { return phase_centre_dec_; }

  /**
   * Deallocate any data that is no longer necessary, but all methods
   * will still return results from the imaging, with the exception of
   * ImageReal/ImageResult().
   */
  virtual void FreeImagingData() {}

  GriddingKernelMode GetGridMode() const { return grid_mode_; }
  void SetGridMode(GriddingKernelMode grid_mode) { grid_mode_ = grid_mode; }

  size_t TrimWidth() const { return trim_width_; }
  size_t TrimHeight() const { return trim_height_; }
  bool HasTrimSize() const { return trim_width_ != 0 || trim_height_ != 0; }
  void SetTrimSize(size_t trim_width, size_t trim_height) {
    trim_width_ = trim_width;
    trim_height_ = trim_height;
  }

  /**
   * @return The normalization factor, which is always equal to
   * the image weight in the current implementation.
   *
   * This interface has separate ImageWeight and NormalizationFactor functions
   * since they are conceptually different and the implementation of
   * NormalizationFactor may change in the future.
   */
  double NormalizationFactor() const { return ImageWeight(); }
  double BeamSize() const { return theoretical_beam_size_; }

  bool HasMetaDataCache() const {
    return (!metadata_cache_->msDataVector.empty());
  }
  void AllocateMetaDataCache(size_t size) {
    metadata_cache_->msDataVector.resize(size);
  }
  MetaDataCache::Entry& GetMetaDataCacheItem(size_t index) {
    return metadata_cache_->msDataVector[index];
  }
  void SetMetaDataCache(std::unique_ptr<MetaDataCache> cache) {
    metadata_cache_ = std::move(cache);
  }
  std::unique_ptr<MetaDataCache> AcquireMetaDataCache() {
    return std::move(metadata_cache_);
  }

  /**
   * Calculate overall metadata that is common to all types of gridders.
   */
  void CalculateOverallMetaData();
  /**
   * Calculate additional gridder specific metadata that was not calculated in
   * @ref CalculateOverallMetaData().
   */
  virtual void CalculateGridderMetaData(){};

  void SetMaxW(double max_w) { max_w_ = max_w; }
  void SetMaxBaseline(double max_baseline) { max_baseline_ = max_baseline; }
  void SetMinW(double min_w) { min_w_ = min_w; }

  // Prevent multiple calls of grid/invert on the same gridder at the same time.
  // Semaphore instead of mutex because we need to block in the scheduler but
  // release in the worker. See @ref
  // MSGridderManager::ExecuteForAllGriddersWithNCores() for more information.
  std::binary_semaphore processing_semaphore_{1};

 protected:
  size_t ActualInversionWidth() const { return actual_inversion_width_; }
  size_t ActualInversionHeight() const { return actual_inversion_height_; }
  double ActualPixelSizeX() const { return actual_pixel_size_x_; }
  double ActualPixelSizeY() const { return actual_pixel_size_y_; }

  size_t GetMsCount() const { return ms_data_vector_.size(); }
  MsProviderCollection::MsData& GetMsData(size_t ms_index) {
    return ms_data_vector_[ms_index];
  }

  virtual size_t GetSuggestedWGridSize() const = 0;

  /**
   * The largest w value present in the data (after applying any selections),
   * in units of number of wavelengths. It is initialized by
   * InitializeMSDataVector() and is undefined beforehand.
   */
  double MaxW() const { return max_w_; }
  double MaxBaseline() const { return max_baseline_; }
  double MinW() const { return min_w_; }

 private:
  bool hasWGridSize() const { return w_grid_size_ != 0; }

  // Contains all MS metadata as well as MS specific gridding data
  std::vector<MsProviderCollection::MsData>& ms_data_vector_;

  std::unique_ptr<MetaDataCache> metadata_cache_;
  size_t actual_inversion_width_ = 0;
  size_t actual_inversion_height_ = 0;
  double actual_pixel_size_x_ = 0.0;
  double actual_pixel_size_y_ = 0.0;
  double phase_centre_ra_ = 0.0;
  double phase_centre_dec_ = 0.0;
  size_t facet_index_ = 0;
  double image_padding_ = 1.0;
  size_t image_width_ = 0.0;
  size_t image_height_ = 0.0;
  size_t trim_width_ = 0;
  size_t trim_height_ = 0;
  size_t w_grid_size_ = 0;
  size_t actual_w_grid_size_ = 0;
  std::string data_column_name_;
  bool small_inversion_ = true;
  double max_w_ = 0.0;
  double min_w_ = 0.0;
  double max_baseline_ = 0.0;
  /// A fractional value that, when non-zero, places a limit on the w-value of
  /// gridded visibilities. Visibilities outside the limit are skipped.
  double w_limit_ = 0.0;
  bool is_complex_ = false;
  WeightMode weighting_ = WeightMode(WeightClass::Uniform);
  bool is_first_task_ = false;
  GriddingKernelMode grid_mode_ = GriddingKernelMode::KaiserBessel;
  double theoretical_beam_size_ = 0.0;
};

}  // namespace wsclean

#endif  // WSCLEAN_GRIDDING_MS_GRIDDER_H_
