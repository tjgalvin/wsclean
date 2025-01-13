#ifndef WSCLEAN_GRIDDING_MS_GRIDDER_DATA_H_
#define WSCLEAN_GRIDDING_MS_GRIDDER_DATA_H_

#include "gridmode.h"

#include <aocommon/banddata.h>
#include <aocommon/image.h>
#include <aocommon/polarization.h>
#include <aocommon/imagecoordinates.h>

#include "../msproviders/msreaders/msreader.h"

#include "h5solutiondata.h"
#include "msprovidercollection.h"
#include "visibilitymodifier.h"
#include "visibilityweightingmode.h"

#include "../main/settings.h"

#include "../scheduling/griddingtaskmanager.h"

#include <aocommon/uvector.h>

namespace wsclean {

enum class PsfMode {
  kNone,    // Not a psf, grid the visibilities in the MS
  kSingle,  // Grid generated visibilities for a point source at the centre of
            // the main image
  kDirectionDependent  // Grid generated visibilities for a point source at the
                       // centre of the current facet
};

namespace internal {

template <size_t PolarizationCount>
inline void CollapseData(
    size_t n_channels, std::complex<float>* buffer,
    [[maybe_unused]] aocommon::PolarizationEnum polarization) {
  if constexpr (PolarizationCount == 2) {
    for (size_t ch = 0; ch != n_channels; ++ch) {
      buffer[ch] = buffer[ch * PolarizationCount] +
                   buffer[(ch * PolarizationCount + (PolarizationCount - 1))];
    }
  } else if constexpr (PolarizationCount == 4) {
    for (size_t ch = 0; ch != n_channels; ++ch) {
      buffer[ch] = aocommon::Polarization::ConvertFromLinear(
          buffer + ch * PolarizationCount, polarization);
    }
  } else
    throw std::runtime_error("Invalid polarization conversion");
}

template <size_t PolarizationCount>
inline void ExpandData(
    size_t n_channels, std::complex<float>* buffer, std::complex<float>* output,
    [[maybe_unused]] aocommon::PolarizationEnum polarization) {
  if constexpr (PolarizationCount == 2) {
    for (size_t ch = 0; ch != n_channels; ++ch) {
      output[ch * 2] = buffer[ch];
      output[ch * 2 + 1] = buffer[ch];
    }
  } else if constexpr (PolarizationCount == 4) {
    for (size_t i = 0; i != n_channels; ++i) {
      const size_t ch = n_channels - 1 - i;
      aocommon::Polarization::ConvertToLinear(buffer[ch], polarization,
                                              &output[ch * 4]);
    }
  } else
    throw std::runtime_error("Invalid polarization conversion");
}

}  // namespace internal

/**
 * MsGridder stores a collection of methods and data required for computing and
 * working with the data required for gridding a measurement set, including
 * Calculation of weights and applying corrections.
 *
 * Data stored in this class is generally specific to a single measurement set,
 * however some more general (used across multiple MS) data is also stored for
 * convenience where it is required by one of the methods supported.
 *
 * MsGridderData is usually used via @ref MsGridder. However when doing shared
 * reads, @ref MsGridderManager instead instantiates a single @ref MsGridderData
 * for multiple gridders, uses it to prepare a single data buffer for gridding
 * and then passes that buffer into the gridders.
 */
class MsGridderData {
 public:
  MsGridderData(const Settings& settings);
  virtual ~MsGridderData() = default;

  /* Copy all member variables relating to task data into another MsGridderData.
   * member variables not related to task data are left as is.
   */
  void CopyTaskData(MsGridderData& other, const H5SolutionData& solution_data,
                    MsProviderCollection::MsData& ms_data) {
    psf_mode_ = other.psf_mode_;
    main_image_dl_ = other.main_image_dl_;
    main_image_dm_ = other.main_image_dm_;
    facet_group_index_ = other.facet_group_index_;
    is_facet_ = other.is_facet_;
    do_subtract_model_ = other.do_subtract_model_;
    polarization_ = other.polarization_;
    l_shift_ = other.l_shift_;
    m_shift_ = other.m_shift_;
    SetImageWeights(other.GetImageWeights());

    if (solution_data.HasData()) {
      visibility_modifier_.SetMSTimes(ms_data.original_ms_index,
                                      ms_data.unique_times);
      visibility_modifier_.SetH5Parm(
          solution_data.GetH5Parms(), solution_data.GetFirstSolutions(),
          solution_data.GetSecondSolutions(), solution_data.GetGainTypes());
    }
  }
  bool WillApplyCorrections() const {
    if (IsFacet() && (GetPsfMode() != PsfMode::kSingle)) {
#ifdef HAVE_EVERYBEAM
      const bool apply_beam =
          settings_.applyFacetBeam || settings_.gridWithBeam;
      if (apply_beam) return true;
#endif
      if (visibility_modifier_.HasH5Parm()) {
        return true;
      }
    }
    return false;
  }
  /**
   * Applies the selected visibility modifier (selected by Mode)
   * solutions to the visibilities and computes the weight corresponding to the
   * combined effect.
   *
   * @tparam Behaviour See @ref ModifierBehaviour and @ref @ref
   * ApplyConjugatedParmResponse for more information
   * @tparam LoadResponse This should always be true unless the calling code
   * knows the response has already been loaded previously, e.g. if we first
   * call `ApplyCorrections<Mode, kSum, true>(...)` we can then afterwards call
   * `ApplyCorrection<Mode, kApply, false>(...)` for the same values
   */
  template <GainMode Mode,
            ModifierBehaviour Behaviour = ModifierBehaviour::kApplyAndSum,
            bool LoadResponse = true, bool UseBufferedOffsets = false>
  void ApplyCorrections(size_t n_antennas, std::complex<float>* visibility_row,
                        const aocommon::BandData& band,
                        const float* weight_buffer,
                        const MSProvider::MetaData& metadata);

  template <GainMode Mode,
            ModifierBehaviour Behaviour = ModifierBehaviour::kApplyAndSum,
            bool LoadResponse = true>
  void ApplyCorrections(size_t n_antennas, std::complex<float>* visibility_row,
                        const aocommon::BandData& band,
                        const float* weight_buffer, double time,
                        size_t field_id, size_t antenna1, size_t antenna2,
                        size_t& time_offset, float* scratch_image_weights);

  /**
   * Initializes MS related data members, i.e. the @c _telescope and the
   * @c _pointResponse data in case a beam is applied on the facets and
   * EveryBeam is available and the @c _predictReader data member in case
   * @c isPredict is true.
   */
  void StartMeasurementSet(size_t ms_count,
                           const MsProviderCollection::MsData& ms_data,
                           bool is_predict);

  /**
   * Write (modelled) visibilities to MS, provides an interface to
   * MSProvider::WriteModel(). Any active facet beam or solution corrections
   * are applied. Method is templated on the number of
   * polarizations (1, 2 or 4). The gain_mode can be used to
   * select an entry or entries from the gain matrix that should be used for the
   * correction.
   * @param buffer should on entry contain n_channels visibilities to be
   * written.
   */
  void WriteCollapsedVisibilities(MSProvider& ms_provider, size_t n_antennas,
                                  const aocommon::BandData& band,
                                  std::complex<float>* buffer,
                                  MSProvider::MetaData& metadata) {
    switch (n_vis_polarizations_) {
      case 1:
        WriteInstrumentalVisibilities(ms_provider, n_antennas, band, buffer,
                                      metadata);
        break;
      case 2:
        internal::ExpandData<2>(band.ChannelCount(), buffer,
                                scratch_model_data_.data(), Polarization());
        WriteInstrumentalVisibilities(ms_provider, n_antennas, band,
                                      scratch_model_data_.data(), metadata);
        break;
      case 4:
        internal::ExpandData<4>(band.ChannelCount(), buffer,
                                scratch_model_data_.data(), Polarization());
        WriteInstrumentalVisibilities(ms_provider, n_antennas, band,
                                      scratch_model_data_.data(), metadata);
        break;
    }
  }

  /**
   * Similar to @ref WriteCollapsedVisibilities(), but assumes the input are
   * instrumental visibilities.
   * @param buffer n_polarizations x n_channels entries, which are the
   * instrumental visibilities.
   */
  void WriteInstrumentalVisibilities(MSProvider& ms_provider, size_t n_antennas,
                                     const aocommon::BandData& band,
                                     std::complex<float>* buffer,
                                     MSProvider::MetaData& metadata);

  bool HasDenormalPhaseCentre() const {
    return l_shift_ != 0.0 || m_shift_ != 0.0;
  }

  double ImageWeight() const {
    return total_weight_ / GetNVisibilities(gain_mode_);
  }

  void ReadPredictMetaData(MSProvider::MetaData& metadata);

  void ResetVisibilityModifierCache(size_t ms_count);

  double PixelSizeX() const { return settings_.pixelScaleX; }
  double PixelSizeY() const { return settings_.pixelScaleY; }

  /**
   * This is the sum of the weights as given by the measurement set, before the
   * image weighting is applied.
   */
  double VisibilityWeightSum() const { return visibility_weight_sum_; }
  /**
   * The number of visibilities that were gridded.
   */
  size_t GriddedVisibilityCount() const { return gridded_visibility_count_; }
  /**
   * The maximum weight, after having applied the imaging weighting.
   */
  double MaxGriddedWeight() const { return max_gridded_weight_; }
  /**
   * The effective number of visibilities, taking into account imaging weighting
   * and visibility weighting. This number is relative to the "best" visibility:
   * if one visibility with a weight of 10 and 5 visibilities with
   * a weight of 4 were gridded, the effective number of visibilities is
   * (10 + 5 x 4) / 10 = 3
   */
  double EffectiveGriddedVisibilityCount() const {
    return ImageWeight() / MaxGriddedWeight();
  }

  void ResetVisibilityCounters() {
    gridded_visibility_count_ = 0;
    total_weight_ = 0.0;
    max_gridded_weight_ = 0.0;
    visibility_weight_sum_ = 0.0;
  }

  bool DoSubtractModel() const { return do_subtract_model_; }
  void SetDoSubtractModel(bool do_subtract_model) {
    do_subtract_model_ = do_subtract_model;
  }
  PsfMode GetPsfMode() const { return psf_mode_; }
  void SetPsfMode(PsfMode psf_mode) { psf_mode_ = psf_mode; }
  VisibilityWeightingMode GetVisibilityWeightingMode() const {
    return visibility_weighting_mode_;
  }
  double MainImageDL() const { return main_image_dl_; }
  void SetMainImageDL(const double main_image_dl) {
    main_image_dl_ = main_image_dl;
  }
  double MainImageDM() const { return main_image_dm_; }
  void SetMainImageDM(const double main_image_dm) {
    main_image_dm_ = main_image_dm;
  }
  aocommon::PolarizationEnum Polarization() const { return polarization_; }
  void SetPolarization(aocommon::PolarizationEnum polarization) {
    polarization_ = polarization;
  }

  VisibilityModifier& GetVisibilityModifier() { return visibility_modifier_; };
  const VisibilityModifier& GetVisibilityModifier() const {
    return visibility_modifier_;
  };

  /**
   * @brief In case of facet-based imaging, the model data in the @param
   * MSProvider is reset to zeros in every major cycle, and predicted data
   * should be add-assigned to the model data (_isFacet = true) rather
   * than overwriting it. For standard imaging (_isFacet = false), the model
   * data should be overwritten.
   */
  void SetIsFacet(bool is_facet) { is_facet_ = is_facet; }
  bool IsFacet() const { return is_facet_; }
  void SetLShift(const double l_shift) { l_shift_ = l_shift; }
  double LShift() const { return l_shift_; }
  void SetMShift(const double m_shift) { m_shift_ = m_shift; }
  double MShift() const { return m_shift_; }

  const ImageWeights* GetImageWeights() const {
    return precalculated_weight_info_;
  }
  void SetImageWeights(const ImageWeights* weights) {
    precalculated_weight_info_ = weights;
  }

  bool StoreImagingWeights() const { return store_imaging_weights_; }
  void SetStoreImagingWeights(bool store_imaging_weights) {
    store_imaging_weights_ = store_imaging_weights;
  }

  void SetWriterLockManager(GriddingTaskManager* writer_lock_manager) {
    writer_lock_manager_ = writer_lock_manager;
  }
  void SetFacetGroupIndex(size_t index) { facet_group_index_ = index; }

  const Settings& GetSettings() const { return settings_; }

  GainMode GetGainMode() const { return gain_mode_; }

  /**
   * The average squared Mueller correction of all applied corrections.
   * This is the weighted sum of squared Mueller matrices, divided by the sum of
   * weights. It is zero if no corrections are applied.
   * @sa VisibilityModifier::TotalCorrectionSum().
   */
  AverageCorrection GetAverageCorrection() const {
    if (ImageWeight() != 0.0) {
      return visibility_modifier_.TotalCorrectionSum() / ImageWeight();
    } else {
      return AverageCorrection();
    }
  }
  /**
   * The average squared Mueller correction. This is the weighted sum of
   * squared Mueller matrices, divided by the sum of weights.
   * It is zero if not both beam and solution corrections are applied.
   * @sa VisibilityModifier::TotalCorrectionSum().
   */
  AverageCorrection GetAverageBeamCorrection() const {
    if (ImageWeight() != 0.0) {
      return visibility_modifier_.BeamCorrectionSum() / ImageWeight();
    } else {
      return AverageCorrection();
    }
  }

  struct InversionRow {
    double uvw[3];
    std::complex<float>* data;
  };

  template <size_t PolarizationCount>
  static void RotateVisibilities(const aocommon::BandData& band,
                                 double shift_factor,
                                 std::complex<float>* data_iter);

 protected:
  /**
   * Read a row of visibility and weights from the msprovider
   *
   * Use this function to correctly populate an InversionRow structure and an
   * accompanying weight_buffer and model_buffer before calling @ref
   * CollapseVisibilities() or @ref ApplyWeightsAndCorrections()
   *
   * @param ms_reader The measurement set provider from which data will be read
   * @param row_data The caller must set this object up to point at the desired
   * portion of an allocated buffer into which the visibilities will be read.
   * After returning from this call the uvw paramater of this object will be
   * populated with the (u/v/w)InM values of `metadata`
   * @param weight_buffer An allocated buffer of size n_chan x n_pol to store
   * intermediate weights in. After returning from the call, these values will
   * hold the weights from `ms_reader`
   * @param model_buffer An allocated buffer of size n_chan x n_pol to store
   * intermediate model data in.
   */
  inline void ReadVisibilities(MSReader& ms_reader,
                               std::complex<float>* row_data,
                               float* weight_buffer,
                               std::complex<float>* model_buffer) {
    if (GetPsfMode() == PsfMode::kNone) {
      ms_reader.ReadData(row_data);
    }
    if (DoSubtractModel()) {
      ms_reader.ReadModel(model_buffer);
    }
    ms_reader.ReadWeights(weight_buffer);
  }

  /**
   * Read a row of visibilities from the msprovider, and apply weights, flags
   * and a-terms.
   *
   * This function applies both the selected method of visibility weighting
   * (i.e. the weights that are normally stored in the WEIGHT_SPECTRUM column)
   * and the imaging weight (coming from uniform or Briggs weighting, etc).
   *
   * To read the data, this function requires scratch weight and model buffers
   * for storing intermediate values. Even if the caller does not need these
   * values, they still need to provide an already allocated buffer. This is to
   * avoid having to allocate memory within this method.
   *
   * This function collapses the visibilities in the polarization direction.
   * Gridders that grid a single polarization should use this method instead of
   * @ref GetInstrumentalVisibilities(). The output is stored in the first
   * n_channel elements of the visibility data buffer in @c row_data.
   *
   * @param ms_reader The measurement set provider from which data will be read
   * @param n_antennas The number of antennas
   * @param row_data The resulting weighted data
   * @param band The spectral band currently being imaged
   * @param weight_buffer An allocated buffer of size n_chan x n_pol to store
   * intermediate weights in. After returning from the call, these values will
   * hold the full applied weight (i.e. visibility weight * imaging weight).
   * @param model_buffer An allocated buffer of size n_chan x n_pol to store
   * intermediate model data in.
   * @param is_selected Per visibility whether that visibility will be gridded
   * in this pass. When the visibility is not gridded, its weight will not be
   * added to the relevant sums (visibility count, weight sum, etc.). This
   * buffer is of size n_chan; i.e. it is not specified per polarization.
   * @param metadata Metadata that has previously been read from a measurement
   * set provider
   */
  inline void GetCollapsedVisibilities(MSReader& ms_reader, size_t n_antennas,
                                       InversionRow& row_data,
                                       const aocommon::BandData& band,
                                       float* weight_buffer,
                                       std::complex<float>* model_buffer,
                                       const bool* is_selected,
                                       const MSProvider::MetaData& metadata) {
    ReadVisibilities(ms_reader, row_data.data, weight_buffer, model_buffer);

    CollapseVisibilities(n_antennas, row_data, band, weight_buffer,
                         model_buffer, is_selected, metadata);

    if (StoreImagingWeights())
      ms_reader.WriteImagingWeights(scratch_image_weights_.data());
  }

  /**
   * Same as @ref GetCollapsedVisibilities(), but without collapsing the
   * polarization direction. This implies that the output visibility buffer in
   * the row_data structure will contain n_channel x n_polarization elements.
   */
  template <size_t PolarizationCount>
  inline void GetInstrumentalVisibilities(
      MSReader& ms_reader, size_t n_antennas, InversionRow& row_data,
      const aocommon::BandData& band, float* weight_buffer,
      std::complex<float>* model_buffer, const bool* is_selected,
      const MSProvider::MetaData& metadata) {
    ReadVisibilities(ms_reader, row_data.data, weight_buffer, model_buffer);

    CalculateWeightsImplementation<PolarizationCount>(
        row_data.uvw, row_data.data, band, weight_buffer, model_buffer,
        is_selected);

    ApplyWeightsAndCorrections(n_antennas, row_data, band, weight_buffer,
                               metadata);

    if (StoreImagingWeights())
      ms_reader.WriteImagingWeights(scratch_image_weights_.data());
  }

  /**
   * @brief Apply corrections as well as visibility and imaging weights
   * Also computes the weight corresponding to the
   * combined effect of the corrections.
   *
   * Requires `scratch_image_weights_` to be populated which is usually done by
   * calling @ref CalculateWeights()
   */
  void ApplyWeightsAndCorrections(size_t n_antennas, InversionRow& row_data,
                                  const aocommon::BandData& band,
                                  float* weight_buffer,
                                  const MSProvider::MetaData& metadata);

  /**
   * Apply weights, flags and a-terms to a row of visibility data that has been
   * read by @ref ReadVisibilities() and collapse in the polarization direction
   *
   * This function applies both the selected method of visibility weighting
   * (i.e. the weights that are normally stored in the WEIGHT_SPECTRUM column)
   * and the imaging weight (coming from uniform or Briggs weighting, etc).
   *
   * To read the data, this function requires scratch weight and model buffers
   * for storing intermediate values. Even if the caller does not need these
   * values, they still need to provide an already allocated buffer. This is to
   * avoid having to allocate memory within this method.
   *
   * This function collapses the visibilities in the polarization direction.
   * Gridders that grid a single polarization should use this method instead of
   * @ref ApplyWeightsAndCorrections(). The output is stored in the first
   * n_channel elements of the visibility data buffer in @c row_data.
   *
   * Normally set to one when imaging a single
   * polarization. It may be set to 2 or 4 for IDG as it images multiple
   * polarizations at once, and it may be set to 2 or 4 when applying
   * solutions.
   * @param row_data The resulting weighted data
   * @param band The spectral band currently being imaged
   * @param weight_buffer An allocated buffer of size n_chan x n_pol to store
   * intermediate weights in. After returning from the call, these values will
   * hold the full applied weight (i.e. visibility weight * imaging weight).
   * @param model_buffer An allocated buffer of size n_chan x n_pol to store
   * intermediate model data in.
   * @param is_selected Per visibility whether that visibility will be gridded
   * in this pass. When the visibility is not gridded, its weight will not be
   * added to the relevant sums (visibility count, weight sum, etc.). This
   * buffer is of size n_chan; i.e. it is not specified per polarization.
   * @param metadata Metadata that has previously been read from a measurement
   * set provider
   */
  inline void CollapseVisibilities(size_t n_antennas, InversionRow& row_data,
                                   const aocommon::BandData& band,
                                   float* weight_buffer,
                                   std::complex<float>* model_buffer,
                                   const bool* is_selected,
                                   const MSProvider::MetaData& metadata) {
    switch (n_vis_polarizations_) {
      case 1:
        CalculateWeightsImplementation<1>(row_data.uvw, row_data.data, band,
                                          weight_buffer, model_buffer,
                                          is_selected);
        ApplyWeightsAndCorrections(n_antennas, row_data, band, weight_buffer,
                                   metadata);
        break;
      case 2:
        CalculateWeightsImplementation<2>(row_data.uvw, row_data.data, band,
                                          weight_buffer, model_buffer,
                                          is_selected);
        ApplyWeightsAndCorrections(n_antennas, row_data, band, weight_buffer,
                                   metadata);
        internal::CollapseData<2>(band.ChannelCount(), row_data.data,
                                  Polarization());
        break;
      case 4:
        CalculateWeightsImplementation<4>(row_data.uvw, row_data.data, band,
                                          weight_buffer, model_buffer,
                                          is_selected);
        ApplyWeightsAndCorrections(n_antennas, row_data, band, weight_buffer,
                                   metadata);
        internal::CollapseData<4>(band.ChannelCount(), row_data.data,
                                  Polarization());
        break;
    }
  }

 private:
  /**
   * @brief Apply visibility and imaging weights
   * Requires `scratch_image_weights_` to be populated which is usually done by
   * calling @ref CalculateWeights()
   */
  template <GainMode Mode>
  void ApplyWeights(std::complex<float>* visibility_row,
                    const size_t channel_count, float* weight_buffer);

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
  void ApplyConjugatedH5Parm(MSReader& ms_reader,
                             const std::vector<std::string>& antenna_names,
                             InversionRow& row_data,
                             const aocommon::BandData& band,
                             const float* weight_buffer,
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
  void ApplyConjugatedFacetBeam(MSReader& ms_reader, InversionRow& row_data,
                                const aocommon::BandData& band,
                                const float* weight_buffer,
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
      MSReader& ms_reader, const std::vector<std::string>& antenna_names,
      InversionRow& row_data, const aocommon::BandData& band,
      const float* weight_buffer, bool apply_forward = false);
#endif  // HAVE_EVERYBEAM

  inline void CalculateWeights(double* uvw_buffer,
                               std::complex<float>* visibility_buffer,
                               const aocommon::BandData& band,
                               float* weight_buffer,
                               std::complex<float>* model_buffer,
                               const bool* is_selected);

  template <size_t PolarizationCount>
  void CalculateWeightsImplementation(double* uvw_buffer,
                                      std::complex<float>* visibility_buffer,
                                      const aocommon::BandData& band,
                                      float* weight_buffer,
                                      std::complex<float>* model_buffer,
                                      const bool* is_selected);

  void InitializePointResponse(const MsProviderCollection::MsData& ms_data);

  template <GainMode Mode>
  void WriteInstrumentalVisibilities(MSProvider& ms_provider, size_t n_antennas,
                                     const aocommon::BandData& band,
                                     std::complex<float>* buffer,
                                     MSProvider::MetaData& metadata);

  const Settings& settings_;

  VisibilityWeightingMode visibility_weighting_mode_ =
      VisibilityWeightingMode::NormalVisibilityWeighting;

  // Reset by the gridders at the start of each inversion, incremented during
  // gridding
  size_t gridded_visibility_count_ = 0;
  double total_weight_ = 0.0;
  double max_gridded_weight_ = 0.0;
  double visibility_weight_sum_ = 0.0;

  // These members are set from the task during InitializeGridderForTask
  bool do_subtract_model_ = false;
  bool is_facet_ = false;  /// @see SetIsFacet()
  bool store_imaging_weights_ = false;
  double main_image_dl_ = 0.0;
  double main_image_dm_ = 0.0;
  double l_shift_ = 0.0;
  double m_shift_ = 0.0;
  aocommon::PolarizationEnum polarization_ = aocommon::Polarization::StokesI;
  const ImageWeights* precalculated_weight_info_ = nullptr;
  PsfMode psf_mode_ = PsfMode::kNone;

  /// These members are initialised from MsData during StartMeasurementSet
  size_t n_vis_polarizations_ = 1;
  size_t original_ms_index_ = 0;
  VisibilityModifier visibility_modifier_;
  GainMode gain_mode_ = GainMode::kTrace;
  /// Used in WriteCollapsedVisibilities() to expand visibilities into.
  aocommon::UVector<std::complex<float>> scratch_model_data_;
  aocommon::UVector<float> scratch_image_weights_;

  /* @p _facetGroupIndex and @p _msIndex in conjunction with the @p GetMsCount()
   * determine the index in the _writerGroupLocks vector, having size
   * FacetGroupCount() * GetMsCount().
   * These variable are only relevant for prediction.
   */
  size_t facet_group_index_ = 0;
  GriddingTaskManager* writer_lock_manager_ = nullptr;
  size_t writer_lock_index_ = 0;
  std::unique_ptr<MSReader> predict_reader_;

  /* @ref MsGridderManager needs to access various methods that we would
   * otherwise need to make public, as only certain parts of the code should
   * access these methods we make @ref MsGridderManager a friend instead.
   */
  friend class MSGridderManager;
};

template <GainMode Mode, ModifierBehaviour Behaviour, bool LoadResponse,
          bool UseBufferedOffsets>
void MsGridderData::ApplyCorrections(size_t n_antennas,
                                     std::complex<float>* visibility_row,
                                     const aocommon::BandData& band,
                                     const float* weight_buffer,
                                     const MSProvider::MetaData& metadata) {
  size_t time_offset = visibility_modifier_.GetTimeOffset(original_ms_index_);
  ApplyCorrections<Mode, Behaviour, LoadResponse>(
      n_antennas, visibility_row, band, weight_buffer, metadata.time,
      metadata.fieldId, metadata.antenna1, metadata.antenna2, time_offset,
      scratch_image_weights_.data());
  visibility_modifier_.SetTimeOffset(original_ms_index_, time_offset);
}

// We can safely pass nullptr for weight buffer and image weights as well as
// 0 for time and field_id because these are unused in
// ModifierBehaviour::kApply mode
template <GainMode Mode, ModifierBehaviour Behaviour, bool LoadResponse>
void MsGridderData::ApplyCorrections(size_t n_antennas,
                                     std::complex<float>* visibility_row,
                                     const aocommon::BandData& band,
                                     const float* weight_buffer, double time,
                                     size_t field_id, size_t antenna1,
                                     size_t antenna2, size_t& time_offset,
                                     float* scratch_image_weights) {
  assert((weight_buffer == nullptr) ==
         (Behaviour == ModifierBehaviour::kApply));

  if (IsFacet() && (GetPsfMode() != PsfMode::kSingle)) {
    const bool apply_beam = settings_.applyFacetBeam || settings_.gridWithBeam;
    const bool apply_forward = GetPsfMode() == PsfMode::kDirectionDependent;
    if (apply_beam && visibility_modifier_.HasH5Parm()) {
#ifdef HAVE_EVERYBEAM
      // Load and apply (in conjugate) both the beam and the h5parm solutions
      if constexpr (LoadResponse) {
        visibility_modifier_.CacheBeamResponse(time, field_id, band);
        visibility_modifier_.CacheParmResponse(time, band, original_ms_index_,
                                               time_offset);
      }
      visibility_modifier_.ApplyConjugatedDual<Behaviour, Mode>(
          visibility_row, weight_buffer, scratch_image_weights,
          band.ChannelCount(), n_antennas, antenna1, antenna2,
          original_ms_index_, apply_forward, time_offset);
    } else if (apply_beam) {
      // Load and apply only the conjugate beam
      if constexpr (LoadResponse) {
        visibility_modifier_.CacheBeamResponse(time, field_id, band);
      }
      visibility_modifier_.ApplyConjugatedBeamResponse<Behaviour, Mode>(
          visibility_row, weight_buffer, scratch_image_weights,
          band.ChannelCount(), antenna1, antenna2, apply_forward);

#endif  // HAVE_EVERYBEAM
    } else if (visibility_modifier_.HasH5Parm()) {
      // Load and apply the h5parm solutions
      if constexpr (LoadResponse) {
        visibility_modifier_.CacheParmResponse(time, band, original_ms_index_,
                                               time_offset);
      }
      visibility_modifier_.ApplyConjugatedParmResponse<Behaviour, Mode>(
          visibility_row, weight_buffer, scratch_image_weights,
          original_ms_index_, band.ChannelCount(), n_antennas, antenna1,
          antenna2, apply_forward, time_offset);
    }
  }
}

template <GainMode Mode>
inline void MsGridderData::ApplyWeights(std::complex<float>* visibility_row,
                                        const size_t channel_count,
                                        float* weight_buffer) {
  const size_t n_pols = GetNVisibilities(Mode);

  for (size_t channel = 0; channel < channel_count; channel++) {
    for (size_t pol = 0; pol < n_pols; pol++) {
      size_t i = channel * n_pols + pol;

      const float cumWeight =
          weight_buffer[i] * scratch_image_weights_[channel];
      // We can use the boolean for computation instead of an if-condition
      // within the loop. This allows the inner part of the loop to be
      // autovectorized more easily.
      const bool has_weight = cumWeight != 0.0;
      if (pol == 0) {
        // Visibility weight sum is the sum of weights excluding imaging weights
        visibility_weight_sum_ += weight_buffer[i] * has_weight;
        max_gridded_weight_ =
            std::max(static_cast<double>(cumWeight), max_gridded_weight_);
        gridded_visibility_count_ += has_weight;
      }
      // Total weight includes imaging weights
      total_weight_ += cumWeight;
      weight_buffer[i] = cumWeight;
      visibility_row[i] *= cumWeight;
    }
  }
}

inline void MsGridderData::CalculateWeights(
    double* uvw_buffer, std::complex<float>* visibility_buffer,
    const aocommon::BandData& band, float* weight_buffer,
    std::complex<float>* model_buffer, const bool* is_selected) {
  switch (n_vis_polarizations_) {
    case 1:
      CalculateWeightsImplementation<1>(uvw_buffer, visibility_buffer, band,
                                        weight_buffer, model_buffer,
                                        is_selected);
      break;
    case 2:
      CalculateWeightsImplementation<2>(uvw_buffer, visibility_buffer, band,
                                        weight_buffer, model_buffer,
                                        is_selected);
      break;
    case 4:
      CalculateWeightsImplementation<4>(uvw_buffer, visibility_buffer, band,
                                        weight_buffer, model_buffer,
                                        is_selected);
      break;
  }
}

template <size_t PolarizationCount>
void MsGridderData::CalculateWeightsImplementation(
    double* uvw_buffer, std::complex<float>* visibility_buffer,
    const aocommon::BandData& band, float* weight_buffer,
    std::complex<float>* model_buffer, const bool* is_selected) {
  const std::size_t data_size = band.ChannelCount() * PolarizationCount;
  if (GetPsfMode() != PsfMode::kNone) {
    // Visibilities for a point source at the phase centre are all ones
    std::fill_n(visibility_buffer, data_size, 1.0);
    double dl = 0.0;
    double dm = 0.0;
    if (GetPsfMode() == PsfMode::kSingle) {
      // The point source is shifted to the centre of the main image
      dl = MainImageDL();
      dm = MainImageDM();
    } else {  // GetPsfMode() == PsfMode::kDirectionDependent
      // The point source is shifted to the centre of the current DdPsf
      // position
      dl = LShift();
      dm = MShift();
    }
    if (dl != 0.0 || dm != 0.0) {
      const double dn = std::sqrt(1.0 - dl * dl - dm * dm) - 1.0;
      const double shift_factor =
          2.0 * M_PI *
          (uvw_buffer[0] * dl + uvw_buffer[1] * dm + uvw_buffer[2] * dn);
      RotateVisibilities<PolarizationCount>(band, shift_factor,
                                            visibility_buffer);
    }
  }

  if (DoSubtractModel()) {
    std::complex<float>* model_iter = model_buffer;
    for (std::complex<float>* iter = visibility_buffer;
         iter != visibility_buffer + data_size; ++iter) {
      *iter -= *model_iter;
      model_iter++;
    }
  }

  // Any visibilities that are not gridded in this pass
  // should not contribute to the weight sum, so set these
  // to have zero weight.
  for (size_t ch = 0; ch != band.ChannelCount(); ++ch) {
    for (size_t p = 0; p != PolarizationCount; ++p) {
      if (!is_selected[ch]) weight_buffer[ch * PolarizationCount + p] = 0.0;
    }
  }

  switch (GetVisibilityWeightingMode()) {
    case VisibilityWeightingMode::NormalVisibilityWeighting:
      // The weight buffer already contains the visibility weights: do nothing
      break;
    case VisibilityWeightingMode::SquaredVisibilityWeighting:
      // Square the visibility weights
      for (size_t i = 0; i != data_size; ++i)
        weight_buffer[i] *= weight_buffer[i];
      break;
    case VisibilityWeightingMode::UnitVisibilityWeighting:
      // Set the visibility weights to one
      for (size_t i = 0; i != data_size; ++i) {
        if (weight_buffer[i] != 0.0) weight_buffer[i] = 1.0f;
      }
      break;
  }

  // Precompute imaging weights
  for (size_t ch = 0; ch != band.ChannelCount(); ++ch) {
    const double u = uvw_buffer[0] / band.ChannelWavelength(ch);
    const double v = uvw_buffer[1] / band.ChannelWavelength(ch);
    scratch_image_weights_[ch] = GetImageWeights()->GetWeight(u, v);
  }
}

// Apply corrections as well as visibility and imaging weights
inline void MsGridderData::ApplyWeightsAndCorrections(
    size_t n_antennas, InversionRow& row_data, const aocommon::BandData& band,
    float* weight_buffer, const MSProvider::MetaData& metadata) {
  switch (gain_mode_) {
    case GainMode::kXX:
      ApplyCorrections<GainMode::kXX>(n_antennas, row_data.data, band,
                                      weight_buffer, metadata);
      ApplyWeights<GainMode::kXX>(row_data.data, band.ChannelCount(),
                                  weight_buffer);
      break;
    case GainMode::kYY:
      ApplyCorrections<GainMode::kYY>(n_antennas, row_data.data, band,
                                      weight_buffer, metadata);
      ApplyWeights<GainMode::kYY>(row_data.data, band.ChannelCount(),
                                  weight_buffer);
      break;
    case GainMode::kTrace:
      ApplyCorrections<GainMode::kTrace>(n_antennas, row_data.data, band,
                                         weight_buffer, metadata);
      ApplyWeights<GainMode::kTrace>(row_data.data, band.ChannelCount(),
                                     weight_buffer);
      break;
    case GainMode::k2VisDiagonal:
      ApplyCorrections<GainMode::k2VisDiagonal>(n_antennas, row_data.data, band,
                                                weight_buffer, metadata);
      ApplyWeights<GainMode::k2VisDiagonal>(row_data.data, band.ChannelCount(),
                                            weight_buffer);
      break;
    case GainMode::kFull:
      ApplyCorrections<GainMode::kFull>(n_antennas, row_data.data, band,
                                        weight_buffer, metadata);
      ApplyWeights<GainMode::kFull>(row_data.data, band.ChannelCount(),
                                    weight_buffer);
      break;
    default:
      throw std::runtime_error(
          "Invalid combination of visibility polarizations and gain mode");
  }
}

template <size_t PolarizationCount>
void MsGridderData::RotateVisibilities(const aocommon::BandData& band,
                                       double shift_factor,
                                       std::complex<float>* data_iter) {
  for (size_t ch = 0; ch != band.ChannelCount(); ++ch) {
    const double wShiftRad = shift_factor / band.ChannelWavelength(ch);
    const std::complex<float> phasor(std::cos(wShiftRad), std::sin(wShiftRad));
    for (size_t p = 0; p != PolarizationCount; ++p) {
      *data_iter *= phasor;
      ++data_iter;
    }
  }
}

inline void MsGridderData::WriteInstrumentalVisibilities(
    MSProvider& ms_provider, size_t n_antennas, const aocommon::BandData& band,
    std::complex<float>* buffer, MSProvider::MetaData& metadata) {
  switch (gain_mode_) {
    case GainMode::kXX:
      WriteInstrumentalVisibilities<GainMode::kXX>(ms_provider, n_antennas,
                                                   band, buffer, metadata);
      break;
    case GainMode::kYY:
      WriteInstrumentalVisibilities<GainMode::kYY>(ms_provider, n_antennas,
                                                   band, buffer, metadata);
      break;
    case GainMode::kTrace:
      WriteInstrumentalVisibilities<GainMode::kTrace>(ms_provider, n_antennas,
                                                      band, buffer, metadata);
      break;
    case GainMode::k2VisDiagonal:
      WriteInstrumentalVisibilities<GainMode::k2VisDiagonal>(
          ms_provider, n_antennas, band, buffer, metadata);
      break;
    case GainMode::kFull:
      WriteInstrumentalVisibilities<GainMode::kFull>(ms_provider, n_antennas,
                                                     band, buffer, metadata);
      break;
  }
}

template <GainMode Mode>
void MsGridderData::WriteInstrumentalVisibilities(
    MSProvider& ms_provider, size_t n_antennas, const aocommon::BandData& band,
    std::complex<float>* buffer, MSProvider::MetaData& metadata) {
  assert(GetPsfMode() == PsfMode::kNone);  // The PSF is never predicted.

#ifdef HAVE_EVERYBEAM
  if (settings_.applyFacetBeam) {
    visibility_modifier_.CacheBeamResponse(metadata.time, metadata.fieldId,
                                           band);

    visibility_modifier_.ApplyBeamResponse<Mode>(
        buffer, band.ChannelCount(), metadata.antenna1, metadata.antenna2);
  }
#endif

  if (visibility_modifier_.HasH5Parm()) {
    assert(!settings_.facetRegionFilename.empty());
    size_t time_offset = visibility_modifier_.GetTimeOffset(original_ms_index_);
    visibility_modifier_.CacheParmResponse(metadata.time, band,
                                           original_ms_index_, time_offset);
    visibility_modifier_.ApplyParmResponse<Mode>(
        buffer, original_ms_index_, band.ChannelCount(), n_antennas,
        metadata.antenna1, metadata.antenna2, time_offset);
    visibility_modifier_.SetTimeOffset(original_ms_index_, time_offset);
  }

  {
    std::unique_ptr<GriddingTaskManager::WriterLock> lock =
        writer_lock_manager_->GetLock(writer_lock_index_);
    ms_provider.WriteModel(buffer, IsFacet());
  }
  ms_provider.NextOutputRow();
}

}  // namespace wsclean

#endif  // WSCLEAN_GRIDDING_MS_GRIDDER_DATA_H_
