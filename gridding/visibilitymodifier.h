#ifndef GRIDDING_VISIBILITY_MODIFIER_H_
#define GRIDDING_VISIBILITY_MODIFIER_H_

#include <stdexcept>

#ifdef HAVE_EVERYBEAM
#include <EveryBeam/beammode.h>
#include <EveryBeam/beamnormalisationmode.h>
#include <EveryBeam/pointresponse/pointresponse.h>
#endif

#include <schaapcommon/h5parm/jonesparameters.h>

#include <aocommon/constants.h>
#include <aocommon/matrix2x2.h>
#include <aocommon/multibanddata.h>
#include <aocommon/polarization.h>
#include <aocommon/vectormap.h>

#include "averagecorrection.h"
#include "gainmode.h"

namespace wsclean {

class SynchronizedMS;
struct ObservationInfo;

/**
 * Control whether to apply the modifier, sum the corrections or both when
 * calling methods that apply a modifier e.g. @ref ApplyConjugatedParmResponse()
 */
enum class ModifierBehaviour { kApply, kSum, kApplyAndSum };

using BeamResponseMap = aocommon::VectorMap<std::vector<aocommon::MC2x2F>>;
using SolutionsResponseMap =
    aocommon::VectorMap<std::vector<std::complex<float>>>;

namespace internal {
/**
 * @brief Apply gains to the visibilities.
 *
 * @tparam Mode Which entry or entries from the gain matrices should be
 * taken into account when correcting the visibilities? See also the
 * documentation of GainMode.
 */
template <GainMode Mode>
void ApplyGain(std::complex<float>* visibilities, const aocommon::MC2x2F& gain1,
               const aocommon::MC2x2F& gain2) {
  if constexpr (Mode == GainMode::kXX) {
    *visibilities = gain1.Get(0) * (*visibilities) * std::conj(gain2.Get(0));
  } else if constexpr (Mode == GainMode::kYY) {
    *visibilities = gain1.Get(3) * (*visibilities) * std::conj(gain2.Get(3));
  } else if constexpr (Mode == GainMode::kTrace) {
    // Stokes-I. Have to calculate v' = G1 x v x G2^H:
    // v' = 0.5 (V_xx + V_yy) with V = v x (G1 x G2^H)
    // V_xx = v x (g1_xx g2_xx* + g1_yx g2_yx*), V_yy = v x (g1_xy g2_xy* +
    // g1_yy g2_yy*). Hence v' = 0.5 * double_dot(G1, G2*)
    *visibilities *= 0.5f * gain1.DoubleDot(gain2.Conjugate());
  } else if constexpr (Mode == GainMode::k2VisDiagonal) {
    visibilities[0] = gain1.Get(0) * visibilities[0] * std::conj(gain2.Get(0)) +
                      gain1.Get(1) * visibilities[1] * std::conj(gain2.Get(1));
    visibilities[1] = gain1.Get(3) * visibilities[1] * std::conj(gain2.Get(3)) +
                      gain1.Get(2) * visibilities[0] * std::conj(gain2.Get(2));
  } else if constexpr (Mode == GainMode::kFull) {
    // All polarizations
    const aocommon::MC2x2F visibilities_mc2x2(visibilities);
    const aocommon::MC2x2F result =
        gain1.Multiply(visibilities_mc2x2).MultiplyHerm(gain2);
    result.AssignTo(visibilities);
  }
}

/**
 * @brief Apply conjugated gains to the visibilities.
 *
 * @tparam Mode Which entry or entries from the gain matrices should be
 * taken into account when correcting the visibilities? See also the
 * documentation of GainMode.
 */
template <GainMode Mode>
void ApplyConjugatedGain(std::complex<float>* visibilities,
                         const aocommon::MC2x2F& gain1,
                         const aocommon::MC2x2F& gain2) {
  if constexpr (Mode == GainMode::kXX) {
    *visibilities = std::conj(gain1.Get(0)) * (*visibilities) * gain2.Get(0);
  } else if constexpr (Mode == GainMode::kYY) {
    *visibilities = std::conj(gain1.Get(3)) * (*visibilities) * gain2.Get(3);
  } else if constexpr (Mode == GainMode::kTrace) {
    // See calculation in ApplyGain() for explanation of double dot.
    *visibilities *= 0.5f * gain2.DoubleDot(gain1.Conjugate());
  } else if constexpr (Mode == GainMode::k2VisDiagonal) {
    visibilities[0] = std::conj(gain1.Get(0)) * visibilities[0] * gain2.Get(0) +
                      std::conj(gain1.Get(2)) * visibilities[1] * gain2.Get(2);
    visibilities[1] = std::conj(gain1.Get(3)) * visibilities[1] * gain2.Get(3) +
                      std::conj(gain1.Get(1)) * visibilities[0] * gain2.Get(1);
  } else if constexpr (Mode == GainMode::kFull) {
    // All polarizations
    const aocommon::MC2x2F visibilities_mc2x2(visibilities);
    const aocommon::MC2x2F result =
        gain1.HermThenMultiply(visibilities_mc2x2).Multiply(gain2);
    result.AssignTo(visibilities);
  }
}
constexpr bool ShouldApplyCorrection(ModifierBehaviour behaviour) {
  return behaviour == ModifierBehaviour::kApply ||
         behaviour == ModifierBehaviour::kApplyAndSum;
}

constexpr bool ShouldSumCorrection(ModifierBehaviour behaviour) {
  return behaviour == ModifierBehaviour::kSum ||
         behaviour == ModifierBehaviour::kApplyAndSum;
}

template <GainMode Mode, typename T>
constexpr decltype(auto) MakeDiagonalIfScalar(T& matrix) {
  if constexpr (AllowScalarCorrection(Mode)) {
    return Diagonal(matrix);
  } else {
    return (matrix);
  }
}
template <size_t NParms, GainMode Mode, typename MatrixArrayType>
auto CreateMatrix2x2OrDiag(MatrixArrayType data, size_t offset) {
  if constexpr (NParms == 2) {
    return aocommon::MC2x2FDiag(data[offset], data[offset + 1]);
  } else {
    if constexpr (AllowScalarCorrection(Mode)) {
      return aocommon::MC2x2FDiag(data[offset], data[offset + 3]);
    } else {
      return aocommon::MC2x2F(&data[offset]);
    }
  }
}
template <size_t NParms, typename MatrixArrayType>
auto CreateMatrix2x2OrDiag(MatrixArrayType data, size_t offset) {
  if constexpr (NParms == 2) {
    return aocommon::MC2x2FDiag(data[offset], data[offset + 1]);
  } else {
    return aocommon::MC2x2F(&data[offset]);
  }
}

template <size_t NParms, typename MatrixArrayType>
auto CreateMatrix2x2(MatrixArrayType data, size_t offset) {
  if constexpr (NParms == 2) {
    return aocommon::MC2x2F(data[offset], 0.0f, 0.0f, data[offset + 1]);
  } else {
    return aocommon::MC2x2F(&data[offset]);
  }
}
}  // namespace internal

/**
 * Cache the beam response for all timesteps in the rows of a chunk, maintain a
 * row to timestep mapping so that the responses are indexable by row number.
 */
struct BeamResponseCacheChunk {
  std::vector<uint32_t> offsets;
  /**
   * The outer vector contains an element for every timestep in the chunk. The
   * elements themselves have a structure like the beam response VectorMaps in
   * VisibilityModifier, see e.g. @ref
   * VisibilityModifier::GetCachedBeamResponse().
   */
  std::vector<BeamResponseMap> responses;

  const BeamResponseMap& GetCachedBeamResponseForRow(size_t row) const {
    assert(row < offsets.size());
    assert(offsets[row] < responses.size());
    return responses[offsets[row]];
  }
};

/**
 * Applies beam and h5parm solutions to visibilities.
 * See the documentation for function @ref ApplyConjugatedParmResponse()
 * for an overview of parameters that hold for most of these functions.
 */
class VisibilityModifier {
 public:
  VisibilityModifier() = default;

  void InitializeBda(bool uses_bda) { uses_bda_ = uses_bda; }

  void InitializePointResponse(SynchronizedMS&& ms,
                               double facet_beam_update_time,
                               const std::string& element_response,
                               const aocommon::MultiBandData& bands,
                               const std::string& data_column,
                               const std::string& mwa_path);

  void InitializeTimeFrequencySmearing(const ObservationInfo& observation_info,
                                       double interval);

  /**
   * A function that initializes this visibility modifier for testing. After
   * calling this function, the Apply* functions can be called.
   * @param beam_response vector of n_channels * n_stations * 4 gain elements.
   * @param parm_response vector of n_channels * n_stations * 2
   */
  void InitializeMockResponse(
      size_t n_channels, size_t n_stations, size_t data_desc_id,
      const std::vector<aocommon::MC2x2F>& beam_response,
      const std::vector<std::complex<float>>& parm_response);

  void SetNoPointResponse() {
#ifdef HAVE_EVERYBEAM
    _pointResponse = nullptr;
    _cachedBeamResponse.Clear();
#endif
  }

  void SetBeamInfo(std::string mode, std::string normalisation) {
#ifdef HAVE_EVERYBEAM
    _beamModeString = std::move(mode);
    _beamNormalisationMode = std::move(normalisation);
#endif
  }

  void ResetCache(size_t n_measurement_sets) {
    _cachedParmResponse.clear();
    _cachedMSTimes.clear();
    time_offsets_.clear();
#ifdef HAVE_EVERYBEAM
    beam_cache_chunks_.clear();
    current_beam_cache_chunk_ = nullptr;
    previous_beam_cache_chunk_ = nullptr;
#endif
  }

  /**
   * @brief Cache the solutions from a h5 solution file
   */
  void InitializeCacheParmResponse(
      const std::vector<std::string>& antenna_names,
      const aocommon::MultiBandData& selected_bands, size_t ms_index);
  /**
   * @brief Calculate how much memory cached h5 solution files are
   * consuming after calling @ref InitializeCacheParmResponse()
   * @return The memory consumption, in bytes.
   */
  size_t GetCacheParmResponseSize() const;

  /**
   * Update the time associated with a cached h5 solution file.
   * Cache must be initialised once per MsData before calling this
   * method, by calling @ref InitializeCacheParmResponse()
   *
   * @param time_offset This value is incrementally updated by each call into
   * cacheParmResponse and epresents an offset into @ref _cachedParmResponse
   * that is needed by functions that will apply the parm response e.g. @ref
   * ApplyParmResponse(). Also used internally as an offset into @ref
   * _cachedMSTimes to avoid searching data already covered by previous calls.
   * See @ref GetTimeOffset() for more information.
   */
  void CacheParmResponse(double time, const aocommon::MultiBandData& bands,
                         size_t ms_index, size_t& time_offset);

  const SolutionsResponseMap& GetCachedParmResponse(size_t ms_index) const {
    return _cachedParmResponse.find(ms_index)->second;
  }

  /**
   * Applies the conjugate (is backward, or imaging direction) h5parm gain
   * to given data.
   * @tparam Behaviour Determines whether we should apply the gain, sum the
   * correction or both.
   * @tparam Mode Gain application mode that defines how the gain is
   * applied, the PolarizationCount is also implied/determined by the gain mode.
   * @param [in,out] data Data array with n_channels x PolarizationCount
   * elements.
   * @param weights Array with for each data value the corresponding weight.
   * @param image_weights Array of size n_channels (polarizations are assumed to
   * have equal imaging weights) with the imaging weighting mode (e.g. from
   * Briggs weighting).
   * @param apply_forward If true, an additional (non-conjugate) forward gain is
   * applied to the data. This is necessary for calculating direction-dependent
   * PSFs.
   */
  template <ModifierBehaviour Behaviour, GainMode Mode, size_t NParms,
            bool ApplyForward>
  void ApplyConjugatedParmResponseForRow(std::complex<float>* data,
                                         const float* weights,
                                         const float* image_weights,
                                         size_t ms_index, size_t n_channels,
                                         size_t data_desc_id, size_t n_antennas,
                                         size_t antenna1, size_t antenna2) {
    ApplyConjugatedParmResponseForRow<Behaviour, Mode, NParms, ApplyForward>(
        data, weights, image_weights, ms_index, n_channels, data_desc_id,
        n_antennas, antenna1, antenna2, GetTimeOffset(ms_index));
  }
  template <ModifierBehaviour Behaviour, GainMode Mode, size_t NParms,
            bool ApplyForward>
  void ApplyConjugatedParmResponseForRow(std::complex<float>* data,
                                         const float* weights,
                                         const float* image_weights,
                                         size_t ms_index, size_t n_channels,
                                         size_t data_desc_id, size_t n_antennas,
                                         size_t antenna1, size_t antenna2,
                                         size_t time_offset) {
    const std::complex<float>* parm_response =
        _cachedParmResponse[ms_index][data_desc_id].data();
    const size_t n_visibilities = GetNVisibilities(Mode);
    for (size_t n_channel = 0; n_channel < n_channels; ++n_channel) {
      ApplyConjugatedParmResponse<Behaviour, Mode, NParms, ApplyForward>(
          parm_response, data, weights, image_weights, n_channel, n_channels,
          n_antennas, antenna1, antenna2, time_offset);
      if constexpr (internal::ShouldApplyCorrection(Behaviour)) {
        data += n_visibilities;
      }
      if constexpr (internal::ShouldSumCorrection(Behaviour)) {
        weights += n_visibilities;
      }
    }
  }
  template <ModifierBehaviour Behaviour, GainMode Mode, size_t NParms,
            bool ApplyForward>
  void ApplyConjugatedParmResponse(const std::complex<float>* parm_response,
                                   std::complex<float>* data,
                                   const float* weights,
                                   const float* image_weights, size_t n_channel,
                                   size_t n_channels, size_t n_antennas,
                                   size_t antenna1, size_t antenna2,
                                   size_t time_offset);

  template <GainMode GainEntry>
  void ApplyParmResponse(std::complex<float>* data, size_t ms_index,
                         size_t n_channels, size_t data_desc_id,
                         size_t n_antennas, size_t antenna1, size_t antenna2) {
    ApplyParmResponseWithTime<GainEntry>(data, ms_index, n_channels,
                                         data_desc_id, n_antennas, antenna1,
                                         antenna2, GetTimeOffset(ms_index));
  }
  template <GainMode GainEntry>
  void ApplyParmResponseWithTime(std::complex<float>* data, size_t ms_index,
                                 size_t n_channels, size_t data_desc_id,
                                 size_t n_antennas, size_t antenna1,
                                 size_t antenna2, size_t time_offset);

  void ApplyTimeFrequencySmearing(std::complex<float>* data, const double* uvw,
                                  const aocommon::BandData& band, int n,
                                  double l, double m);

  void SetMSTimes(size_t ms_index, std::shared_ptr<std::vector<double>> times) {
    _cachedMSTimes[ms_index] = std::move(times);
  }

  void SetH5Parm(
      const std::vector<schaapcommon::h5parm::H5Parm>& h5parms,
      const std::vector<schaapcommon::h5parm::SolTab*>& first_solutions,
      const std::vector<schaapcommon::h5parm::SolTab*>& second_solutions,
      const std::vector<schaapcommon::h5parm::GainType>& gain_types) {
    _h5parms = &h5parms;
    _firstSolutions = &first_solutions;
    _secondSolutions = &second_solutions;
    _gainTypes = &gain_types;
  }

  bool HasH5Parm() const { return _h5parms && !_h5parms->empty(); }

  /**
   * Return the current time offset for the MS corresponding to `ms_index`
   * The time offset is a value that is incrementally updated inside @ref
   * CacheParmResponse() for every @ref ApplyCorrection() call and represents an
   * offset from the start of @ref _cachedMSTimes that can be used to speed up
   * subsequent searches in subsequent calls to @ref CacheParmResponse() as well
   * as an index into @ref _cachedParmResponse that is needed by functions that
   * apply the parm response e.g. @ref ApplyParmResponse()
   */
  size_t GetTimeOffset(size_t ms_index) { return time_offsets_[ms_index]; }
  void SetTimeOffset(size_t ms_index, size_t time_offset) {
    time_offsets_[ms_index] = time_offset;
  }

  /**
   * Return the beam response cache for a given procesing chunk. Internal
   * references are released and ownership given to the caller.
   */
  std::shared_ptr<BeamResponseCacheChunk> TakeCachedBeamResponse(size_t chunk) {
#ifdef HAVE_EVERYBEAM
    assert(beam_cache_chunks_[chunk]);
    assert(beam_cache_chunks_[chunk].get() != current_beam_cache_chunk_);
    std::shared_ptr<BeamResponseCacheChunk> beam_response =
        beam_cache_chunks_[chunk];
    beam_cache_chunks_[chunk] = nullptr;
    return beam_response;
#else
    return nullptr;
#endif
  }

  /**
   * Start processing a new chunk.
   * @ref FinishProcessingChunk() must be called at end of processing the chunk.
   * @ref FinishChunkedProcessing() must be called after processing all chunks.
   *
   * @param check_for_empty Should always be set to true when called externally.
   * It exists for internal usage only.
   */
  void StartProcessingChunk(bool check_for_empty = true) {
#ifdef HAVE_EVERYBEAM
    // Temporary access is required to the previous cache to compute the first
    // response in the new cache.
    //
    // Only add a new chunk if it will be the first chunk (check_for_empty=true)
    // Internally when called from FinishChunkedProcessing() always add a chunk
    // (check_for_empty=false)
    // This behaviour is necessary to prevent a data race with
    // TakeCachedBeamResponse() on debug builds.
    if (!check_for_empty || beam_cache_chunks_.empty()) {
      current_beam_cache_chunk_ =
          beam_cache_chunks_
              .emplace_back(std::make_shared<BeamResponseCacheChunk>())
              .get();
    }
#endif
  }

  /**
   * Finalise processing of a chunk started by @ref StartProcessingChunk().
   * NB! @ref FinishChunkedProcessing() must also be called after processing the
   * last chunk.
   */
  void FinishProcessingChunk() {
#ifdef HAVE_EVERYBEAM
    assert(!beam_cache_chunks_.empty());
    // Temporary access is required to the previous cache to compute the first
    // response in the new cache.
    previous_beam_cache_chunk_ = beam_cache_chunks_.back();
    StartProcessingChunk(false);
#endif
  }

  /**
   * Finalise any internal state after processing all data.
   * Must be called once at end of processing when using @ref
   * StartProcessingChunk()
   */
  void FinishChunkedProcessing() {
#ifdef HAVE_EVERYBEAM
    current_beam_cache_chunk_ = nullptr;
    previous_beam_cache_chunk_ = nullptr;
#endif
  }

  /**
   * Get the cached beam response data for the most recently processed timestep.
   */
  BeamResponseMap& GetCachedBeamResponse() {
#ifdef HAVE_EVERYBEAM
    return _cachedBeamResponse;
#else
    static BeamResponseMap empty_map;
    return empty_map;
#endif
  }

#ifdef HAVE_EVERYBEAM
  /**
   * Compute and cache the beam response for the provided time; if it is not
   * already cached.
   */
  void CacheBeamResponse(double time, size_t field_id,
                         const aocommon::MultiBandData& bands,
                         bool cache_entire_beam) {
    if (cache_entire_beam) {
      CacheBeamResponseWithCacheChunk(time, field_id, bands);
    } else {
      UpdateCachedBeamResponseForRow(time, field_id, bands,
                                     _cachedBeamResponse);
    }
  }
  /**
   * Compute the beam response for the provided time; if it is not already
   * cached.
   * Update the cache chunk with the new response or index to existing response.
   */
  void CacheBeamResponseWithCacheChunk(double time, size_t field_id,
                                       const aocommon::MultiBandData& bands) {
    // Add a new response to the cache if its a new time interval.
    if (UpdateCachedBeamResponseForRow(time, field_id, bands,
                                       _cachedBeamResponse)) {
      current_beam_cache_chunk_->responses.push_back(_cachedBeamResponse);
    }
    // If we have not entered a new time interval, but are processing the first
    // row of a new chunk, there is no previous reponse to reference.
    // In this specific case it becomes necessarry to retrieve/copy the last
    // response of the previous chunks cache.
    if (previous_beam_cache_chunk_) {
      if (current_beam_cache_chunk_->responses.empty()) {
        current_beam_cache_chunk_->responses.push_back(
            previous_beam_cache_chunk_->responses.back());
      }
      previous_beam_cache_chunk_ = nullptr;
    }
    // Current row indexes the most recent response.
    // Either the freshly generated response if this is a new time interval,
    // otherwise the previously cached response.
    current_beam_cache_chunk_->offsets.push_back(
        current_beam_cache_chunk_->responses.size() - 1);
  }
  /**
   * Compute the beam response for the current time if the time is in a new time
   * interval. Return true if new beam response computed; otherwise false.
   */
  bool UpdateCachedBeamResponseForRow(double time, size_t field_id,
                                      const aocommon::MultiBandData& bands,
                                      BeamResponseMap& cached_beam_response);

  template <GainMode Mode>
  void ApplyBeamResponse(std::complex<float>* data, size_t n_channels,
                         size_t data_desc_id, size_t antenna1, size_t antenna2);

  template <ModifierBehaviour Behaviour, GainMode Mode, bool ApplyForward>
  void ApplyConjugatedBeamResponseForRow(std::complex<float>* data,
                                         const float* weights,
                                         const float* image_weights,
                                         size_t n_channels, size_t data_desc_id,
                                         size_t antenna1, size_t antenna2) {
    const size_t n_visibilities = GetNVisibilities(Mode);
    const BeamResponseMap& cached_beam_response = GetCachedBeamResponse();
    for (size_t n_channel = 0; n_channel < n_channels; ++n_channel) {
      ApplyConjugatedBeamResponse<Behaviour, Mode, ApplyForward>(
          data, weights, image_weights, n_channel, n_channels, data_desc_id,
          antenna1, antenna2, cached_beam_response);
      if constexpr (internal::ShouldApplyCorrection(Behaviour)) {
        data += n_visibilities;
      }
      if constexpr (internal::ShouldSumCorrection(Behaviour)) {
        weights += n_visibilities;
      }
    }
  }
  template <ModifierBehaviour Behaviour, GainMode Mode, bool ApplyForward>
  void ApplyConjugatedBeamResponse(std::complex<float>* data,
                                   const float* weights,
                                   const float* image_weights, size_t n_channel,
                                   size_t n_channels, size_t data_desc_id,
                                   size_t antenna1, size_t antenna2,
                                   const BeamResponseMap& cached_beam_response);

  /**
   * Correct the data for both the conjugated beam and the
   * conjugated h5parm solutions.
   */
  template <ModifierBehaviour Behaviour, GainMode Mode, size_t NParms>
  void ApplyConjugatedDualForRow(std::complex<float>* data,
                                 const float* weights,
                                 const float* image_weights, size_t n_channels,
                                 size_t data_desc_id, size_t antenna1,
                                 size_t antenna2, size_t ms_index,
                                 bool apply_forward) {
    ApplyConjugatedDualForRowWithTime<Behaviour, Mode, NParms>(
        data, weights, image_weights, n_channels, data_desc_id, antenna1,
        antenna2, ms_index, apply_forward, GetTimeOffset(ms_index));
  }

  template <ModifierBehaviour Behaviour, GainMode Mode, size_t NParms>
  void ApplyConjugatedDualForRowWithTime(std::complex<float>* data,
                                         const float* weights,
                                         const float* image_weights,
                                         size_t n_channels, size_t data_desc_id,
                                         size_t antenna1, size_t antenna2,
                                         size_t ms_index, bool apply_forward,
                                         size_t time_offset) {
    const std::complex<float>* parm_response =
        _cachedParmResponse[ms_index][data_desc_id].data();
    const size_t n_visibilities = GetNVisibilities(Mode);
    const BeamResponseMap& cached_beam_response = GetCachedBeamResponse();
    for (size_t n_channel = 0; n_channel < n_channels; ++n_channel) {
      ApplyConjugatedDual<Behaviour, Mode, NParms>(
          parm_response, data, weights, image_weights, n_channel, n_channels,
          data_desc_id, antenna1, antenna2, apply_forward, time_offset,
          cached_beam_response);
      if constexpr (internal::ShouldApplyCorrection(Behaviour)) {
        data += n_visibilities;
      }
      if constexpr (internal::ShouldSumCorrection(Behaviour)) {
        weights += n_visibilities;
      }
    }
  }
  template <ModifierBehaviour Behaviour, GainMode Mode, size_t NParms>
  void ApplyConjugatedDual(const std::complex<float>* parm_response,
                           std::complex<float>* data, const float* weights,
                           const float* image_weights, size_t n_channel,
                           size_t n_channels, size_t data_desc_id,
                           size_t antenna1, size_t antenna2, bool apply_forward,
                           size_t time_offset,
                           const BeamResponseMap& cached_beam_response);
#endif

  void SetFacetDirection(double ra, double dec) {
    _facetDirectionRA = ra;
    _facetDirectionDec = dec;
  }
  double FacetDirectionRA() const { return _facetDirectionRA; }
  double FacetDirectionDec() const { return _facetDirectionDec; }
  void ResetSums() {
    correction_sum_ = AverageCorrection();
    beam_correction_sum_ = AverageCorrection();
  }
  /**
   * Sum of the full corrections applied to the visibilities. In case both
   * beam and h5parm solutions are applied, this is the weighed sum over the
   * (squared) product of both. Otherwise it is over the (squared) contribution
   * of either the beam or solutions.
   *
   * It is not an average yet, because this class doesn't store the sum of
   * weights. It would be a redundant calculation, because the gridder
   * already calculates the sum of weights.
   */
  const AverageCorrection& TotalCorrectionSum() const {
    return correction_sum_;
  }
  /**
   * In case both beam and solution gains are applied, this represents the beam
   * part of the corrections. In that case, it is the weighed sum of the squared
   * beam matrices. This is used in the final primary beam correction to
   * (approximately) separate the beam and solution parts to be able to apply a
   * smooth beam correction. If only one correction is applied, this value
   * remains zero.
   *
   * Like @ref TotalCorrectionSum(), this should be divided by the sum of
   * weights to turn it into the average correction.
   */
  const AverageCorrection& BeamCorrectionSum() const {
    return beam_correction_sum_;
  }

  /**
   * Number of complex values per solution: 4 in the case of fulljones,
   * otherwise 2 (scalar solutions are already converted to dual solutions
   * by the h5parm reader).
   */
  inline size_t NValuesPerSolution(size_t ms_index) const {
    using schaapcommon::h5parm::GainType;
    const size_t solution_index = (*_gainTypes).size() == 1 ? 0 : ms_index;
    return (*_gainTypes)[solution_index] == GainType::kFullJones ? 4 : 2;
  }

 private:
#ifdef HAVE_EVERYBEAM
  void GetBeamResponse(const aocommon::BandData& band, size_t field_id,
                       std::vector<aocommon::MC2x2F>& result) const;

  // _telescope attribute needed to keep the telecope in _pointResponse alive
  std::unique_ptr<everybeam::telescope::Telescope> _telescope;
  std::unique_ptr<everybeam::pointresponse::PointResponse> _pointResponse;
  /**
   * The beam response for the currently processed timestep. The VectorMap
   * is indexed by data_desc_id. The vector is of size n_channels x
   * n_stations, where n_stations is the fastest changing index.
   */
  BeamResponseMap _cachedBeamResponse;

  /** Hold the caches for all processed chunks in memory until they are no
   * longer required.
   * Calling @ref TakeCachedBeamResponse() will set the corresponding entry to
   * null and pass ownership to the caller via shared pointer.
   * Shared pointer is used in place of unique pointer because
   * `previous_beam_cache_chunk_` also references this data; while this
   * reference will usually be short lived it is not impossible that it outlives
   * the pointer returned by @ref TakeCachedBeamResponse(). */
  std::vector<std::shared_ptr<BeamResponseCacheChunk>> beam_cache_chunks_;
  BeamResponseCacheChunk* current_beam_cache_chunk_ = nullptr;
  /* The previous cache is kept in memory until we have processed the first row
   * of the current cache. This is because the last response of the previous
   * cache needs to be copied in some instances. */
  std::shared_ptr<BeamResponseCacheChunk> previous_beam_cache_chunk_ = nullptr;

  everybeam::BeamMode _beamMode = everybeam::BeamMode::kNone;
#endif
  /**
   * This value is retrieved from the EveryBeam telescope class, and only
   * available when the beam is applied.
   */
  size_t n_stations_ = 0;
  std::string _beamModeString;
  std::string _beamNormalisationMode;
  /**
   * Element ms_index is a VectorMap of vectors, where an element is
   * indexed by data_desc_id and contains the complex gains of
   * size n_times x n_channels x n_stations x n_parameters, where n_parameters
   * is the fastest changing. n_parameters is 2 (for diagonal) or
   * 4 (for full jones).
   *
   * The ms_index may not be a consecutive index because this gridder
   * may not have to grid all measurement sets that were specified to
   * wsclean.
   */
  std::map<size_t, SolutionsResponseMap> _cachedParmResponse;
  /**
   * Each element holds a vector with the measurement set times. The map
   * is indexed by a (non-consecutive) ms_index.
   */
  std::map<size_t, std::shared_ptr<std::vector<double>>> _cachedMSTimes;
  /**
   * Each element holds the current offset position into the _cachedParmResponse
   * and _cachedMSTimes elements of the same ms_index. The map is indexed by a
   * (non-consecutive) ms_index.
   */
  std::map<size_t, size_t> time_offsets_;
  /**
   * Optional pointers to vectors with h5parm solution objects.
   * The GriddingTaskManager, which always outlives GriddingTasks and their
   * VisibilityModifier, owns the objects.
   * If all measurement sets use the same solution, the vectors have one
   * element. Otherwise, they have one element for each ms.
   * The second solution table is optional. If both tables are used, the first
   * table has the amplitude part and the second table has the phase part.
   * @{
   */
  const std::vector<schaapcommon::h5parm::H5Parm>* _h5parms = nullptr;
  const std::vector<schaapcommon::h5parm::SolTab*>* _firstSolutions = nullptr;
  const std::vector<schaapcommon::h5parm::SolTab*>* _secondSolutions = nullptr;
  const std::vector<schaapcommon::h5parm::GainType>* _gainTypes = nullptr;
  /** @} */
  double _facetDirectionRA = 0.0;
  double _facetDirectionDec = 0.0;
  AverageCorrection correction_sum_;
  AverageCorrection beam_correction_sum_;
  double scaled_ncp_uvw_[3];
  bool uses_bda_ = false;
};

#ifdef HAVE_EVERYBEAM
template <ModifierBehaviour Behaviour, GainMode Mode, bool ApplyForward>
inline void VisibilityModifier::ApplyConjugatedBeamResponse(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channel, size_t n_channels, size_t data_desc_id, size_t antenna1,
    size_t antenna2,
    const aocommon::VectorMap<std::vector<aocommon::MC2x2F>>&
        cached_beam_response) {
  using internal::MakeDiagonalIfScalar;
  const size_t offset = n_channel * n_stations_;
  const size_t offset1 = offset + antenna1;
  const size_t offset2 = offset + antenna2;

  const aocommon::MC2x2F& gain1(cached_beam_response[data_desc_id][offset1]);
  const aocommon::MC2x2F& gain2(cached_beam_response[data_desc_id][offset2]);
  if constexpr (internal::ShouldApplyCorrection(Behaviour)) {
    if constexpr (ApplyForward) {
      internal::ApplyGain<Mode>(data, gain1, gain2);
    }
    internal::ApplyConjugatedGain<Mode>(data, gain1, gain2);
  }
  if constexpr (internal::ShouldSumCorrection(Behaviour)) {
    // This assumes that the weights of the polarizations are the same
    correction_sum_.Add<Mode>(MakeDiagonalIfScalar<Mode>(gain1),
                              MakeDiagonalIfScalar<Mode>(gain2),
                              image_weights[n_channel] * weights[0]);
  }
}

template <ModifierBehaviour Behaviour, GainMode Mode, size_t NParms>
inline void VisibilityModifier::ApplyConjugatedDual(
    const std::complex<float>* parm_response, std::complex<float>* data,
    const float* weights, const float* image_weights, size_t n_channel,
    size_t n_channels, size_t data_desc_id, size_t antenna1, size_t antenna2,
    bool apply_forward, size_t time_offset,
    const aocommon::VectorMap<std::vector<aocommon::MC2x2F>>&
        cached_beam_response) {
  using internal::CreateMatrix2x2OrDiag;
  using internal::MakeDiagonalIfScalar;

  // Get facet beam
  const size_t beam_offset = n_channel * n_stations_;
  const size_t beam_offset1 = beam_offset + antenna1;
  const size_t beam_offset2 = beam_offset + antenna2;

  const aocommon::MC2x2F& gain_b_1(
      cached_beam_response[data_desc_id][beam_offset1]);
  const aocommon::MC2x2F& gain_b_2(
      cached_beam_response[data_desc_id][beam_offset2]);

  // Get h5 solution
  // Column major indexing
  const size_t h5_offset =
      (time_offset * n_channels + n_channel) * n_stations_ * NParms;
  const size_t h5_offset1 = h5_offset + antenna1 * NParms;
  const size_t h5_offset2 = h5_offset + antenna2 * NParms;
  const auto gain_h5_1 =
      CreateMatrix2x2OrDiag<NParms>(parm_response, h5_offset1);
  const auto gain_h5_2 =
      CreateMatrix2x2OrDiag<NParms>(parm_response, h5_offset2);

  // Combine H5parm and beam. The beam is applied first on the data,
  // and therefore needs to be the last in the multiplication.
  const aocommon::MC2x2F gain_combined_1 = gain_h5_1 * gain_b_1;
  const aocommon::MC2x2F gain_combined_2 = gain_h5_2 * gain_b_2;

  if constexpr (internal::ShouldApplyCorrection(Behaviour)) {
    if (apply_forward) {
      internal::ApplyGain<Mode>(data, gain_combined_1, gain_combined_2);
    }
    internal::ApplyConjugatedGain<Mode>(data, gain_combined_1, gain_combined_2);
  }
  if constexpr (internal::ShouldSumCorrection(Behaviour)) {
    beam_correction_sum_.Add<Mode>(MakeDiagonalIfScalar<Mode>(gain_b_1),
                                   MakeDiagonalIfScalar<Mode>(gain_b_2),
                                   weights[0] * image_weights[n_channel]);
    correction_sum_.Add<Mode>(MakeDiagonalIfScalar<Mode>(gain_combined_1),
                              MakeDiagonalIfScalar<Mode>(gain_combined_2),
                              weights[0] * image_weights[n_channel]);
  }
}
#endif  // HAVE_EVERYBEAM

template <ModifierBehaviour Behaviour, GainMode Mode, size_t NParms,
          bool ApplyForward>
inline void VisibilityModifier::ApplyConjugatedParmResponse(
    const std::complex<float>* parm_response, std::complex<float>* data,
    const float* weights, const float* image_weights, size_t n_channel,
    size_t n_channels, size_t n_antennas, size_t antenna1, size_t antenna2,
    size_t time_offset) {
  using internal::CreateMatrix2x2;
  using internal::CreateMatrix2x2OrDiag;
  using internal::MakeDiagonalIfScalar;

  // Column major indexing
  const size_t offset =
      (time_offset * n_channels + n_channel) * n_antennas * NParms;
  const size_t offset1 = offset + antenna1 * NParms;
  const size_t offset2 = offset + antenna2 * NParms;
  if constexpr (internal::ShouldApplyCorrection(Behaviour)) {
    const aocommon::MC2x2F gain1 =
        CreateMatrix2x2<NParms>(parm_response, offset1);
    const aocommon::MC2x2F gain2 =
        CreateMatrix2x2<NParms>(parm_response, offset2);
    if constexpr (ApplyForward) {
      internal::ApplyGain<Mode>(data, gain1, gain2);
    }
    internal::ApplyConjugatedGain<Mode>(data, gain1, gain2);
  }
  if constexpr (internal::ShouldSumCorrection(Behaviour)) {
    // Assumes that the weights of the polarizations are the same
    const auto gain1 =
        CreateMatrix2x2OrDiag<NParms, Mode>(parm_response, offset1);
    const auto gain2 =
        CreateMatrix2x2OrDiag<NParms, Mode>(parm_response, offset2);
    correction_sum_.Add<Mode>(gain1, gain2,
                              image_weights[n_channel] * weights[0]);
  }
}

namespace internal {
inline double Sinc(double x) {
  return std::abs(x) < 1e-6 ? 1.0 : std::sin(x) / x;
}
}  // namespace internal

inline void VisibilityModifier::ApplyTimeFrequencySmearing(
    std::complex<float>* data, const double* uvw,
    const aocommon::BandData& band, int n_pol_per_vis, double dl, double dm) {
  using internal::Sinc;
  if (dl != 0.0 || dm != 0.0) {
    const double dn = std::sqrt(1.0 - dl * dl - dm * dm) - 1.0;
    const double smearing_factor_frequency =
        Sinc(M_PI * (uvw[0] * dl + uvw[1] * dm + uvw[2] * dn) /
             aocommon::kSpeedOfLight * band.ChannelWidth(0));

    // Below terms with scaled_ncp_uvw_[0] are commented out because for epoch
    // J2000 this element is zero. In case the initialization of
    // scaled_ncp_uvw_ is modified to the current epoch,
    // these terms need to be uncommented, because then they are non-zero.
    const double delta_time =
        M_PI *
        (dl * (uvw[1] * scaled_ncp_uvw_[2] - uvw[2] * scaled_ncp_uvw_[1]) +
         dm * (/* uvw[2] * scaled_ncp_uvw_[0] */ -uvw[0] * scaled_ncp_uvw_[2]) +
         dn * (uvw[0] *
               scaled_ncp_uvw_[1] /* - uvw[1] * scaled_ncp_uvw_[0] */)) /
        aocommon::kSpeedOfLight;

    const size_t n_channels = band.ChannelCount();
    const size_t n = n_pol_per_vis * n_channels;

    for (size_t i = 0; i < n; ++i) {
      data[i] *= smearing_factor_frequency *
                 Sinc(delta_time * band.ChannelFrequency(i));
    }
  }
}

}  // namespace wsclean

#endif
