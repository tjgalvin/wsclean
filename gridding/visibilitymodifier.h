#ifndef GRIDDING_VISIBILITY_MODIFIER_H_
#define GRIDDING_VISIBILITY_MODIFIER_H_

#include <stdexcept>

#ifdef HAVE_EVERYBEAM
#include <EveryBeam/beammode.h>
#include <EveryBeam/beamnormalisationmode.h>
#include <EveryBeam/pointresponse/pointresponse.h>
#endif

#include <schaapcommon/h5parm/jonesparameters.h>

#include <aocommon/banddata.h>
#include <aocommon/matrix2x2.h>
#include <aocommon/polarization.h>
#include <aocommon/uvector.h>

#include "averagecorrection.h"
#include "gainmode.h"

namespace wsclean {

class SynchronizedMS;

/**
 * Control whether to apply the modifier, sum the corrections or both when
 * calling methods that apply a modifier e.g. @ref ApplyConjugatedParmResponse()
 */
enum class ModifierBehaviour { kApply, kSum, kApplyAndSum };

/**
 * Applies beam and h5parm solutions to visibilities.
 * See the documentation for function @ref ApplyConjugatedParmResponse()
 * for an overview of parameters that hold for most of these functions.
 */
class VisibilityModifier {
 public:
  VisibilityModifier() = default;

  void InitializePointResponse(SynchronizedMS&& ms,
                               double facet_beam_update_time,
                               const std::string& element_response,
                               size_t n_channels,
                               const std::string& data_column,
                               const std::string& mwa_path);

  /**
   * A function that initializes this visibility modifier for testing. After
   * calling this function, the Apply* functions can be called.
   * @param beam_response vector of n_channels * n_stations * 4 gain elements.
   * @param parm_response vector of n_channels * n_stations * 2
   */
  void InitializeMockResponse(
      size_t n_channels, size_t n_stations,
      const std::vector<std::complex<double>>& beam_response,
      const std::vector<std::complex<float>>& parm_response);

  void SetNoPointResponse() {
#ifdef HAVE_EVERYBEAM
    _pointResponse = nullptr;
    _cachedBeamResponse.clear();
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
  }

  /**
   * @brief Cache the solutions from a h5 solution file
   */
  void InitializeCacheParmResponse(const std::vector<std::string>& antennaNames,
                                   const aocommon::BandData& band,
                                   size_t ms_index);
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
  void CacheParmResponse(double time, const aocommon::BandData& band,
                         size_t ms_index, size_t& time_offset);

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
  template <ModifierBehaviour Behaviour, GainMode Mode>
  void ApplyConjugatedParmResponse(std::complex<float>* data,
                                   const float* weights,
                                   const float* image_weights, size_t ms_index,
                                   size_t n_channels, size_t n_antennas,
                                   size_t antenna1, size_t antenna2,
                                   bool apply_forward) {
    ApplyConjugatedParmResponse<Behaviour, Mode>(
        data, weights, image_weights, ms_index, n_channels, n_antennas,
        antenna1, antenna2, apply_forward, GetTimeOffset(ms_index));
  }
  template <ModifierBehaviour Behaviour, GainMode Mode>
  void ApplyConjugatedParmResponse(std::complex<float>* data,
                                   const float* weights,
                                   const float* image_weights, size_t ms_index,
                                   size_t n_channels, size_t n_antennas,
                                   size_t antenna1, size_t antenna2,
                                   bool apply_forward, size_t time_offset);

  template <GainMode GainEntry>
  void ApplyParmResponse(std::complex<float>* data, size_t ms_index,
                         size_t n_channels, size_t n_antennas, size_t antenna1,
                         size_t antenna2) {
    ApplyParmResponse<GainEntry>(data, ms_index, n_channels, n_antennas,
                                 antenna1, antenna2, GetTimeOffset(ms_index));
  }
  template <GainMode GainEntry>
  void ApplyParmResponse(std::complex<float>* data, size_t ms_index,
                         size_t n_channels, size_t n_antennas, size_t antenna1,
                         size_t antenna2, size_t time_offset);

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

  /*
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

#ifdef HAVE_EVERYBEAM
  /**
   * @brief Compute and cache the beam response if no cached response
   * present for the provided time.
   */
  void CacheBeamResponse(double time, size_t fieldId,
                         const aocommon::BandData& band);

  template <GainMode Mode>
  void ApplyBeamResponse(std::complex<float>* data, size_t n_channels,
                         size_t antenna1, size_t antenna2);

  template <ModifierBehaviour Behaviour, GainMode Mode>
  void ApplyConjugatedBeamResponse(std::complex<float>* data,
                                   const float* weights,
                                   const float* image_weights,
                                   size_t n_channels, size_t antenna1,
                                   size_t antenna2, bool apply_forward);

  /**
   * Correct the data for both the conjugated beam and the
   * conjugated h5parm solutions.
   */
  template <ModifierBehaviour Behaviour, GainMode Mode>
  void ApplyConjugatedDual(std::complex<float>* data, const float* weights,
                           const float* image_weights, size_t n_channels,
                           size_t n_stations, size_t antenna1, size_t antenna2,
                           size_t ms_index, bool apply_forward) {
    ApplyConjugatedDual<Behaviour, Mode>(
        data, weights, image_weights, n_channels, n_stations, antenna1,
        antenna2, ms_index, apply_forward, GetTimeOffset(ms_index));
  }
  template <ModifierBehaviour Behaviour, GainMode Mode>
  void ApplyConjugatedDual(std::complex<float>* data, const float* weights,
                           const float* image_weights, size_t n_channels,
                           size_t n_stations, size_t antenna1, size_t antenna2,
                           size_t ms_index, bool apply_forward,
                           size_t time_offset);
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
  // _telescope attribute needed to keep the telecope in _pointResponse alive
  std::unique_ptr<everybeam::telescope::Telescope> _telescope;
  std::unique_ptr<everybeam::pointresponse::PointResponse> _pointResponse;
  /**
   * The beam response for the currently processed timestep.
   * It's of size n_channels x _pointResponseBufferSize, which equals
   * n_channels x n_stations x n_elements(=4), where n_elements is the fastest
   * changing index.
   */
  aocommon::UVector<std::complex<float>> _cachedBeamResponse;
  everybeam::BeamMode _beamMode = everybeam::BeamMode::kNone;
#endif
  std::string _beamModeString;
  std::string _beamNormalisationMode;
  /**
   * Element ms_index is a vector of complex gains of
   * size n_times x n_channels x n_stations x n_parameters, where n_parameters
   * is the fastest changing. n_parameters is 2 (for diagonal) or
   * 4 (for full jones).
   *
   * The ms_index may not be a consecutive index because this gridder
   * may not have to grid all measurement sets that were specified to
   * wsclean.
   */
  std::map<size_t, std::vector<std::complex<float>>> _cachedParmResponse;
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
  size_t _pointResponseBufferSize = 0;
  double _facetDirectionRA = 0.0;
  double _facetDirectionDec = 0.0;
  AverageCorrection correction_sum_;
  AverageCorrection beam_correction_sum_;
};

namespace internal {
/**
 * @brief Apply gains to the visibilities.
 *
 * @tparam PolarizationCount polarization count, 2 or 4 for IDG, 1 for all other
 * gridders.
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
}  // namespace internal

#ifdef HAVE_EVERYBEAM
template <ModifierBehaviour Behaviour, GainMode Mode>
inline void VisibilityModifier::ApplyConjugatedBeamResponse(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t antenna1, size_t antenna2, bool apply_forward) {
  using internal::MakeDiagonalIfScalar;
  for (size_t ch = 0; ch < n_channels; ++ch) {
    const size_t offset = ch * _pointResponseBufferSize;
    const size_t offset1 = offset + antenna1 * 4u;
    const size_t offset2 = offset + antenna2 * 4u;

    const aocommon::MC2x2F gain1(&_cachedBeamResponse[offset1]);
    const aocommon::MC2x2F gain2(&_cachedBeamResponse[offset2]);
    if constexpr (internal::ShouldApplyCorrection(Behaviour)) {
      if (apply_forward) {
        internal::ApplyGain<Mode>(data, gain1, gain2);
      }
      internal::ApplyConjugatedGain<Mode>(data, gain1, gain2);
      data += GetNVisibilities(Mode);
    }
    if constexpr (internal::ShouldSumCorrection(Behaviour)) {
      // This assumes that the weights of the polarizations are the same
      correction_sum_.Add<Mode>(MakeDiagonalIfScalar<Mode>(gain1),
                                MakeDiagonalIfScalar<Mode>(gain2),
                                image_weights[ch] * weights[0]);
      weights += GetNVisibilities(Mode);
    }
  }
}

template <ModifierBehaviour Behaviour, GainMode Mode>
inline void VisibilityModifier::ApplyConjugatedDual(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t n_channels, size_t n_stations, size_t antenna1, size_t antenna2,
    size_t ms_index, bool apply_forward, size_t time_offset) {
  using internal::MakeDiagonalIfScalar;
  const size_t nparms = NValuesPerSolution(ms_index);
  const std::vector<std::complex<float>>& parm_response =
      _cachedParmResponse[ms_index];

  if (nparms == 2) {
    for (size_t ch = 0; ch < n_channels; ++ch) {
      // Compute facet beam
      const size_t beam_offset = ch * _pointResponseBufferSize;
      const size_t beam_offset1 = beam_offset + antenna1 * 4u;
      const size_t beam_offset2 = beam_offset + antenna2 * 4u;

      const aocommon::MC2x2F gain_b_1(&_cachedBeamResponse[beam_offset1]);
      const aocommon::MC2x2F gain_b_2(&_cachedBeamResponse[beam_offset2]);

      // Get H5 solutions
      // Column major indexing
      const size_t h5_offset =
          (time_offset * n_channels + ch) * n_stations * nparms;
      const size_t h5_offset1 = h5_offset + antenna1 * nparms;
      const size_t h5_offset2 = h5_offset + antenna2 * nparms;
      const aocommon::MC2x2FDiag gain_h5_1(parm_response[h5_offset1],
                                           parm_response[h5_offset1 + 1]);
      const aocommon::MC2x2FDiag gain_h5_2(parm_response[h5_offset2],
                                           parm_response[h5_offset2 + 1]);

      // Combine H5parm and beam. The beam is applied first on the data,
      // and therefore needs to be the last in the multiplication.
      const aocommon::MC2x2F gain_combined_1 = gain_h5_1 * gain_b_1;
      const aocommon::MC2x2F gain_combined_2 = gain_h5_2 * gain_b_2;

      if constexpr (internal::ShouldApplyCorrection(Behaviour)) {
        if (apply_forward) {
          internal::ApplyGain<Mode>(data, gain_combined_1, gain_combined_2);
        }
        internal::ApplyConjugatedGain<Mode>(data, gain_combined_1,
                                            gain_combined_2);
        data += GetNVisibilities(Mode);
      }
      if constexpr (internal::ShouldSumCorrection(Behaviour)) {
        beam_correction_sum_.Add<Mode>(MakeDiagonalIfScalar<Mode>(gain_b_1),
                                       MakeDiagonalIfScalar<Mode>(gain_b_2),
                                       weights[0] * image_weights[ch]);
        correction_sum_.Add<Mode>(MakeDiagonalIfScalar<Mode>(gain_combined_1),
                                  MakeDiagonalIfScalar<Mode>(gain_combined_2),
                                  weights[0] * image_weights[ch]);
        weights += GetNVisibilities(Mode);
      }
    }
  } else {
    // This branch handles full jones H5 parm files (nparms == 4)
    for (size_t ch = 0; ch < n_channels; ++ch) {
      // Get facet beam
      const size_t beam_offset = ch * _pointResponseBufferSize;
      const size_t beam_offset1 = beam_offset + antenna1 * 4u;
      const size_t beam_offset2 = beam_offset + antenna2 * 4u;

      const aocommon::MC2x2F gain_b_1(&_cachedBeamResponse[beam_offset1]);
      const aocommon::MC2x2F gain_b_2(&_cachedBeamResponse[beam_offset2]);

      // Get h5 solution
      // Column major indexing
      const size_t offset_h5 =
          (time_offset * n_channels + ch) * n_stations * nparms;
      const size_t offset_h5_1 = offset_h5 + antenna1 * nparms;
      const size_t offset_h5_2 = offset_h5 + antenna2 * nparms;
      const aocommon::MC2x2F gain_h5_1(&parm_response[offset_h5_1]);
      const aocommon::MC2x2F gain_h5_2(&parm_response[offset_h5_2]);

      // Combine H5parm and beam. The beam is applied first on the data,
      // and therefore needs to be the last in the multiplication.
      const aocommon::MC2x2F gain_combined_1 = gain_h5_1 * gain_b_1;
      const aocommon::MC2x2F gain_combined_2 = gain_h5_2 * gain_b_2;

      if constexpr (internal::ShouldApplyCorrection(Behaviour)) {
        internal::ApplyConjugatedGain<Mode>(data, gain_combined_1,
                                            gain_combined_2);
        if (apply_forward) {
          internal::ApplyGain<Mode>(data, gain_combined_1, gain_combined_2);
        }
        data += GetNVisibilities(Mode);
      }
      if constexpr (internal::ShouldSumCorrection(Behaviour)) {
        beam_correction_sum_.Add<Mode>(MakeDiagonalIfScalar<Mode>(gain_b_1),
                                       MakeDiagonalIfScalar<Mode>(gain_b_2),
                                       weights[0] * image_weights[ch]);
        correction_sum_.Add<Mode>(MakeDiagonalIfScalar<Mode>(gain_combined_1),
                                  MakeDiagonalIfScalar<Mode>(gain_combined_2),
                                  weights[0] * image_weights[ch]);
        weights += GetNVisibilities(Mode);
      }
    }
  }
}
#endif  // HAVE_EVERYBEAM

template <ModifierBehaviour Behaviour, GainMode Mode>
inline void VisibilityModifier::ApplyConjugatedParmResponse(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, size_t n_channels, size_t n_antennas, size_t antenna1,
    size_t antenna2, bool apply_forward, size_t time_offset) {
  using internal::MakeDiagonalIfScalar;
  const size_t nparms = NValuesPerSolution(ms_index);
  const std::vector<std::complex<float>>& parm_response =
      _cachedParmResponse[ms_index];

  // Conditional could be templated once C++ supports partial function
  // specialization
  if (nparms == 2) {
    for (size_t ch = 0; ch < n_channels; ++ch) {
      // Column major indexing
      const size_t offset =
          (time_offset * n_channels + ch) * n_antennas * nparms;
      const size_t offset1 = offset + antenna1 * nparms;
      const size_t offset2 = offset + antenna2 * nparms;
      if constexpr (internal::ShouldApplyCorrection(Behaviour)) {
        const aocommon::MC2x2F gain1(parm_response[offset1], 0.0f, 0.0f,
                                     parm_response[offset1 + 1]);
        const aocommon::MC2x2F gain2(parm_response[offset2], 0.0f, 0.0f,
                                     parm_response[offset2 + 1]);
        if (apply_forward) {
          internal::ApplyGain<Mode>(data, gain1, gain2);
        }
        internal::ApplyConjugatedGain<Mode>(data, gain1, gain2);
        data += GetNVisibilities(Mode);
      }
      if constexpr (internal::ShouldSumCorrection(Behaviour)) {
        const aocommon::MC2x2FDiag gain1(parm_response[offset1],
                                         parm_response[offset1 + 1]);
        const aocommon::MC2x2FDiag gain2(parm_response[offset2],
                                         parm_response[offset2 + 1]);
        correction_sum_.Add<Mode>(gain1, gain2, image_weights[ch] * weights[0]);
        weights += GetNVisibilities(Mode);
      }
    }
  } else {
    for (size_t ch = 0; ch < n_channels; ++ch) {
      // Column major indexing
      const size_t offset =
          (time_offset * n_channels + ch) * n_antennas * nparms;
      const size_t offset1 = offset + antenna1 * nparms;
      const size_t offset2 = offset + antenna2 * nparms;
      const aocommon::MC2x2F gain1(&parm_response[offset1]);
      const aocommon::MC2x2F gain2(&parm_response[offset2]);
      if constexpr (internal::ShouldApplyCorrection(Behaviour)) {
        if (apply_forward) {
          internal::ApplyGain<Mode>(data, gain1, gain2);
        }
        internal::ApplyConjugatedGain<Mode>(data, gain1, gain2);
        data += GetNVisibilities(Mode);
      }
      if constexpr (internal::ShouldSumCorrection(Behaviour)) {
        // Assumes that the weights of the polarizations are the same
        correction_sum_.Add<Mode>(MakeDiagonalIfScalar<Mode>(gain1),
                                  MakeDiagonalIfScalar<Mode>(gain2),
                                  image_weights[ch] * weights[0]);
        weights += GetNVisibilities(Mode);
      }
    }
  }
}

}  // namespace wsclean

#endif
