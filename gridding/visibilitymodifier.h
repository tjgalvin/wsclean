#ifndef GRIDDING_VISIBILITY_MODIFIER_H_
#define GRIDDING_VISIBILITY_MODIFIER_H_

#ifdef HAVE_EVERYBEAM
#include <EveryBeam/beammode.h>
#include <EveryBeam/beamnormalisationmode.h>
#include <EveryBeam/pointresponse/pointresponse.h>
#endif

#include <schaapcommon/h5parm/jonesparameters.h>

#include <aocommon/banddata.h>
#include <aocommon/uvector.h>

class SynchronizedMS;

/**
 * Enum for selecting the entry or entries from the direction dependent gain
 * matrix that are to be used to correct the visibilities during the reading
 * and/or writing operations.
 */
enum class DDGainMatrix { kXX, kYY, kTrace, kFull };

class VisibilityModifier {
 public:
  VisibilityModifier() = default;

  void InitializePointResponse(SynchronizedMS&& ms,
                               double facet_beam_update_time,
                               const std::string& element_response,
                               size_t n_channels,
                               const std::string& data_column,
                               const std::string& mwa_path);

  void SetNoPointResponse() {
#ifdef HAVE_EVERYBEAM
    _pointResponse = nullptr;
    _cachedBeamResponse.resize(0);
#endif
  }

  void SetBeamInfo(std::string mode, std::string normalisation) {
#ifdef HAVE_EVERYBEAM
    _beamModeString = std::move(mode);
    _beamNormalisationMode = std::move(normalisation);
#endif
  }

  void ResetCache(size_t n_measurement_sets,
                  const std::vector<std::string>& solutionFiles,
                  const std::vector<std::string>& solutionTables) {
    // Assign, rather than a resize here to make sure that
    // caches are re-initialized - even in the case an MSGridderBase
    // object would be re-used for a multiple gridding tasks.
    _cachedParmResponse.assign(n_measurement_sets, {});
    _cachedMSTimes.assign(n_measurement_sets, {});
    _timeOffset.assign(n_measurement_sets, 0u);
    SetH5Parms(n_measurement_sets, solutionFiles, solutionTables);
  }

  /**
   * @brief Cache the solutions from a h5 solution file and update the
   * associated time.
   */
  void CacheParmResponse(double time,
                         const std::vector<std::string>& antennaNames,
                         const aocommon::BandData& band, size_t ms_index);

  template <size_t PolarizationCount, DDGainMatrix GainEntry>
  void CorrectParmResponse(std::complex<float>* data, size_t ms_index,
                           const aocommon::BandData& band, size_t n_antennas,
                           size_t antenna1, size_t antenna2);

  template <size_t PolarizationCount, DDGainMatrix GainEntry>
  float CorrectConjugatedParmResponse(std::complex<float>* data,
                                      const float* weights,
                                      const float* image_weights,
                                      size_t ms_index,
                                      const aocommon::BandData& band,
                                      size_t n_antennas, size_t antenna1,
                                      size_t antenna2, bool apply_forward);

  void SetMSTimes(size_t ms_index, std::vector<double>&& times) {
    _cachedMSTimes[ms_index] = std::move(times);
  }

  bool HasH5Parm() const { return !_h5parms.empty(); }

  struct DualResult {
    float h5Sum = 0.0;
    float correctionSum = 0.0;
  };

#ifdef HAVE_EVERYBEAM
  /**
   * @brief Compute and cache the beam response if no cached response
   * present for the provided time.
   */
  void CacheBeamResponse(double time, size_t fieldId,
                         const aocommon::BandData& band);

  template <size_t PolarizationCount, DDGainMatrix GainEntry>
  void ApplyBeamResponse(std::complex<float>* data,
                         const aocommon::BandData& band, size_t antenna1,
                         size_t antenna2);

  template <size_t PolarizationCount, DDGainMatrix GainEntry>
  float ApplyConjugatedBeamResponse(std::complex<float>* data,
                                    const float* weights,
                                    const float* image_weights,
                                    const aocommon::BandData& band,
                                    size_t antenna1, size_t antenna2,
                                    bool apply_forward);

  /**
   * Correct the data for both the conjugated beam and the
   * conjugated h5parm solutions.
   */
  template <size_t PolarizationCount, DDGainMatrix GainEntry>
  DualResult ApplyConjugatedDual(std::complex<float>* data,
                                 const float* weights,
                                 const float* image_weights,
                                 const aocommon::BandData& band,
                                 const std::vector<std::string>& antennaNames,
                                 size_t antenna1, size_t antenna2,
                                 size_t ms_index, bool apply_forward);
#endif

  void SetFacetDirection(double ra, double dec) {
    _facetDirectionRA = ra;
    _facetDirectionDec = dec;
  }
  double FacetDirectionRA() const { return _facetDirectionRA; }
  double FacetDirectionDec() const { return _facetDirectionDec; }

 private:
  void SetH5Parms(size_t n_measurement_sets,
                  const std::vector<std::string>& solutionFiles,
                  const std::vector<std::string>& solutionTables);

#ifdef HAVE_EVERYBEAM
  // _telescope attribute needed to keep the telecope in _pointResponse alive
  std::unique_ptr<everybeam::telescope::Telescope> _telescope;
  std::unique_ptr<everybeam::pointresponse::PointResponse> _pointResponse;
  aocommon::UVector<std::complex<float>> _cachedBeamResponse;
  everybeam::BeamMode _beamMode = everybeam::BeamMode::kNone;
#endif
  std::string _beamModeString;
  std::string _beamNormalisationMode;
  std::vector<std::vector<std::complex<float>>> _cachedParmResponse;
  std::vector<std::unique_ptr<schaapcommon::h5parm::H5Parm>> _h5parms;
  std::vector<
      std::pair<schaapcommon::h5parm::SolTab*, schaapcommon::h5parm::SolTab*>>
      _h5SolTabs;
  std::vector<schaapcommon::h5parm::JonesParameters::CorrectType> _correctType;
  std::vector<std::vector<double>> _cachedMSTimes;
  std::vector<size_t> _timeOffset;
  double _facetDirectionRA = 0.0;
  double _facetDirectionDec = 0.0;
};

#endif
