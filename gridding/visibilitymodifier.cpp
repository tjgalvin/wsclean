#include "visibilitymodifier.h"

#include "../msproviders/synchronizedms.h"

#ifdef HAVE_EVERYBEAM
#include <EveryBeam/load.h>
#include <EveryBeam/aterms/atermconfig.h>
#include <EveryBeam/pointresponse/phasedarraypoint.h>
#endif

// Only needed for EB/H5Parm related options
#include "../io/findmwacoefffile.h"

#include <aocommon/matrix2x2.h>

namespace {
void setNonFiniteToZero(std::vector<std::complex<float>>& values) {
  for (std::complex<float>& v : values) {
    if (!std::isfinite(v.real()) || !std::isfinite(v.imag())) {
      v = 0.0;
    }
  }
}

/**
 * @brief Compute the gain based from the given gain matrices.
 *
 * @tparam GainEntry Which entry or entries from the gain matrices should be
 * taken into account? See GainMode for further documentation.
 */
template <GainMode GainEntry>
std::complex<float> ComputeGain(const aocommon::MC2x2F& gain1,
                                const aocommon::MC2x2F& gain2);

template <>
std::complex<float> ComputeGain<GainMode::kXX>(const aocommon::MC2x2F& gain1,
                                               const aocommon::MC2x2F& gain2) {
  return gain2[0] * std::conj(gain1[0]);
}

template <>
std::complex<float> ComputeGain<GainMode::kYY>(const aocommon::MC2x2F& gain1,
                                               const aocommon::MC2x2F& gain2) {
  return gain2[3] * std::conj(gain1[3]);
}

template <>
std::complex<float> ComputeGain<GainMode::kDiagonal>(
    const aocommon::MC2x2F& gain1, const aocommon::MC2x2F& gain2) {
  return 0.5f * gain2.DoubleDot(gain1.Conjugate());
}

template <>
std::complex<float> ComputeGain<GainMode::kFull>(
    [[maybe_unused]] const aocommon::MC2x2F& gain1,
    [[maybe_unused]] const aocommon::MC2x2F& gain2) {
  throw std::runtime_error("ComputeGain<GainMode::kFull> not implemented!");
}

/**
 * @brief Apply conjugated gains to the visibilities.
 *
 * @tparam PolarizationCount polarization count, 2 or 4 for IDG, 1 for all other
 * gridders.
 * @tparam GainEntry Which entry or entries from the gain matrices should be
 * taken into account when correcting the visibilities? See also the
 * documentation of GainMode.
 */
template <size_t PolarizationCount, GainMode GainEntry>
void ApplyConjugatedGain(std::complex<float>* visibilities,
                         const aocommon::MC2x2F& gain1,
                         const aocommon::MC2x2F& gain2);

template <>
void ApplyConjugatedGain<1, GainMode::kXX>(std::complex<float>* visibilities,
                                           const aocommon::MC2x2F& gain1,
                                           const aocommon::MC2x2F& gain2) {
  *visibilities = std::conj(gain1[0]) * (*visibilities) * gain2[0];
}

template <>
void ApplyConjugatedGain<1, GainMode::kYY>(std::complex<float>* visibilities,
                                           const aocommon::MC2x2F& gain1,
                                           const aocommon::MC2x2F& gain2) {
  *visibilities = std::conj(gain1[3]) * (*visibilities) * gain2[3];
}

template <>
void ApplyConjugatedGain<1, GainMode::kDiagonal>(
    std::complex<float>* visibilities, const aocommon::MC2x2F& gain1,
    const aocommon::MC2x2F& gain2) {
  // Stokes-I
  *visibilities *= 0.5f * gain2.DoubleDot(gain1.Conjugate());
}

template <>
void ApplyConjugatedGain<2, GainMode::kFull>(std::complex<float>*,
                                             const aocommon::MC2x2F&,
                                             const aocommon::MC2x2F&) {
  throw std::runtime_error(
      "ApplyConjugatedGain<2, GainMode::kFull> not implemented");
}

template <>
void ApplyConjugatedGain<4, GainMode::kFull>(std::complex<float>* visibilities,
                                             const aocommon::MC2x2F& gain1,
                                             const aocommon::MC2x2F& gain2) {
  // All polarizations
  const aocommon::MC2x2F visibilities_mc2x2(visibilities);
  const aocommon::MC2x2F result =
      gain1.HermThenMultiply(visibilities_mc2x2).Multiply(gain2);
  result.AssignTo(visibilities);
}

/**
 * @brief Apply gains to the visibilities.
 *
 * @tparam PolarizationCount polarization count, 2 or 4 for IDG, 1 for all other
 * gridders.
 * @tparam GainEntry Which entry or entries from the gain matrices should be
 * taken into account when correcting the visibilities? See also the
 * documentation of GainMode.
 */
template <size_t PolarizationCount, GainMode GainEntry>
void ApplyGain(std::complex<float>* visibilities, const aocommon::MC2x2F& gain1,
               const aocommon::MC2x2F& gain2);

template <>
void ApplyGain<1, GainMode::kXX>(std::complex<float>* visibilities,
                                 const aocommon::MC2x2F& gain1,
                                 const aocommon::MC2x2F& gain2) {
  *visibilities = gain1[0] * (*visibilities) * std::conj(gain2[0]);
}

template <>
void ApplyGain<1, GainMode::kYY>(std::complex<float>* visibilities,
                                 const aocommon::MC2x2F& gain1,
                                 const aocommon::MC2x2F& gain2) {
  *visibilities = gain1[3] * (*visibilities) * std::conj(gain2[3]);
}

template <>
void ApplyGain<1, GainMode::kDiagonal>(std::complex<float>* visibilities,
                                       const aocommon::MC2x2F& gain1,
                                       const aocommon::MC2x2F& gain2) {
  // Stokes-I.
  *visibilities *= 0.5f * gain1.DoubleDot(gain2.Conjugate());
}

template <>
void ApplyGain<2, GainMode::kFull>(std::complex<float>*,
                                   const aocommon::MC2x2F&,
                                   const aocommon::MC2x2F&) {
  throw std::runtime_error("ApplyGain<2, GainMode::kFull> not implemented");
}

template <>
void ApplyGain<4, GainMode::kFull>(std::complex<float>* visibilities,
                                   const aocommon::MC2x2F& gain1,
                                   const aocommon::MC2x2F& gain2) {
  // All polarizations
  const aocommon::MC2x2F visibilities_mc2x2(visibilities);
  const aocommon::MC2x2F result =
      gain1.Multiply(visibilities_mc2x2).MultiplyHerm(gain2);
  result.AssignTo(visibilities);
}

}  // namespace

void VisibilityModifier::InitializePointResponse(
    SynchronizedMS&& ms, double facet_beam_update_time,
    const std::string& element_response, size_t n_channels,
    const std::string& data_column, const std::string& mwa_path) {
#ifdef HAVE_EVERYBEAM
  // Hard-coded for now
  const bool frequency_interpolation = true;
  const bool use_channel_frequency = true;

  // Get path to coefficient file for MWA telescope
  everybeam::TelescopeType telescope_type = everybeam::GetTelescopeType(*ms);
  const std::string coeff_path =
      (telescope_type == everybeam::TelescopeType::kMWATelescope)
          ? wsclean::mwa::FindCoeffFile(mwa_path)
          : "";

  everybeam::ATermSettings aterm_settings;
  aterm_settings.coeff_path = coeff_path;
  aterm_settings.data_column_name = data_column;

  everybeam::Options options =
      everybeam::aterms::ATermConfig::ConvertToEBOptions(
          *ms, aterm_settings, frequency_interpolation, _beamNormalisationMode,
          use_channel_frequency, element_response, _beamModeString);
  _beamMode = options.beam_mode;

  _telescope = everybeam::Load(*ms, options);
  // Initialize with 0.0 time to make sure first call to UpdateTime()
  // will fill the beam response cache.
  _pointResponse = _telescope->GetPointResponse(0.0);
  _pointResponse->SetUpdateInterval(facet_beam_update_time);
  _cachedBeamResponse.resize(n_channels *
                             _pointResponse->GetAllStationsBufferSize());
#endif
}

void VisibilityModifier::CacheParmResponse(
    double time, const std::vector<std::string>& antennaNames,
    const aocommon::BandData& band, size_t ms_index) {
  // Only extract DD solutions if the corresponding cache entry is empty.
  if (_cachedParmResponse[ms_index].empty()) {
    const size_t nparms =
        (_correctType[ms_index] ==
         schaapcommon::h5parm::JonesParameters::CorrectType::FULLJONES)
            ? 4
            : 2;
    const std::vector<double> freqs(band.begin(), band.end());
    const size_t responseSize = _cachedMSTimes[ms_index].size() * freqs.size() *
                                antennaNames.size() * nparms;
    const std::string dirName = _h5parms[ms_index]->GetNearestSource(
        _facetDirectionRA, _facetDirectionDec);
    const size_t dirIndex = _h5SolTabs[ms_index].first->GetDirIndex(dirName);
    schaapcommon::h5parm::JonesParameters jonesParameters(
        freqs, _cachedMSTimes[ms_index], antennaNames, _correctType[ms_index],
        schaapcommon::h5parm::JonesParameters::InterpolationType::NEAREST,
        dirIndex, _h5SolTabs[ms_index].first, _h5SolTabs[ms_index].second,
        false, 0.0f, 0u,
        schaapcommon::h5parm::JonesParameters::MissingAntennaBehavior::kUnit);
    // parms (Casacore::Cube) is column major
    const casacore::Cube<std::complex<float>>& parms =
        jonesParameters.GetParms();
    _cachedParmResponse[ms_index].assign(&parms(0, 0, 0),
                                         &parms(0, 0, 0) + responseSize);
    setNonFiniteToZero(_cachedParmResponse[ms_index]);
  }

  auto it = std::find(_cachedMSTimes[ms_index].begin() + _timeOffset[ms_index],
                      _cachedMSTimes[ms_index].end(), time);
  if (it != _cachedMSTimes[ms_index].end()) {
    // Update _timeOffset value with index
    _timeOffset[ms_index] = std::distance(_cachedMSTimes[ms_index].begin(), it);
  } else {
    throw std::runtime_error(
        "Time not found in cached times. A potential reason could be that the "
        "time values in the provided MS are not in ascending order.");
  }
}

#ifdef HAVE_EVERYBEAM
void VisibilityModifier::CacheBeamResponse(double time, size_t fieldId,
                                           const aocommon::BandData& band) {
  _pointResponse->UpdateTime(time);
  if (_pointResponse->HasTimeUpdate()) {
    for (size_t ch = 0; ch < band.ChannelCount(); ++ch) {
      _pointResponse->ResponseAllStations(
          _beamMode,
          &_cachedBeamResponse[ch * _pointResponse->GetAllStationsBufferSize()],
          _facetDirectionRA, _facetDirectionDec, band.ChannelFrequency(ch),
          fieldId);
    }
  }
}

template <size_t PolarizationCount, GainMode GainEntry>
void VisibilityModifier::ApplyBeamResponse(std::complex<float>* data,
                                           const aocommon::BandData& band,
                                           size_t antenna1, size_t antenna2) {
  for (size_t ch = 0; ch < band.ChannelCount(); ++ch) {
    const size_t offset = ch * _pointResponse->GetAllStationsBufferSize();
    const size_t offset1 = offset + antenna1 * 4u;
    const size_t offset2 = offset + antenna2 * 4u;

    const aocommon::MC2x2F gain1(&_cachedBeamResponse[offset1]);
    const aocommon::MC2x2F gain2(&_cachedBeamResponse[offset2]);
    ApplyGain<PolarizationCount, GainEntry>(data, gain1, gain2);
    data += PolarizationCount;
  }
}

template void VisibilityModifier::ApplyBeamResponse<1, GainMode::kXX>(
    std::complex<float>* data, const aocommon::BandData& band, size_t antenna1,
    size_t antenna2);

template void VisibilityModifier::ApplyBeamResponse<1, GainMode::kYY>(
    std::complex<float>* data, const aocommon::BandData& band, size_t antenna1,
    size_t antenna2);

template void VisibilityModifier::ApplyBeamResponse<1, GainMode::kDiagonal>(
    std::complex<float>* data, const aocommon::BandData& band, size_t antenna1,
    size_t antenna2);

template void VisibilityModifier::ApplyBeamResponse<2, GainMode::kFull>(
    std::complex<float>* data, const aocommon::BandData& band, size_t antenna1,
    size_t antenna2);

template void VisibilityModifier::ApplyBeamResponse<4, GainMode::kFull>(
    std::complex<float>* data, const aocommon::BandData& band, size_t antenna1,
    size_t antenna2);

template <size_t PolarizationCount, GainMode GainEntry>
float VisibilityModifier::ApplyConjugatedBeamResponse(
    std::complex<float>* data, const float* weights, const float* image_weights,
    const aocommon::BandData& band, size_t antenna1, size_t antenna2,
    bool apply_forward) {
  float correction_sum = 0.0;
  for (size_t ch = 0; ch < band.ChannelCount(); ++ch) {
    const size_t offset = ch * _pointResponse->GetAllStationsBufferSize();
    const size_t offset1 = offset + antenna1 * 4u;
    const size_t offset2 = offset + antenna2 * 4u;

    const aocommon::MC2x2F gain1(&_cachedBeamResponse[offset1]);
    const aocommon::MC2x2F gain2(&_cachedBeamResponse[offset2]);
    if (apply_forward) {
      ApplyGain<PolarizationCount, GainEntry>(data, gain1, gain2);
    }
    ApplyConjugatedGain<PolarizationCount, GainEntry>(data, gain1, gain2);
    const std::complex<float> g = ComputeGain<GainEntry>(gain1, gain2);

    const float weight = *weights * image_weights[ch];
    correction_sum += (conj(g) * weight * g).real();

    // Only admissible PolarizationCount for applying the facet beam is 1.
    data += PolarizationCount;
    weights += PolarizationCount;
  }
  return correction_sum;
}

template float
VisibilityModifier::ApplyConjugatedBeamResponse<1, GainMode::kXX>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    const aocommon::BandData& band, size_t antenna1, size_t antenna2,
    bool apply_forward);

template float
VisibilityModifier::ApplyConjugatedBeamResponse<1, GainMode::kYY>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    const aocommon::BandData& band, size_t antenna1, size_t antenna2,
    bool apply_forward);

template float
VisibilityModifier::ApplyConjugatedBeamResponse<1, GainMode::kDiagonal>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    const aocommon::BandData& band, size_t antenna1, size_t antenna2,
    bool apply_forward);

template float
VisibilityModifier::ApplyConjugatedBeamResponse<2, GainMode::kFull>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    const aocommon::BandData& band, size_t antenna1, size_t antenna2,
    bool apply_forward);

template float
VisibilityModifier::ApplyConjugatedBeamResponse<4, GainMode::kFull>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    const aocommon::BandData& band, size_t antenna1, size_t antenna2,
    bool apply_forward);

#endif

void VisibilityModifier::SetH5Parms(
    size_t n_measurement_sets, const std::vector<std::string>& solutionFiles,
    const std::vector<std::string>& solutionTables) {
  _h5parms.resize(n_measurement_sets);
  _h5SolTabs.resize(n_measurement_sets);
  _correctType.resize(n_measurement_sets);

  for (size_t i = 0; i < _h5parms.size(); ++i) {
    if (solutionFiles.size() > 1) {
      _h5parms[i].reset(new schaapcommon::h5parm::H5Parm(solutionFiles[i]));
    } else {
      _h5parms[i].reset(new schaapcommon::h5parm::H5Parm(solutionFiles[0]));
    }

    if (solutionTables.size() == 1) {
      _h5SolTabs[i] =
          std::make_pair(&_h5parms[i]->GetSolTab(solutionTables[0]), nullptr);
      _correctType[i] =
          schaapcommon::h5parm::JonesParameters::StringToCorrectType(
              _h5SolTabs[i].first->GetType());
    } else if (solutionTables.size() == 2) {
      const std::array<std::string, 2> solTabTypes{
          _h5parms[i]->GetSolTab(solutionTables[0]).GetType(),
          _h5parms[i]->GetSolTab(solutionTables[1]).GetType()};

      const std::array<std::string, 2>::const_iterator amplitude_iter =
          std::find(solTabTypes.begin(), solTabTypes.end(), "amplitude");
      const std::array<std::string, 2>::const_iterator phase_iter =
          std::find(solTabTypes.begin(), solTabTypes.end(), "phase");

      if (amplitude_iter == solTabTypes.end() ||
          phase_iter == solTabTypes.end()) {
        throw std::runtime_error(
            "WSClean expects solution tables with name 'amplitude' and "
            "'phase', but received " +
            solTabTypes[0] + " and " + solTabTypes[1]);
      } else {
        const size_t amplitude_index =
            std::distance(solTabTypes.begin(), amplitude_iter);
        const size_t phase_index =
            std::distance(solTabTypes.begin(), phase_iter);
        _h5SolTabs[i] = std::make_pair(
            &_h5parms[i]->GetSolTab(solutionTables[amplitude_index]),
            &_h5parms[i]->GetSolTab(solutionTables[phase_index]));
      }

      const size_t n_amplitude_pol =
          _h5SolTabs[i].first->HasAxis("pol")
              ? _h5SolTabs[i].first->GetAxis("pol").size
              : 1;
      const size_t n_phase_pol = _h5SolTabs[i].second->HasAxis("pol")
                                     ? _h5SolTabs[i].second->GetAxis("pol").size
                                     : 1;
      if (n_amplitude_pol == 1 && n_phase_pol == 1) {
        _correctType[i] =
            schaapcommon::h5parm::JonesParameters::CorrectType::SCALARGAIN;
      } else if (n_amplitude_pol == 2 && n_phase_pol == 2) {
        _correctType[i] =
            schaapcommon::h5parm::JonesParameters::CorrectType::GAIN;
      } else if (n_amplitude_pol == 4 && n_phase_pol == 4) {
        _correctType[i] =
            schaapcommon::h5parm::JonesParameters::CorrectType::FULLJONES;
      } else {
        throw std::runtime_error(
            "Incorrect or mismatching number of polarizations in the "
            "provided amplitude and phase soltabs. The number of polarizations "
            "should both be either 1, 2 or 4, but received " +
            std::to_string(n_amplitude_pol) + " for amplitude and " +
            std::to_string(n_phase_pol) + " for phase");
      }
    } else {
      throw std::runtime_error(
          "Specify the solution table name(s) with "
          "-soltab-names=soltabname1[OPTIONAL,soltabname2]");
    }
  }
}

template <size_t PolarizationCount, GainMode GainEntry>
void VisibilityModifier::CorrectParmResponse(std::complex<float>* data,
                                             size_t ms_index,
                                             const aocommon::BandData& band,
                                             size_t n_antennas, size_t antenna1,
                                             size_t antenna2) {
  const size_t nparms =
      (_correctType[ms_index] ==
       schaapcommon::h5parm::JonesParameters::CorrectType::FULLJONES)
          ? 4
          : 2;
  const size_t nchannels = band.ChannelCount();
  if (nparms == 2) {
    for (size_t ch = 0; ch < nchannels; ++ch) {
      // Column major indexing
      const size_t offset =
          (_timeOffset[ms_index] * nchannels + ch) * n_antennas * nparms;
      const size_t offset1 = offset + antenna1 * nparms;
      const size_t offset2 = offset + antenna2 * nparms;
      const aocommon::MC2x2F gain1(_cachedParmResponse[ms_index][offset1], 0, 0,
                                   _cachedParmResponse[ms_index][offset1 + 1]);
      const aocommon::MC2x2F gain2(_cachedParmResponse[ms_index][offset2], 0, 0,
                                   _cachedParmResponse[ms_index][offset2 + 1]);
      ApplyGain<PolarizationCount, GainEntry>(data, gain1, gain2);
      data += PolarizationCount;
    }
  } else {
    for (size_t ch = 0; ch < nchannels; ++ch) {
      // Column major indexing
      const size_t offset =
          (_timeOffset[ms_index] * nchannels + ch) * n_antennas * nparms;
      const size_t offset1 = offset + antenna1 * nparms;
      const size_t offset2 = offset + antenna2 * nparms;
      const aocommon::MC2x2F gain1(&_cachedParmResponse[ms_index][offset1]);
      const aocommon::MC2x2F gain2(&_cachedParmResponse[ms_index][offset2]);
      ApplyGain<PolarizationCount, GainEntry>(data, gain1, gain2);
      data += PolarizationCount;
    }
  }
}

template void VisibilityModifier::CorrectParmResponse<1, GainMode::kXX>(
    std::complex<float>* data, size_t ms_index, const aocommon::BandData& band,
    size_t n_antennas, size_t antenna1, size_t antenna2);

template void VisibilityModifier::CorrectParmResponse<1, GainMode::kYY>(
    std::complex<float>* data, size_t ms_index, const aocommon::BandData& band,
    size_t n_antennas, size_t antenna1, size_t antenna2);

template void VisibilityModifier::CorrectParmResponse<1, GainMode::kDiagonal>(
    std::complex<float>* data, size_t ms_index, const aocommon::BandData& band,
    size_t n_antennas, size_t antenna1, size_t antenna2);

template void VisibilityModifier::CorrectParmResponse<2, GainMode::kFull>(
    std::complex<float>* data, size_t ms_index, const aocommon::BandData& band,
    size_t n_antennas, size_t antenna1, size_t antenna2);

template void VisibilityModifier::CorrectParmResponse<4, GainMode::kFull>(
    std::complex<float>* data, size_t ms_index, const aocommon::BandData& band,
    size_t n_antennas, size_t antenna1, size_t antenna2);

template <size_t PolarizationCount, GainMode GainEntry>
float VisibilityModifier::CorrectConjugatedParmResponse(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, const aocommon::BandData& band, size_t n_antennas,
    size_t antenna1, size_t antenna2, bool apply_forward) {
  const size_t nparms =
      (_correctType[ms_index] ==
       schaapcommon::h5parm::JonesParameters::CorrectType::FULLJONES)
          ? 4
          : 2;
  const size_t nchannels = band.ChannelCount();

  float correctionSum = 0.0;

  // Conditional could be templated once C++ supports partial function
  // specialization
  if (nparms == 2) {
    for (size_t ch = 0; ch < nchannels; ++ch) {
      // Column major indexing
      const size_t offset =
          (_timeOffset[ms_index] * nchannels + ch) * n_antennas * nparms;
      const size_t offset1 = offset + antenna1 * nparms;
      const size_t offset2 = offset + antenna2 * nparms;
      const aocommon::MC2x2F gain1(_cachedParmResponse[ms_index][offset1], 0, 0,
                                   _cachedParmResponse[ms_index][offset1 + 1]);
      const aocommon::MC2x2F gain2(_cachedParmResponse[ms_index][offset2], 0, 0,
                                   _cachedParmResponse[ms_index][offset2 + 1]);
      if (apply_forward) {
        ApplyGain<PolarizationCount, GainEntry>(data, gain1, gain2);
      }
      ApplyConjugatedGain<PolarizationCount, GainEntry>(data, gain1, gain2);
      const std::complex<float> g = ComputeGain<GainEntry>(gain1, gain2);

      const float weight = *weights * image_weights[ch];
      correctionSum += (conj(g) * weight * g).real();

      // Only admissible PolarizationCount for applying gains from solution
      // file is 1.
      data += PolarizationCount;
      weights += PolarizationCount;
    }
  } else {
    for (size_t ch = 0; ch < nchannels; ++ch) {
      // Column major indexing
      const size_t offset =
          (_timeOffset[ms_index] * nchannels + ch) * n_antennas * nparms;
      const size_t offset1 = offset + antenna1 * nparms;
      const size_t offset2 = offset + antenna2 * nparms;
      const aocommon::MC2x2F gain1(&_cachedParmResponse[ms_index][offset1]);
      const aocommon::MC2x2F gain2(&_cachedParmResponse[ms_index][offset2]);
      if (apply_forward) {
        ApplyGain<PolarizationCount, GainEntry>(data, gain1, gain2);
      }
      ApplyConjugatedGain<PolarizationCount, GainEntry>(data, gain1, gain2);
      const std::complex<float> g = ComputeGain<GainEntry>(gain1, gain2);

      const float weight = *weights * image_weights[ch];
      correctionSum += (conj(g) * weight * g).real();

      // Only admissible PolarizationCount for applying gains from solution
      // file is 1.
      data += PolarizationCount;
      weights += PolarizationCount;
    }
  }
  return correctionSum;
}

template float
VisibilityModifier::CorrectConjugatedParmResponse<1, GainMode::kXX>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, const aocommon::BandData& band, size_t n_antennas,
    size_t antenna1, size_t antenna2, bool apply_forward);

template float
VisibilityModifier::CorrectConjugatedParmResponse<1, GainMode::kYY>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, const aocommon::BandData& band, size_t n_antennas,
    size_t antenna1, size_t antenna2, bool apply_forward);

template float
VisibilityModifier::CorrectConjugatedParmResponse<1, GainMode::kDiagonal>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, const aocommon::BandData& band, size_t n_antennas,
    size_t antenna1, size_t antenna2, bool apply_forward);

template float
VisibilityModifier::CorrectConjugatedParmResponse<2, GainMode::kFull>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, const aocommon::BandData& band, size_t n_antennas,
    size_t antenna1, size_t antenna2, bool apply_forward);

template float
VisibilityModifier::CorrectConjugatedParmResponse<4, GainMode::kFull>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    size_t ms_index, const aocommon::BandData& band, size_t n_antennas,
    size_t antenna1, size_t antenna2, bool apply_forward);

#ifdef HAVE_EVERYBEAM
template <size_t PolarizationCount, GainMode GainEntry>
VisibilityModifier::DualResult VisibilityModifier::ApplyConjugatedDual(
    std::complex<float>* data, const float* weights, const float* image_weights,
    const aocommon::BandData& band,
    const std::vector<std::string>& antennaNames, size_t antenna1,
    size_t antenna2, size_t ms_index, bool apply_forward) {
  DualResult result;
  const size_t nparms =
      (_correctType[ms_index] ==
       schaapcommon::h5parm::JonesParameters::CorrectType::FULLJONES)
          ? 4
          : 2;
  const size_t nchannels = band.ChannelCount();

  // Conditional could be templated once C++ supports partial function
  // specialization
  if (nparms == 2) {
    for (size_t ch = 0; ch < nchannels; ++ch) {
      // Apply facet beam
      const size_t beam_offset =
          ch * _pointResponse->GetAllStationsBufferSize();
      const size_t beam_offset1 = beam_offset + antenna1 * 4u;
      const size_t beam_offset2 = beam_offset + antenna2 * 4u;

      const aocommon::MC2x2F gain_b_1(&_cachedBeamResponse[beam_offset1]);
      const aocommon::MC2x2F gain_b_2(&_cachedBeamResponse[beam_offset2]);
      if (apply_forward) {
        ApplyGain<PolarizationCount, GainEntry>(data, gain_b_1, gain_b_2);
      }
      ApplyConjugatedGain<PolarizationCount, GainEntry>(data, gain_b_1,
                                                        gain_b_2);
      const std::complex<float> g_b =
          ComputeGain<GainEntry>(gain_b_1, gain_b_2);

      // Apply H5 solutions
      // Column major indexing
      const size_t h5_offset = (_timeOffset[ms_index] * nchannels + ch) *
                               antennaNames.size() * nparms;
      const size_t h5_offset1 = h5_offset + antenna1 * nparms;
      const size_t h5_offset2 = h5_offset + antenna2 * nparms;
      const aocommon::MC2x2F gain_h5_1(
          _cachedParmResponse[ms_index][h5_offset1], 0, 0,
          _cachedParmResponse[ms_index][h5_offset1 + 1]);
      const aocommon::MC2x2F gain_h5_2(
          _cachedParmResponse[ms_index][h5_offset2], 0, 0,
          _cachedParmResponse[ms_index][h5_offset2 + 1]);
      if (apply_forward) {
        ApplyGain<PolarizationCount, GainEntry>(data, gain_h5_1, gain_h5_2);
      }
      ApplyConjugatedGain<PolarizationCount, GainEntry>(data, gain_h5_1,
                                                        gain_h5_2);
      const std::complex<float> g_h5 =
          ComputeGain<GainEntry>(gain_h5_1, gain_h5_2);

      const float weight = *weights * image_weights[ch];
      // h5Sum needed for the primary beam correction
      // This is equivalent to: real(conj(g_h5) * weight * g_h5)
      result.h5Sum += std::norm(g_h5) * weight;
      // This is equivalent to: real(conj(g_h5) * conj(g_b) * weight * g_b *
      // g_h5)
      result.correctionSum += std::norm(g_h5) * std::norm(g_b) * weight;

      // Only admissible PolarizationCount for applying gains from solution
      // file is 1.
      data += PolarizationCount;
      weights += PolarizationCount;
    }
  } else {
    for (size_t ch = 0; ch < nchannels; ++ch) {
      // Apply facet beam
      const size_t beam_offset =
          ch * _pointResponse->GetAllStationsBufferSize();
      const size_t beam_offset1 = beam_offset + antenna1 * 4u;
      const size_t beam_offset2 = beam_offset + antenna2 * 4u;

      const aocommon::MC2x2F gain_b_1(&_cachedBeamResponse[beam_offset1]);
      const aocommon::MC2x2F gain_b_2(&_cachedBeamResponse[beam_offset2]);
      if (apply_forward) {
        ApplyGain<PolarizationCount, GainEntry>(data, gain_b_1, gain_b_2);
      }
      ApplyConjugatedGain<PolarizationCount, GainEntry>(data, gain_b_1,
                                                        gain_b_2);
      const std::complex<float> g_b =
          ComputeGain<GainEntry>(gain_b_1, gain_b_2);

      // Apply h5
      // Column major indexing
      const size_t offset_h5 = (_timeOffset[ms_index] * nchannels + ch) *
                               antennaNames.size() * nparms;
      const size_t offset_h5_1 = offset_h5 + antenna1 * nparms;
      const size_t offset_h5_2 = offset_h5 + antenna2 * nparms;
      const aocommon::MC2x2F gain_h5_1(
          &_cachedParmResponse[ms_index][offset_h5_1]);
      const aocommon::MC2x2F gain_h5_2(
          &_cachedParmResponse[ms_index][offset_h5_2]);
      if (apply_forward) {
        ApplyGain<PolarizationCount, GainEntry>(data, gain_h5_1, gain_h5_2);
      }
      ApplyConjugatedGain<PolarizationCount, GainEntry>(data, gain_h5_1,
                                                        gain_h5_2);
      const std::complex<float> g_h5 =
          ComputeGain<GainEntry>(gain_h5_1, gain_h5_2);

      // Compute weight
      const float weight = *weights * image_weights[ch];
      // h5Sum needed for the primary beam correction
      // This is equivalent to real(conj(g_h5) * weight * g_h5)
      result.h5Sum += std::norm(g_h5) * weight;
      // This is equivalent to real(conj(g_h5) * conj(g_b) * weight * g_b *
      // g_h5)
      result.correctionSum += std::norm(g_h5) * std::norm(g_b) * weight;

      // Only admissible PolarizationCount for applying gains from solution
      // file is 1.
      data += PolarizationCount;
      weights += PolarizationCount;
    }
  }
  return result;
}

template VisibilityModifier::DualResult
VisibilityModifier::ApplyConjugatedDual<1, GainMode::kXX>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    const aocommon::BandData& band,
    const std::vector<std::string>& antennaNames, size_t antenna1,
    size_t antenna2, size_t ms_index, bool apply_forward);

template VisibilityModifier::DualResult
VisibilityModifier::ApplyConjugatedDual<1, GainMode::kYY>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    const aocommon::BandData& band,
    const std::vector<std::string>& antennaNames, size_t antenna1,
    size_t antenna2, size_t ms_index, bool apply_forward);

template VisibilityModifier::DualResult
VisibilityModifier::ApplyConjugatedDual<1, GainMode::kDiagonal>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    const aocommon::BandData& band,
    const std::vector<std::string>& antennaNames, size_t antenna1,
    size_t antenna2, size_t ms_index, bool apply_forward);

template VisibilityModifier::DualResult
VisibilityModifier::ApplyConjugatedDual<2, GainMode::kFull>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    const aocommon::BandData& band,
    const std::vector<std::string>& antennaNames, size_t antenna1,
    size_t antenna2, size_t ms_index, bool apply_forward);

template VisibilityModifier::DualResult
VisibilityModifier::ApplyConjugatedDual<4, GainMode::kFull>(
    std::complex<float>* data, const float* weights, const float* image_weights,
    const aocommon::BandData& band,
    const std::vector<std::string>& antennaNames, size_t antenna1,
    size_t antenna2, size_t ms_index, bool apply_forward);

#endif  // HAVE_EVERYBEAM
