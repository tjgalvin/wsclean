#include "visibilitymodifier.h"

#include "../msproviders/synchronizedms.h"
#include "../structures/observationinfo.h"

#ifdef HAVE_EVERYBEAM
#include <EveryBeam/load.h>
#include <EveryBeam/aterms/atermconfig.h>
#include <EveryBeam/pointresponse/phasedarraypoint.h>
#endif

// Only needed for EB/H5Parm related options
#include "../io/findmwacoefffile.h"

#include "solutionchannelinterpolation.h"

using aocommon::MC2x2F;

namespace wsclean {

namespace {
void setNonFiniteToZero(std::vector<std::complex<float>>& values) {
  for (std::complex<float>& v : values) {
    if (!std::isfinite(v.real()) || !std::isfinite(v.imag())) {
      v = 0.0;
    }
  }
}
}  // namespace

void VisibilityModifier::InitializeTimeFrequencySmearing(
    const ObservationInfo& observation_info, double interval) {
  /**
   * Length of a sidereal day in seconds.
   */
  constexpr double kSiderealDay = 86164.0905;

  const double angular_speed = 2 * M_PI * interval / kSiderealDay;

  // scaled_ncp_uvw_ is initialized with the NCP for epoch J2000.
  // The expression is simpler than for the the current epoch (of the
  // observation). The first element is zero and that is currently assumed where
  // scaled_ncp_uvw_ is used. If the initialization here is modified to the more
  // accurate current epoch with a non-zero first element, please update the
  // usage as well
  scaled_ncp_uvw_[0] = 0.0;
  scaled_ncp_uvw_[1] =
      angular_speed * std::cos(observation_info.phaseCentreDec);
  scaled_ncp_uvw_[2] =
      angular_speed * std::sin(observation_info.phaseCentreDec);
}

void VisibilityModifier::InitializePointResponse(
    SynchronizedMS&& ms, double facet_beam_update_time,
    const std::string& element_response, const aocommon::MultiBandData& bands,
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
  n_stations_ = _telescope->GetNrStations();
  for (size_t data_desc_id : bands.DataDescIds()) {
    std::vector<aocommon::MC2x2F>& band_response =
        _cachedBeamResponse.AlwaysEmplace(data_desc_id);
    band_response.resize(bands[data_desc_id].ChannelCount() * n_stations_);
  }
#endif
}

void VisibilityModifier::InitializeMockResponse(
    size_t n_channels, size_t n_stations, size_t data_desc_id,
    const std::vector<aocommon::MC2x2F>& beam_response,
    const std::vector<std::complex<float>>& parm_response) {
  n_stations_ = n_stations;
  assert(beam_response.size() == n_channels * n_stations);
  assert(parm_response.size() == n_channels * n_stations * 2 ||
         parm_response.size() == n_channels * n_stations * 4);
#ifdef HAVE_EVERYBEAM
  std::vector<aocommon::MC2x2F>& band_beam_response =
      _cachedBeamResponse.AlwaysEmplace(data_desc_id);
  band_beam_response.assign(beam_response.begin(), beam_response.end());
#endif
  SolutionsResponseMap& band_parm_response =
      _cachedParmResponse.emplace(0, SolutionsResponseMap()).first->second;
  band_parm_response.AlwaysEmplace(data_desc_id) = parm_response;
  time_offsets_ = {std::pair(0, 0)};
}

void VisibilityModifier::InitializeCacheParmResponse(
    const std::vector<std::string>& antennaNames,
    const aocommon::MultiBandData& selected_bands, size_t ms_index) {
  using schaapcommon::h5parm::JonesParameters;

  const size_t solution_index = (*_h5parms).size() == 1 ? 0 : ms_index;

  // Only extract DD solutions if the corresponding cache entry is empty.
  SolutionsResponseMap& parm_response = _cachedParmResponse[ms_index];
  if (parm_response.Empty()) {
    const size_t nparms = NValuesPerSolution(ms_index);
    const size_t n_times = _cachedMSTimes[ms_index]->size();
    const std::string dir_name = (*_h5parms)[solution_index].GetNearestSource(
        _facetDirectionRA, _facetDirectionDec);
    schaapcommon::h5parm::SolTab* const first_solution =
        (*_firstSolutions)[solution_index];
    schaapcommon::h5parm::SolTab* const second_solution =
        _secondSolutions->empty() ? nullptr
                                  : (*_secondSolutions)[solution_index];
    const size_t dir_index = first_solution->GetDirIndex(dir_name);

    if (uses_bda_) {
      // In BDA mode, the channels are interpolated from the data_desc_id with
      // the most channels, to avoid rereading the file many times.
      const size_t max_channel_id = selected_bands.DataDescIdWithMaxChannels();
      const aocommon::BandData& band = selected_bands[max_channel_id];
      const std::vector<double> freqs(band.begin(), band.end());
      const size_t response_size =
          n_times * freqs.size() * antennaNames.size() * nparms;
      JonesParameters jones_parameters(
          freqs, *_cachedMSTimes[ms_index], antennaNames,
          (*_gainTypes)[solution_index],
          JonesParameters::InterpolationType::NEAREST, dir_index,
          first_solution, second_solution, false, 0u,
          JonesParameters::MissingAntennaBehavior::kUnit);
      // parms (Casacore::Cube) is column major
      const casacore::Cube<std::complex<float>>& parms =
          jones_parameters.GetParms();

      std::vector<std::complex<float>> full_channel_response(
          &parms(0, 0, 0), &parms(0, 0, 0) + response_size);
      for (size_t data_desc_id : selected_bands.DataDescIds()) {
        std::vector<std::complex<float>>& band_response =
            parm_response.AlwaysEmplace(data_desc_id);
        InterpolateChannels(full_channel_response, band, band_response,
                            selected_bands[data_desc_id], n_times,
                            antennaNames.size() * nparms);
        setNonFiniteToZero(band_response);
      }
    } else {
      for (size_t data_desc_id : selected_bands.DataDescIds()) {
        const aocommon::BandData& band = selected_bands[data_desc_id];
        const std::vector<double> freqs(band.begin(), band.end());
        const size_t response_size =
            n_times * freqs.size() * antennaNames.size() * nparms;
        JonesParameters jones_parameters(
            freqs, *_cachedMSTimes[ms_index], antennaNames,
            (*_gainTypes)[solution_index],
            JonesParameters::InterpolationType::NEAREST, dir_index,
            first_solution, second_solution, false, 0u,
            JonesParameters::MissingAntennaBehavior::kUnit);

        const casacore::Cube<std::complex<float>>& parms =
            jones_parameters.GetParms();
        std::vector<std::complex<float>>& band_response =
            parm_response.AlwaysEmplace(data_desc_id);
        band_response.assign(&parms(0, 0, 0), &parms(0, 0, 0) + response_size);
        setNonFiniteToZero(band_response);
      }
    }
  }
}

size_t VisibilityModifier::GetCacheParmResponseSize() const {
  size_t size = 0;
  for (const auto& ms_info : _cachedParmResponse) {
    for (const std::vector<std::complex<float>>& band_solutions :
         ms_info.second) {
      size += band_solutions.capacity();
    }
  }
  return size * sizeof(std::complex<float>);
}

void VisibilityModifier::CacheParmResponse(double time,
                                           const aocommon::MultiBandData& bands,
                                           size_t ms_index,
                                           size_t& time_offset) {
  const auto it = std::find(_cachedMSTimes[ms_index]->begin() + time_offset,
                            _cachedMSTimes[ms_index]->end(), time);
  if (it != _cachedMSTimes[ms_index]->end()) {
    // Update time_offsets_ value with index
    time_offset = std::distance(_cachedMSTimes[ms_index]->begin(), it);
  } else {
    throw std::runtime_error(
        "Time not found in cached times. A potential reason could be that the "
        "time values in the provided MS are not in ascending order. "
        "Error occurred with ms index = " +
        std::to_string(ms_index) +
        ", cache "
        "contained " +
        std::to_string(_cachedMSTimes[ms_index]->size()) + " elements.\n");
  }
}

#ifdef HAVE_EVERYBEAM

void VisibilityModifier::GetBeamResponse(
    const aocommon::BandData& band, size_t field_id,
    std::vector<aocommon::MC2x2F>& result) const {
  std::vector<double> frequencies(band.ChannelCount());
  std::vector<aocommon::MC2x2F> response(n_stations_ * band.ChannelCount());
  for (size_t ch = 0; ch < band.ChannelCount(); ++ch) {
    frequencies[ch] = band.ChannelFrequency(ch);
  }
  _pointResponse->ResponseAllStations(_beamMode, response.data(),
                                      _facetDirectionRA, _facetDirectionDec,
                                      frequencies, field_id);
  // Transpose the structure
  for (size_t station = 0; station != n_stations_; ++station) {
    for (size_t ch = 0; ch < band.ChannelCount(); ++ch) {
      result[ch * n_stations_ + station] =
          response[station * band.ChannelCount() + ch];
    }
  }
}

bool VisibilityModifier::UpdateCachedBeamResponseForRow(
    double time, size_t field_id, const aocommon::MultiBandData& bands,
    BeamResponseMap& cached_beam_response) {
  _pointResponse->UpdateTime(time);
  if (_pointResponse->HasTimeUpdate()) {
    if (uses_bda_) {
      // In BDA mode, interpolate the channels from the data_desc_id with the
      // highest nr of channels.
      const size_t max_channels_id = bands.DataDescIdWithMaxChannels();
      const aocommon::BandData& band = bands[max_channels_id];
      std::vector<aocommon::MC2x2F>& full_cached_response =
          cached_beam_response[max_channels_id];
      GetBeamResponse(band, field_id, full_cached_response);
      for (size_t data_desc_id : bands.DataDescIds()) {
        std::vector<aocommon::MC2x2F>& destination =
            cached_beam_response[data_desc_id];
        const aocommon::BandData& destination_band = bands[data_desc_id];
        InterpolateChannels(full_cached_response, band, destination,
                            destination_band, 1, n_stations_);
      }
    } else {
      for (size_t data_desc_id : bands.DataDescIds()) {
        const aocommon::BandData& band = bands[data_desc_id];
        std::vector<aocommon::MC2x2F>& response =
            cached_beam_response[data_desc_id];
        GetBeamResponse(band, field_id, response);
      }
    }
    return true;
  }
  return false;
}

template <GainMode Mode>
void VisibilityModifier::ApplyBeamResponse(std::complex<float>* data,
                                           size_t n_channels,
                                           size_t data_desc_id, size_t antenna1,
                                           size_t antenna2) {
  for (size_t ch = 0; ch < n_channels; ++ch) {
    const size_t offset = ch * n_stations_;
    const size_t offset1 = offset + antenna1;
    const size_t offset2 = offset + antenna2;

    const MC2x2F& gain1(_cachedBeamResponse[data_desc_id][offset1]);
    const MC2x2F& gain2(_cachedBeamResponse[data_desc_id][offset2]);
    internal::ApplyGain<Mode>(data, gain1, gain2);
    data += GetNVisibilities(Mode);
  }
}

template void VisibilityModifier::ApplyBeamResponse<GainMode::kXX>(
    std::complex<float>* data, size_t n_channels, size_t data_desc_id,
    size_t antenna1, size_t antenna2);

template void VisibilityModifier::ApplyBeamResponse<GainMode::kYY>(
    std::complex<float>* data, size_t n_channels, size_t data_desc_id,
    size_t antenna1, size_t antenna2);

template void VisibilityModifier::ApplyBeamResponse<GainMode::kTrace>(
    std::complex<float>* data, size_t n_channels, size_t data_desc_id,
    size_t antenna1, size_t antenna2);

template void VisibilityModifier::ApplyBeamResponse<GainMode::k2VisDiagonal>(
    std::complex<float>* data, size_t n_channels, size_t data_desc_id,
    size_t antenna1, size_t antenna2);

template void VisibilityModifier::ApplyBeamResponse<GainMode::kFull>(
    std::complex<float>* data, size_t n_channels, size_t data_desc_id,
    size_t antenna1, size_t antenna2);
#endif

template <GainMode Mode>
void VisibilityModifier::ApplyParmResponseWithTime(
    std::complex<float>* data, size_t ms_index, size_t n_channels,
    size_t data_desc_id, size_t n_antennas, size_t antenna1, size_t antenna2,
    size_t time_offset) {
  const size_t nparms = NValuesPerSolution(ms_index);
  const std::vector<std::complex<float>>& parm_response =
      _cachedParmResponse[ms_index][data_desc_id];
  if (nparms == 2) {
    for (size_t ch = 0; ch < n_channels; ++ch) {
      // Column major indexing
      const size_t offset =
          (time_offset * n_channels + ch) * n_antennas * nparms;
      const size_t offset1 = offset + antenna1 * nparms;
      const size_t offset2 = offset + antenna2 * nparms;
      const MC2x2F gain1(parm_response[offset1], 0, 0,
                         parm_response[offset1 + 1]);
      const MC2x2F gain2(parm_response[offset2], 0, 0,
                         parm_response[offset2 + 1]);
      internal::ApplyGain<Mode>(data, gain1, gain2);
      data += GetNVisibilities(Mode);
    }
  } else {
    for (size_t ch = 0; ch < n_channels; ++ch) {
      // Column major indexing
      const size_t offset =
          (time_offset * n_channels + ch) * n_antennas * nparms;
      const size_t offset1 = offset + antenna1 * nparms;
      const size_t offset2 = offset + antenna2 * nparms;
      const MC2x2F gain1(&parm_response[offset1]);
      const MC2x2F gain2(&parm_response[offset2]);
      internal::ApplyGain<Mode>(data, gain1, gain2);
      data += GetNVisibilities(Mode);
    }
  }
}

template void VisibilityModifier::ApplyParmResponseWithTime<GainMode::kXX>(
    std::complex<float>* data, size_t ms_index, size_t n_channels,
    size_t data_desc_id, size_t n_antennas, size_t antenna1, size_t antenna2,
    size_t time_offset);

template void VisibilityModifier::ApplyParmResponseWithTime<GainMode::kYY>(
    std::complex<float>* data, size_t ms_index, size_t n_channels,
    size_t data_desc_id, size_t n_antennas, size_t antenna1, size_t antenna2,
    size_t time_offset);

template void VisibilityModifier::ApplyParmResponseWithTime<GainMode::kTrace>(
    std::complex<float>* data, size_t ms_index, size_t n_channels,
    size_t data_desc_id, size_t n_antennas, size_t antenna1, size_t antenna2,
    size_t time_offset);

template void
VisibilityModifier::ApplyParmResponseWithTime<GainMode::k2VisDiagonal>(
    std::complex<float>* data, size_t ms_index, size_t n_channels,
    size_t data_desc_id, size_t n_antennas, size_t antenna1, size_t antenna2,
    size_t time_offset);

template void VisibilityModifier::ApplyParmResponseWithTime<GainMode::kFull>(
    std::complex<float>* data, size_t ms_index, size_t n_channels,
    size_t data_desc_id, size_t n_antennas, size_t antenna1, size_t antenna2,
    size_t time_offset);

}  // namespace wsclean
