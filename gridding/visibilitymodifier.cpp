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

void internal::TransposeAntennaChannelDimensions(const std::complex<float>* source,
                                                 std::complex<float>* destination, size_t n_times,
                                                 size_t n_channels, size_t n_antennas,
                                                 size_t n_parms) {
  // Transpose from (n_times * n_chan * n_ant * n_parms) to (n_times * n_ant * n_chan * n_parms)
  const std::complex<float>* source_ptr = source;
  for (size_t time = 0; time != n_times; ++time) {
    for (size_t channel = 0; channel != n_channels; ++channel) {
      for (size_t antenna = 0; antenna != n_antennas; ++antenna) {
        const size_t dest_offset = ((time * n_antennas + antenna) * n_channels + channel) * n_parms;
        std::copy_n(source_ptr, n_parms, destination + dest_offset);
        source_ptr += n_parms;
      }
    }
  }
}

void VisibilityModifier::InitializeTimeFrequencySmearing(const ObservationInfo& observation_info,
                                                         double l_shift, double m_shift) {
  cos_phase_centre_dec_ = std::cos(observation_info.phaseCentreDec);
  sin_phase_centre_dec_ = std::sin(observation_info.phaseCentreDec);

  dlmn_[0] = l_shift;
  dlmn_[1] = m_shift;
  dlmn_[2] = std::sqrt(1.0 - l_shift * l_shift - m_shift * m_shift) - 1.0;
}

void VisibilityModifier::InitializePointResponse(SynchronizedMS&& ms, double facet_beam_update_time,
                                                 const std::string& element_response,
                                                 const aocommon::MultiBandData& bands,
                                                 const std::string& data_column,
                                                 const std::string& mwa_path) {
#ifdef HAVE_EVERYBEAM
  // Hard-coded for now
  constexpr bool frequency_interpolation = true;
  constexpr bool use_channel_frequency = true;

  // Get path to coefficient file for MWA telescope
  everybeam::TelescopeType telescope_type = everybeam::GetTelescopeType(*ms);
  const std::string coeff_path = (telescope_type == everybeam::TelescopeType::kMWATelescope)
                                     ? wsclean::mwa::FindCoeffFile(mwa_path)
                                     : "";

  everybeam::ATermSettings aterm_settings;
  aterm_settings.coeff_path = coeff_path;
  aterm_settings.data_column_name = data_column;

  everybeam::Options options = everybeam::aterms::ATermConfig::ConvertToEBOptions(
      *ms, aterm_settings, frequency_interpolation, _beamNormalisationMode, use_channel_frequency,
      element_response, _beamModeString);
  _beamMode = options.beam_mode;

  _telescope = everybeam::Load(*ms, options);
  // Initialize with 0.0 time to make sure first call to UpdateTime()
  // will fill the beam response cache.
  constexpr double initial_time = 0.0;
  time_cache_ =
      std::make_unique<everybeam::TimeCache>(everybeam::CreateTimeCache(initial_time, *_telescope));
  time_cache_->SetUpdateInterval(facet_beam_update_time);
  n_stations_ = _telescope->GetNrStations();
  if (uses_frequency_bda_) {
    std::vector<aocommon::MC2x2F>& band_response = _cachedBeamResponse.AlwaysEmplace(0);
    band_response.resize(bands.MaxBandChannels() * n_stations_);
    interpolated_beam_scratch_.Resize(bands.MaxBandChannels());
  } else {
    for (size_t data_desc_id : bands.DataDescIds()) {
      std::vector<aocommon::MC2x2F>& band_response =
          _cachedBeamResponse.AlwaysEmplace(data_desc_id);
      band_response.resize(bands[data_desc_id].ChannelCount() * n_stations_);
    }
  }
#endif

  if (channel_interpolation_table_.Empty()) {
    channel_interpolation_table_ = MakeChannelInterpolationTable(bands);
    n_full_resolution_channels_ = bands.MaxBandChannels();
  }
}

void VisibilityModifier::InitializeMockResponse(
    size_t n_stations, const aocommon::MultiBandData& bands,
    const std::vector<aocommon::MC2x2F>& beam_response,
    const std::vector<std::complex<float>>& parm_response) {
  n_stations_ = n_stations;
  size_t n_channels = bands.MaxBandChannels();
  assert(beam_response.size() == n_channels * n_stations);
  assert(parm_response.size() == n_channels * n_stations * 2 ||
         parm_response.size() == n_channels * n_stations * 4);
  const size_t n_parms = parm_response.size() / (n_channels * n_stations);
  constexpr size_t kMSIndex = 0;
  SolutionsResponseMap& ms_parm_response =
      _cachedParmResponse.emplace(kMSIndex, SolutionsResponseMap()).first->second;

  std::vector<std::vector<uint32_t>> interpolation_table(bands.HighestBandId() + 1);
  uses_frequency_bda_ = bands.BandCount() > 1;
  const size_t data_desc_id = uses_frequency_bda_ ? 0 : bands.DataDescIdWithMaxChannels();
#ifdef HAVE_EVERYBEAM
  std::vector<aocommon::MC2x2F>& band_beam_response =
      _cachedBeamResponse.AlwaysEmplace(data_desc_id);
  band_beam_response.assign(beam_response.begin(), beam_response.end());
#endif
  ms_parm_response.AlwaysEmplace(data_desc_id) = parm_response;

  channel_interpolation_table_ = MakeChannelInterpolationTable(bands);
  time_offsets_ = {std::pair(0, 0)};

  n_full_resolution_channels_ = n_channels;
  interpolated_beam_scratch_.Resize(n_full_resolution_channels_);
  interpolated_solutions_scratch_.Resize(n_full_resolution_channels_ * n_parms);
}

void VisibilityModifier::InitializeCacheParmResponseForDataDescId(
    std::vector<std::complex<float>>& band_response, const std::vector<std::string>& antenna_names,
    const aocommon::BandData& band, size_t ms_index, size_t data_desc_id) {
  using schaapcommon::h5parm::JonesParameters;
  const size_t n_parms = NValuesPerSolution(ms_index);
  const size_t n_times = _cachedMSTimes[ms_index]->size();
  const size_t n_antenna = antenna_names.size();
  const size_t solution_index = _h5parms->size() == 1 ? 0 : ms_index;
  const std::vector<double> freqs(band.begin(), band.end());
  const size_t response_size = n_times * freqs.size() * n_antenna * n_parms;
  const std::string dir_name =
      (*_h5parms)[solution_index].GetNearestSource(_facetDirectionRA, _facetDirectionDec);
  schaapcommon::h5parm::SolTab* first_solution = (*_firstSolutions)[solution_index];
  schaapcommon::h5parm::SolTab* second_solution =
      _secondSolutions->empty() ? nullptr : (*_secondSolutions)[solution_index];
  const size_t dir_index = first_solution->GetDirIndex(dir_name);
  JonesParameters jones_parameters(
      freqs, *_cachedMSTimes[ms_index], antenna_names, (*_gainTypes)[solution_index],
      JonesParameters::InterpolationType::NEAREST, dir_index, first_solution, second_solution,
      false, 0u, JonesParameters::MissingAntennaBehavior::kUnit);
  const casacore::Cube<std::complex<float>>& parms = jones_parameters.GetParms();
  assert(static_cast<size_t>(parms.shape()[0]) == n_parms);
  assert(static_cast<size_t>(parms.shape()[1]) == n_antenna);
  assert(static_cast<size_t>(parms.shape()[2]) == n_times * freqs.size());
  band_response.resize(response_size);
  internal::TransposeAntennaChannelDimensions(parms.data(), band_response.data(), n_times,
                                              freqs.size(), n_antenna, n_parms);
  setNonFiniteToZero(band_response);
}

void VisibilityModifier::InitializeCacheParmResponse(const std::vector<std::string>& antenna_names,
                                                     const aocommon::MultiBandData& selected_bands,
                                                     size_t ms_index) {
  // Only extract DD solutions if the corresponding cache entry is empty.
  SolutionsResponseMap& parm_response = _cachedParmResponse[ms_index];
  if (parm_response.Empty()) {
    if (uses_frequency_bda_) {
      // In BDA mode, the channels are interpolated from the data_desc_id with
      // the most channels.
      const size_t max_channel_id = selected_bands.DataDescIdWithMaxChannels();
      const aocommon::BandData& band(selected_bands[max_channel_id]);
      const size_t n_channels = band.ChannelCount();
      std::vector<std::complex<float>>& band_response = parm_response.AlwaysEmplace(0);
      InitializeCacheParmResponseForDataDescId(band_response, antenna_names, band, ms_index,
                                               max_channel_id);

      const size_t n_parms = NValuesPerSolution(ms_index);
      interpolated_solutions_scratch_.Resize(n_channels * n_parms);
    } else {
      for (size_t data_desc_id : selected_bands.DataDescIds()) {
        const aocommon::BandData& band = selected_bands[data_desc_id];
        std::vector<std::complex<float>>& band_response = parm_response.AlwaysEmplace(data_desc_id);
        InitializeCacheParmResponseForDataDescId(band_response, antenna_names, band, ms_index,
                                                 data_desc_id);
      }
    }

    if (channel_interpolation_table_.Empty()) {
      channel_interpolation_table_ = MakeChannelInterpolationTable(selected_bands);
      n_full_resolution_channels_ = selected_bands.MaxBandChannels();
    }
  }
}

size_t VisibilityModifier::GetCacheParmResponseSize() const {
  size_t size = 0;
  for (const auto& ms_info : _cachedParmResponse) {
    for (const std::vector<std::complex<float>>& band_solutions : ms_info.second) {
      size += band_solutions.capacity();
    }
  }
  return size * sizeof(std::complex<float>);
}

void VisibilityModifier::CacheParmResponse(double time, const aocommon::MultiBandData& bands,
                                           size_t ms_index, size_t& time_offset) {
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

void VisibilityModifier::GetBeamResponse(const aocommon::BandData& band, size_t field_id,
                                         std::vector<aocommon::MC2x2F>& result) const {
  std::vector<double> frequencies(band.ChannelCount());
  for (size_t ch = 0; ch < band.ChannelCount(); ++ch) {
    frequencies[ch] = band.ChannelFrequency(ch);
  }

  everybeam::pointresponse::PointResponse& point_response = time_cache_->GetPointResponse();
  point_response.ResponseAllStations(_beamMode, result.data(), _facetDirectionRA,
                                     _facetDirectionDec, frequencies, field_id);
}

bool VisibilityModifier::UpdateCachedBeamResponseForRow(double time, size_t field_id,
                                                        const aocommon::MultiBandData& bands,
                                                        BeamResponseMap& cached_beam_response) {
  if (ShouldUpdateCachedBeamResponseForTime(time)) {
    if (uses_frequency_bda_) {
      // In BDA mode, interpolate the channels from the data_desc_id with the
      // highest nr of channels.
      const size_t max_channels_id = bands.DataDescIdWithMaxChannels();
      const aocommon::BandData& interpolation_band = bands[max_channels_id];
      std::vector<aocommon::MC2x2F>& full_cached_response = cached_beam_response[0];
      GetBeamResponse(interpolation_band, field_id, full_cached_response);
    } else {
      for (size_t data_desc_id : bands.DataDescIds()) {
        const aocommon::BandData& band = bands[data_desc_id];
        std::vector<aocommon::MC2x2F>& response = cached_beam_response[data_desc_id];
        GetBeamResponse(band, field_id, response);
      }
    }
    return true;
  }
  return false;
}

template <GainMode Mode>
void VisibilityModifier::ApplyBeamResponse(std::complex<float>* data,
                                           const aocommon::BandData& band, size_t data_desc_id,
                                           size_t antenna1, size_t antenna2) {
  const size_t n_channels = band.ChannelCount();
  const std::array<const aocommon::MC2x2F*, 2> responses =
      GetBeamResponseForRow(data_desc_id, antenna1, antenna2);
  for (size_t ch = 0; ch < n_channels; ++ch) {
    internal::ApplyGain<Mode>(data, responses[0][ch], responses[1][ch]);
    data += GetNVisibilities(Mode);
  }
}

template void VisibilityModifier::ApplyBeamResponse<GainMode::kXX>(std::complex<float>* data,
                                                                   const aocommon::BandData& band,
                                                                   size_t data_desc_id,
                                                                   size_t antenna1,
                                                                   size_t antenna2);

template void VisibilityModifier::ApplyBeamResponse<GainMode::kYY>(std::complex<float>* data,
                                                                   const aocommon::BandData& band,
                                                                   size_t data_desc_id,
                                                                   size_t antenna1,
                                                                   size_t antenna2);

template void VisibilityModifier::ApplyBeamResponse<GainMode::kTrace>(
    std::complex<float>* data, const aocommon::BandData& band, size_t data_desc_id, size_t antenna1,
    size_t antenna2);

template void VisibilityModifier::ApplyBeamResponse<GainMode::k2VisDiagonal>(
    std::complex<float>* data, const aocommon::BandData& band, size_t data_desc_id, size_t antenna1,
    size_t antenna2);

template void VisibilityModifier::ApplyBeamResponse<GainMode::kFull>(std::complex<float>* data,
                                                                     const aocommon::BandData& band,
                                                                     size_t data_desc_id,
                                                                     size_t antenna1,
                                                                     size_t antenna2);
#endif

template <GainMode Mode>
void VisibilityModifier::ApplyParmResponseForRowWithTime(std::complex<float>* data, size_t ms_index,
                                                         size_t n_channels, size_t data_desc_id,
                                                         size_t n_antenna, size_t antenna1,
                                                         size_t antenna2, size_t time_offset) {
  const size_t n_parms = NValuesPerSolution(ms_index);
  auto [parm1, parm2] =
      GetParmResponseForRow(ms_index, data_desc_id, n_antenna, antenna1, antenna2, time_offset);
  if (n_parms == 2) {
    for (size_t ch = 0; ch < n_channels; ++ch) {
      const MC2x2F gain1(parm1[0], 0, 0, parm1[1]);
      const MC2x2F gain2(parm2[0], 0, 0, parm2[1]);
      internal::ApplyGain<Mode>(data, gain1, gain2);
      data += GetNVisibilities(Mode);
      parm1 += 2;
      parm2 += 2;
    }
  } else {
    for (size_t ch = 0; ch < n_channels; ++ch) {
      const MC2x2F gain1(parm1);
      const MC2x2F gain2(parm2);
      internal::ApplyGain<Mode>(data, gain1, gain2);
      data += GetNVisibilities(Mode);
      parm1 += 4;
      parm2 += 4;
    }
  }
}

template void VisibilityModifier::ApplyParmResponseForRowWithTime<GainMode::kXX>(
    std::complex<float>* data, size_t ms_index, size_t n_channels, size_t data_desc_id,
    size_t n_antennas, size_t antenna1, size_t antenna2, size_t time_offset);

template void VisibilityModifier::ApplyParmResponseForRowWithTime<GainMode::kYY>(
    std::complex<float>* data, size_t ms_index, size_t n_channels, size_t data_desc_id,
    size_t n_antennas, size_t antenna1, size_t antenna2, size_t time_offset);

template void VisibilityModifier::ApplyParmResponseForRowWithTime<GainMode::kTrace>(
    std::complex<float>* data, size_t ms_index, size_t n_channels, size_t data_desc_id,
    size_t n_antennas, size_t antenna1, size_t antenna2, size_t time_offset);

template void VisibilityModifier::ApplyParmResponseForRowWithTime<GainMode::k2VisDiagonal>(
    std::complex<float>* data, size_t ms_index, size_t n_channels, size_t data_desc_id,
    size_t n_antennas, size_t antenna1, size_t antenna2, size_t time_offset);

template void VisibilityModifier::ApplyParmResponseForRowWithTime<GainMode::kFull>(
    std::complex<float>* data, size_t ms_index, size_t n_channels, size_t data_desc_id,
    size_t n_antennas, size_t antenna1, size_t antenna2, size_t time_offset);

}  // namespace wsclean
