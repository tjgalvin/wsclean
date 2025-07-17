#include "visibilitymodifier.h"

#include "../msproviders/synchronizedms.h"

#ifdef HAVE_EVERYBEAM
#include <EveryBeam/load.h>
#include <EveryBeam/aterms/atermconfig.h>
#include <EveryBeam/pointresponse/phasedarraypoint.h>
#endif

// Only needed for EB/H5Parm related options
#include "../io/findmwacoefffile.h"

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
  _pointResponseBufferSize = _pointResponse->GetAllStationsBufferSize();
  for (size_t data_desc_id : bands.DataDescIds()) {
    aocommon::UVector<std::complex<float>>& band_response =
        _cachedBeamResponse.AlwaysEmplace(data_desc_id);
    band_response.resize(bands[data_desc_id].ChannelCount() *
                         _pointResponseBufferSize);
  }
#endif
}

void VisibilityModifier::InitializeMockResponse(
    size_t n_channels, size_t n_stations, size_t data_desc_id,
    const std::vector<std::complex<double>>& beam_response,
    const std::vector<std::complex<float>>& parm_response) {
  _pointResponseBufferSize = n_stations * 4;
  assert(beam_response.size() == n_channels * _pointResponseBufferSize);
  assert(parm_response.size() == n_channels * n_stations * 2 ||
         parm_response.size() == n_channels * n_stations * 4);
#ifdef HAVE_EVERYBEAM
  aocommon::UVector<std::complex<float>>& band_beam_response =
      _cachedBeamResponse.AlwaysEmplace(data_desc_id);
  band_beam_response.assign(beam_response.begin(), beam_response.end());
#endif
  aocommon::VectorMap<std::vector<std::complex<float>>>& band_parm_response =
      _cachedParmResponse
          .emplace(0, aocommon::VectorMap<std::vector<std::complex<float>>>())
          .first->second;
  band_parm_response.AlwaysEmplace(data_desc_id) = parm_response;
  time_offsets_ = {std::pair(0, 0)};
}

void VisibilityModifier::InitializeCacheParmResponse(
    const std::vector<std::string>& antennaNames,
    const aocommon::MultiBandData& selected_bands, size_t ms_index) {
  using schaapcommon::h5parm::JonesParameters;

  const size_t solution_index = (*_h5parms).size() == 1 ? 0 : ms_index;

  // Only extract DD solutions if the corresponding cache entry is empty.
  aocommon::VectorMap<std::vector<std::complex<float>>>& parm_response =
      _cachedParmResponse[ms_index];
  if (parm_response.Empty()) {
    for (size_t data_desc_id : selected_bands.DataDescIds()) {
      const size_t nparms = NValuesPerSolution(ms_index);
      const aocommon::BandData& band = selected_bands[data_desc_id];
      const std::vector<double> freqs(band.begin(), band.end());
      const size_t response_size = _cachedMSTimes[ms_index]->size() *
                                   freqs.size() * antennaNames.size() * nparms;
      const std::string dir_name = (*_h5parms)[solution_index].GetNearestSource(
          _facetDirectionRA, _facetDirectionDec);
      schaapcommon::h5parm::SolTab* const first_solution =
          (*_firstSolutions)[solution_index];
      schaapcommon::h5parm::SolTab* const second_solution =
          _secondSolutions->empty() ? nullptr
                                    : (*_secondSolutions)[solution_index];
      const size_t dir_index = first_solution->GetDirIndex(dir_name);
      JonesParameters jones_parameters(
          freqs, *_cachedMSTimes[ms_index], antennaNames,
          (*_gainTypes)[solution_index],
          JonesParameters::InterpolationType::NEAREST, dir_index,
          first_solution, second_solution, false, 0u,
          JonesParameters::MissingAntennaBehavior::kUnit);
      // parms (Casacore::Cube) is column major
      const casacore::Cube<std::complex<float>>& parms =
          jones_parameters.GetParms();

      std::vector<std::complex<float>>& band_response =
          parm_response.AlwaysEmplace(data_desc_id);
      band_response.assign(&parms(0, 0, 0), &parms(0, 0, 0) + response_size);
      setNonFiniteToZero(band_response);
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
template <GainMode Mode>
void VisibilityModifier::ApplyBeamResponse(std::complex<float>* data,
                                           size_t n_channels,
                                           size_t data_desc_id, size_t antenna1,
                                           size_t antenna2) {
  for (size_t ch = 0; ch < n_channels; ++ch) {
    const size_t offset = ch * _pointResponseBufferSize;
    const size_t offset1 = offset + antenna1 * 4u;
    const size_t offset2 = offset + antenna2 * 4u;

    const MC2x2F gain1(&_cachedBeamResponse[data_desc_id][offset1]);
    const MC2x2F gain2(&_cachedBeamResponse[data_desc_id][offset2]);
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
