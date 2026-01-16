#include "msgridderdata.h"

#include "msprovidercollection.h"

namespace wsclean {

MsGridderData::MsGridderData(const Settings& settings)
    : settings_(settings),
      visibility_weighting_mode_(settings.visibilityWeightingMode) {
#ifdef HAVE_EVERYBEAM
  visibility_modifier_.SetBeamInfo(settings.beamMode,
                                   settings.beamNormalisationMode);
#endif
}

void MsGridderData::StartMeasurementSet(
    size_t ms_count, const MsProviderCollection::MsData& ms_data,
    bool is_predict) {
  visibility_modifier_.InitializeBda(ms_data.ms_provider->HasFrequencyBda());
  InitializePointResponse(ms_data);
  if (settings_.applyTimeFrequencySmearing) {
    visibility_modifier_.InitializeTimeFrequencySmearing(
        ms_data.ms_provider->GetObservationInfo(),
        ms_data.ms_provider->Interval());
  }
  if (visibility_modifier_.HasH5Parm()) {
    visibility_modifier_.InitializeCacheParmResponse(
        ms_data.antenna_names, ms_data.ms_provider->SelectedBands(),
        ms_data.original_ms_index);
  }
  original_ms_index_ = ms_data.original_ms_index;
  writer_lock_index_ =
      facet_group_index_ * ms_count + ms_data.original_ms_index;
  n_vis_polarizations_ = ms_data.ms_provider->NPolarizations();
  gain_mode_ = SelectGainMode(Polarization(), n_vis_polarizations_);
  const size_t n_channels = ms_data.ms_provider->NMaxChannels();
  scratch_image_weights_.resize(n_channels);
  if (is_predict) {
    scratch_model_data_.resize(n_channels *
                               ms_data.ms_provider->NPolarizations());
    predict_reader_ = ms_data.ms_provider->MakeReader();
  }
}

void MsGridderData::ReadPredictMetaData(MSProvider::MetaData& metadata) {
  predict_reader_->ReadMeta(metadata);
  predict_reader_->NextInputRow();
}

void MsGridderData::ResetVisibilityModifierCache(size_t ms_count) {
  visibility_modifier_.ResetCache(ms_count);
  visibility_modifier_.ResetSums();
}

void MsGridderData::InitializePointResponse(
    const MsProviderCollection::MsData& ms_data) {
#ifdef HAVE_EVERYBEAM
  if (settings_.applyFacetBeam || settings_.gridWithBeam) {
    const std::string element_response_string =
        !settings_.beamModel.empty() ? settings_.beamModel : "DEFAULT";
    std::optional<SynchronizedMS> optional_ms =
        ms_data.ms_provider->MsIfAvailable();
    if (!optional_ms) {
      throw std::runtime_error(
          "Can not compute the beam without an on-disk measurement set");
    }
    visibility_modifier_.InitializePointResponse(
        std::move(*optional_ms), settings_.facetBeamUpdateTime,
        element_response_string, ms_data.ms_provider->SelectedBands(),
        settings_.dataColumnName, settings_.mwaPath);
  } else {
    visibility_modifier_.SetNoPointResponse();
  }
#else
  if (settings_.applyFacetBeam || settings_.gridWithBeam) {
    throw std::runtime_error(
        "-apply-facet-beam or -grid-with-beam was set, but wsclean was not "
        "compiled "
        "with EveryBeam. Please compile wsclean with EveryBeam to "
        "use the Facet Beam functionality");
  }
#endif
}

}  // namespace wsclean
