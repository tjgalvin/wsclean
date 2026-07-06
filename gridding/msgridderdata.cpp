#include "msgridderdata.h"

#include "msprovidercollection.h"

namespace wsclean {

MsGridderData::MsGridderData(const Settings& settings)
    : settings_(settings), visibility_weighting_mode_(settings.visibilityWeightingMode) {
#ifdef HAVE_EVERYBEAM
  visibility_modifier_.SetBeamInfo(settings.beamMode, settings.beamNormalisationMode);
#endif
}

void MsGridderData::StartMeasurementSet(size_t ms_count,
                                        const MsProviderCollection::MsData& ms_data,
                                        bool is_predict) {
  visibility_modifier_.InitializeBda(ms_data.ms_provider->HasFrequencyBda());
  // Clear beam cache from previous measurement set to prevent stale beam data.
  visibility_modifier_.ResetBeamCache();
  InitializePointResponse(ms_data);
  if (settings_.applyTimeFrequencySmearing) {
    // For smearing calculations, the facet target direction is used, such that the target
    // correction is accurately corrected and sources away from the facet target direction
    // receive a slightly inaccurate correction. The target direction is typically the
    // direction toward the strongest source. It is not the same as the facet's centre
    // position (given by the facet's l/m shift), which is the centre of the bounding box
    // of the facet.
    const double facet_ra = visibility_modifier_.FacetDirectionRA();
    const double facet_dec = visibility_modifier_.FacetDirectionDec();
    const double phase_centre_ra = ms_data.ms_provider->GetObservationInfo().phaseCentreRA;
    const double phase_centre_dec = ms_data.ms_provider->GetObservationInfo().phaseCentreDec;
    double l;
    double m;
    aocommon::ImageCoordinates::RaDecToLM(facet_ra, facet_dec, phase_centre_ra, phase_centre_dec, l,
                                          m);
    visibility_modifier_.InitializeTimeFrequencySmearing(ms_data.ms_provider->GetObservationInfo(),
                                                         l, m);
  }
  if (visibility_modifier_.HasH5Parm()) {
    visibility_modifier_.InitializeCacheParmResponse(
        ms_data.antenna_names, ms_data.ms_provider->SelectedBands(), ms_data.original_ms_index);
  }
  original_ms_index_ = ms_data.original_ms_index;
  writer_lock_index_ = facet_group_index_ * ms_count + ms_data.original_ms_index;
  n_vis_polarizations_ = ms_data.ms_provider->NPolarizations();
  gain_mode_ = SelectGainMode(Polarization(), n_vis_polarizations_);
  const size_t n_channels = ms_data.ms_provider->NMaxChannels();
  scratch_image_weights_.resize(n_channels);
  if (is_predict) {
    scratch_model_data_.resize(n_channels * ms_data.ms_provider->NPolarizations());
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

void MsGridderData::InitializePointResponse(const MsProviderCollection::MsData& ms_data) {
#ifdef HAVE_EVERYBEAM
  if (settings_.applyFacetBeam || settings_.gridWithBeam) {
    const std::string element_response_string =
        !settings_.beamModel.empty() ? settings_.beamModel : "DEFAULT";
    std::optional<SynchronizedMS> optional_ms = ms_data.ms_provider->MsIfAvailable();
    if (!optional_ms) {
      throw std::runtime_error("Can not compute the beam without an on-disk measurement set");
    }
    visibility_modifier_.InitializePointResponse(
        std::move(*optional_ms), settings_.facetBeamUpdateTime, element_response_string,
        ms_data.ms_provider->SelectedBands(), settings_.dataColumnName, settings_.mwaPath);
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
