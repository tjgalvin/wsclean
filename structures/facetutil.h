#ifndef STRUCTURES_DDPSF_H_
#define STRUCTURES_DDPSF_H_

#include <vector>

#include <schaapcommon/facets/facet.h>

#include "../main/settings.h"
#include "../structures/observationinfo.h"

schaapcommon::facets::Facet::InitializationData CreateFacetInitializationData(
    const Settings& settings, const ObservationInfo& observation_info);

std::vector<schaapcommon::facets::Facet> CreateFacetGrid(
    const schaapcommon::facets::Facet::InitializationData& facet_data,
    size_t grid_width, size_t grid_height);

#endif  // STRUCTURES_DDPSF_H_