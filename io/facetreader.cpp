
#include "facetreader.h"

#include <schaapcommon/facets/ds9facetfile.h>

#include "../structures/facetutil.h"

using schaapcommon::facets::DS9FacetFile;
using schaapcommon::facets::Facet;

std::vector<std::shared_ptr<Facet>> FacetReader::ReadFacets(
    const Settings& settings, const ObservationInfo& observation_info) {
  const Facet::InitializationData data =
      CreateFacetInitializationData(settings, observation_info);

  std::vector<std::shared_ptr<Facet>> facets;
  if (!settings.facetRegionFilename.empty()) {
    facets = DS9FacetFile(settings.facetRegionFilename).ReadShared(data);

    if (facets.empty()) {
      throw std::runtime_error("No facets found in " +
                               settings.facetRegionFilename);
    }
  }

  return facets;
}

std::size_t FacetReader::CountFacets(const std::string& filename) {
  return filename.empty() ? 0 : DS9FacetFile(filename).Count();
}
