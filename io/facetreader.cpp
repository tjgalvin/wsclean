
#include "facetreader.h"

#include <schaapcommon/facets/ds9facetfile.h>

std::vector<schaapcommon::facets::Facet> FacetReader::ReadFacets(
    std::string facetRegionFilename) {
  std::vector<schaapcommon::facets::Facet> facets;
  if (!facetRegionFilename.empty()) {
    facets = schaapcommon::facets::DS9FacetFile(facetRegionFilename).Read();

    if (facets.empty())
      throw std::runtime_error("No facets found in " + facetRegionFilename);
  }

  return facets;
}
