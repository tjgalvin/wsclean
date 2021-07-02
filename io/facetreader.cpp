
#include "facetreader.h"

#include <schaapcommon/facets/ds9facetfile.h>

std::vector<std::shared_ptr<schaapcommon::facets::Facet>>
FacetReader::ReadFacets(std::string facetRegionFilename) {
  std::vector<std::shared_ptr<schaapcommon::facets::Facet>> facets;
  if (!facetRegionFilename.empty()) {
    facets =
        schaapcommon::facets::DS9FacetFile(facetRegionFilename).ReadShared();

    if (facets.empty())
      throw std::runtime_error("No facets found in " + facetRegionFilename);
  }

  return facets;
}
