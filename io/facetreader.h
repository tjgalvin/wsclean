#ifndef FACET_READER_H
#define FACET_READER_H

#include <schaapcommon/facets/facet.h>
#include "../main/settings.h"

class FacetReader {
 public:
  // Reading facets requires the scale and size so do it after those settings
  // are validated, and validate the facet settings here.
  static std::vector<schaapcommon::facets::Facet> ReadFacets(
      std::string facetRegionFilename);
};
#endif
