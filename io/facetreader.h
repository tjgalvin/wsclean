#ifndef FACET_READER_H
#define FACET_READER_H

#include <memory>
#include <string>
#include <vector>

#include <schaapcommon/facets/facet.h>

#include "../main/settings.h"
#include "../structures/observationinfo.h"

class FacetReader {
 public:
  // Reading facets requires the scale and size so do it after those settings
  // are validated, and validate the facet settings here.
  static std::vector<std::shared_ptr<schaapcommon::facets::Facet>> ReadFacets(
      const Settings& settings, const ObservationInfo& observation_info);

  static std::size_t CountFacets(const std::string& facet_region_filename);
};
#endif
