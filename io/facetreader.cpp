
#include "facetreader.h"

#include <schaapcommon/facets/ds9facetfile.h>

using schaapcommon::facets::DS9FacetFile;
using schaapcommon::facets::Facet;

std::vector<std::shared_ptr<Facet>> FacetReader::ReadFacets(
    const Settings& settings, const ObservationInfo& observation_info) {
  Facet::InitializationData data(settings.pixelScaleX, settings.pixelScaleY,
                                 settings.trimmedImageWidth,
                                 settings.trimmedImageHeight);
  data.phase_centre.ra = observation_info.phaseCentreRA;
  data.phase_centre.dec = observation_info.phaseCentreDec;
  data.shift_l = observation_info.shiftL;
  data.shift_m = observation_info.shiftM;
  data.padding = settings.imagePadding;
  data.align = 2;
  data.make_square = settings.gridderType == GridderType::IDG;

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
