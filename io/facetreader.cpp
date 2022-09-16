#include "facetreader.h"

#include <algorithm>

#include <schaapcommon/facets/ds9facetfile.h>

#include <aocommon/logger.h>

#include "../structures/facetutil.h"

using schaapcommon::facets::DS9FacetFile;
using schaapcommon::facets::Facet;

std::vector<std::shared_ptr<Facet>> FacetReader::ReadFacets(
    std::string filename, double width, double height, double pixelScaleX,
    double pixelScaleY, double phaseCentreRA, double phaseCentreDec,
    double shiftL, double shiftM, double imagePadding, bool make_square) {
  const Facet::InitializationData data = CreateFacetInitializationData(
      width, height, pixelScaleX, pixelScaleY, phaseCentreRA, phaseCentreDec,
      shiftL, shiftM, imagePadding, make_square);

  std::vector<std::shared_ptr<Facet>> facets;
  if (!filename.empty()) {
    facets = DS9FacetFile(filename).ReadShared(data);

    if (facets.empty()) {
      throw std::runtime_error("No facets found in " + filename);
    }

    std::vector<std::shared_ptr<Facet>>::const_iterator to_remove =
        std::remove_if(
            facets.begin(), facets.end(),
            [](std::shared_ptr<Facet> facet) { return facet->Empty(); });
    if (to_remove != facets.end()) {
      const size_t n = facets.end() - to_remove;
      aocommon::Logger::Warn << n
                             << " facets fall outside the image boundaries.\n";
      facets.erase(to_remove, facets.end());
      if (facets.empty()) {
        throw std::runtime_error("All facets fall outside the image");
      }
    }
  }

  return facets;
}

std::size_t FacetReader::CountFacets(const std::string& filename) {
  return filename.empty() ? 0 : DS9FacetFile(filename).Count();
}
