#include "facetutil.h"

#include <aocommon/imagecoordinates.h>

using schaapcommon::facets::Facet;

Facet::InitializationData CreateFacetInitializationData(
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
  return data;
}

std::vector<Facet> CreateFacetGrid(const Facet::InitializationData& facet_data,
                                   size_t grid_width, size_t grid_height) {
  std::vector<Facet> facets;
  facets.reserve(grid_height * grid_width);

  for (int grid_y = 0; grid_y < static_cast<int>(grid_height); ++grid_y) {
    const int facet_start_y = grid_y * facet_data.image_height / grid_height;
    const int facet_end_y =
        (grid_y + 1) * facet_data.image_height / grid_height;

    for (int grid_x = 0; grid_x < static_cast<int>(grid_width); ++grid_x) {
      const int facet_start_x = grid_x * facet_data.image_width / grid_width;
      const int facet_end_x =
          (grid_x + 1) * facet_data.image_width / grid_width;
      const schaapcommon::facets::BoundingBox box(
          {{facet_start_x, facet_start_y}, {facet_end_x, facet_end_y}});

      facets.emplace_back(facet_data, box);
      // add a name label for this box
      facets.back().SetDirectionLabel(std::to_string(grid_x) + ", " +
                                      std::to_string(grid_y));
    }
  }

  return facets;
}