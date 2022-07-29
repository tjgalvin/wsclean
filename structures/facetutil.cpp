#include "facetutil.h"

#include <aocommon/imagecoordinates.h>

using schaapcommon::facets::Facet;

Facet::InitializationData CreateFacetInitializationData(
    double width, double height, double pixelScaleX, double pixelScaleY,
    double phaseCentreRA, double phaseCentreDec, double shiftL, double shiftM,
    double imagePadding, bool make_square) {
  Facet::InitializationData data(pixelScaleX, pixelScaleY, width, height);
  data.phase_centre.ra = phaseCentreRA;
  data.phase_centre.dec = phaseCentreDec;
  data.shift_l = shiftL;
  data.shift_m = shiftM;
  data.padding = imagePadding;
  data.align = 2;
  data.make_square = make_square;
  return data;
}

std::vector<std::shared_ptr<Facet>> CreateFacetGrid(
    const Facet::InitializationData& facet_data, size_t grid_width,
    size_t grid_height) {
  std::vector<std::shared_ptr<Facet>> facets;
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

      facets.emplace_back(std::make_shared<Facet>(facet_data, box));
      // add a name label for this box
      facets.back()->SetDirectionLabel(std::to_string(grid_x) + ", " +
                                       std::to_string(grid_y));
    }
  }

  return facets;
}