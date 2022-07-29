#include "../io/cachedimageset.h"

#include <aocommon/fits/fitswriter.h>
#include <aocommon/image.h>
#include <aocommon/polarization.h>
#include <aocommon/uvector.h>
#include <schaapcommon/facets/facet.h>
#include <schaapcommon/facets/facetimage.h>

#include <boost/test/unit_test.hpp>

#include <math.h>
#include <limits>

using aocommon::PolarizationEnum;
using schaapcommon::facets::Facet;

BOOST_AUTO_TEST_SUITE(cachedimageset)

// NOTE: only check whether point x,y collides with the bounding box of a facet
size_t FacetCollision(const std::vector<std::shared_ptr<Facet>>& facets, int x,
                      int y) {
  size_t facetCollision = std::numeric_limits<size_t>::max();
  for (size_t f = 0; f < facets.size(); ++f) {
    schaapcommon::facets::BoundingBox bbox(facets[f]->GetPixels());
    if (y >= bbox.Min().y && y < bbox.Max().y && x >= bbox.Min().x &&
        x < bbox.Max().x) {
      facetCollision = f;
      break;
    }
  }
  return facetCollision;
}

BOOST_AUTO_TEST_CASE(store_and_load_facet) {
  std::string prefix = "facettest";
  size_t image_width = 8;
  size_t image_height = 8;
  double dl_dm = 0.0125;

  aocommon::FitsWriter writer;
  writer.SetImageDimensions(image_width, image_height, dl_dm, dl_dm);

  // Make two 4x4 facets
  std::vector<schaapcommon::facets::Coord> coords0{
      {0.05, -0.05}, {0.0, -0.05}, {0.0, 0.0}, {0.05, 0.0}};

  // Second facet (i=1) is mirrored in origin
  std::vector<schaapcommon::facets::Coord> coords1;
  for (schaapcommon::facets::Coord coordinate : coords0) {
    coords1.emplace_back(-1.0 * coordinate.ra, -1.0 * coordinate.dec);
  }

  Facet::InitializationData facet_data(writer.PixelSizeX(), writer.PixelSizeY(),
                                       writer.Width(), writer.Height());
  facet_data.phase_centre.ra = writer.RA();
  facet_data.phase_centre.dec = writer.Dec();
  facet_data.shift_l = writer.PhaseCentreDL();
  facet_data.shift_m = writer.PhaseCentreDM();
  // The bounding box is padded such that it is partially outside the main image
  facet_data.padding = 1.5;

  std::vector<std::shared_ptr<Facet>> facets{
      std::make_shared<Facet>(facet_data, std::move(coords0)),
      std::make_shared<Facet>(facet_data, std::move(coords1))};
  std::vector<aocommon::Image> facet_images;
  facet_images.reserve(facets.size());

  for (size_t i = 0; i < facets.size(); ++i) {
    facet_images.emplace_back(facets[i]->GetTrimmedBoundingBox().Width(),
                              facets[i]->GetTrimmedBoundingBox().Height(),
                              static_cast<float>(i + 1));
  }

  CachedImageSet cSet;
  cSet.Initialize(writer, 2, 1, facets.size(), prefix);

  std::vector<aocommon::PolarizationEnum> polarizations{
      aocommon::Polarization::XX, aocommon::Polarization::YY};

  for (const auto& polarization : polarizations) {
    for (size_t facet_idx = 0; facet_idx < facets.size(); ++facet_idx) {
      cSet.StoreFacet(facet_images[facet_idx], polarization, 1, facet_idx,
                      facets[facet_idx], false);
    }
  }

  // Retrieve the cached tmp filenames
  std::vector<std::string> storedNames(cSet.GetStoredNames().begin(),
                                       cSet.GetStoredNames().end());

  schaapcommon::facets::FacetImage imageStorage(image_width, image_height, 1);
  for (size_t pol_idx = 0; pol_idx < polarizations.size(); ++pol_idx) {
    aocommon::Image imageMain(image_width, image_height, 0.0f);
    for (size_t facet_idx = 0; facet_idx < facets.size(); ++facet_idx) {
      // Offset in file list
      size_t offset = pol_idx * facets.size() + facet_idx;
      imageStorage.SetFacet(*facets[facet_idx], true);
      BOOST_CHECK_EQUAL(storedNames[offset],
                        prefix + "-" +
                            aocommon::Polarization::TypeToShortString(
                                polarizations[pol_idx]) +
                            "-f000" + std::to_string(facet_idx) + "-tmp.fits");

      size_t num_facet_pixels =
          facets[facet_idx]->GetTrimmedBoundingBox().Width() *
          facets[facet_idx]->GetTrimmedBoundingBox().Height();
      cSet.LoadFacet(imageStorage.Data(0), polarizations[pol_idx], 1, facet_idx,
                     facets[facet_idx], false);
      imageStorage.AddToImage({imageMain.Data()});
      BOOST_CHECK_EQUAL_COLLECTIONS(
          imageStorage.Data(0), imageStorage.Data(0) + num_facet_pixels,
          facet_images[facet_idx].begin(), facet_images[facet_idx].end());
    }

    // Check whether data in imageMain is correct
    for (size_t y = 0; y != image_height; ++y) {
      for (size_t x = 0; x != image_width; ++x) {
        size_t offset = y * image_width + x;
        size_t facetCollision =
            FacetCollision(facets, static_cast<int>(x), static_cast<int>(y));
        BOOST_CHECK_EQUAL(imageMain[offset],
                          (facetCollision == std::numeric_limits<size_t>::max())
                              ? 0.0f
                              : facetCollision + 1.0f);
      }
    }
    cSet.Store(imageMain.Data(), polarizations[pol_idx], 1, false);
    BOOST_CHECK(
        cSet.GetStoredNames().find(
            prefix + "-" +
            aocommon::Polarization::TypeToShortString(polarizations[pol_idx]) +
            "-tmp.fits") != cSet.GetStoredNames().end());
  }
}

BOOST_AUTO_TEST_SUITE_END()