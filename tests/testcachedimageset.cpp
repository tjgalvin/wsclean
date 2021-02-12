#include "../io/cachedimageset.h"

#include <aocommon/polarization.h>
#include <aocommon/uvector.h>
#include <schaapcommon/facets/facet.h>

#include <boost/test/unit_test.hpp>

#include <math.h>

using aocommon::PolarizationEnum;

BOOST_AUTO_TEST_SUITE(cachedimageset)

BOOST_AUTO_TEST_CASE(store_and_load_facet) {
  std::string prefix = "facettest";
  size_t image_width = 8;
  size_t image_height = 8;
  double dl_dm = 0.0125;

  FitsWriter writer;
  writer.SetImageDimensions(image_width, image_height, dl_dm, dl_dm);

  // Make two 2x2 facets
  std::vector<std::pair<double, double>> coords{
      {0.05, -0.05}, {0.0, -0.05}, {0.0, 0.0}, {0.05, 0.0}};

  // Do not change num_facets!
  const size_t num_facets = 2;
  std::vector<schaapcommon::facets::Facet> facets(num_facets);
  std::vector<aocommon::UVector<float>> facets_data(num_facets);

  for (size_t i = 0; i < num_facets; ++i) {
    for (const auto& coord : coords) {
      // Second facet (i=1) is mirrored in origin
      facets[i].AddVertex(std::pow(-1, i) * coord.first,
                          std::pow(-1, i) * coord.second);
    }
    // dl (and dm) should be retrievable from FitsWriter
    facets[i].CalculatePixelPositions(
        writer.RA(), writer.Dec(), dl_dm, dl_dm, writer.Width(),
        writer.Height(), writer.PhaseCentreDL(), writer.PhaseCentreDM(), false);
    facets_data[i].assign(facets[i].GetBoundingBox().Width() *
                              facets[i].GetBoundingBox().Height(),
                          static_cast<float>(i + 1));
  }

  CachedImageSet cSet;
  cSet.Initialize(writer, 2, 1, facets.size(), prefix);

  std::vector<aocommon::PolarizationEnum> polarizations{
      aocommon::Polarization::XX, aocommon::Polarization::YY};

  for (auto& polarization : polarizations) {
    for (size_t facet_idx = 0; facet_idx < facets.size(); ++facet_idx) {
      cSet.StoreFacet(facets_data[facet_idx].data(), polarization, 1, facet_idx,
                      &facets[facet_idx], false);
    }
  }

  // Retrieve the cached filenames
  std::vector<std::string> storedNames(cSet.GetStoredNames().begin(),
                                       cSet.GetStoredNames().end());
  for (size_t pol_idx = 0; pol_idx < polarizations.size(); ++pol_idx) {
    for (size_t facet_idx = 0; facet_idx < facets.size(); ++facet_idx) {
      size_t offset = pol_idx * facets.size() + facet_idx;
      BOOST_CHECK_EQUAL(storedNames[offset],
                        prefix + "-" +
                            aocommon::Polarization::TypeToShortString(
                                polarizations[pol_idx]) +
                            "-f000" + std::to_string(facet_idx) + "-tmp.fits");

      size_t num_facet_pixels = facets[facet_idx].GetBoundingBox().Width() *
                                facets[facet_idx].GetBoundingBox().Height();
      aocommon::UVector<float> read_buffer(num_facet_pixels);
      cSet.LoadFacet(read_buffer.data(), polarizations[pol_idx], 1, facet_idx,
                     &facets[facet_idx], false);
      BOOST_CHECK_EQUAL_COLLECTIONS(
          read_buffer.data(), read_buffer.data() + num_facet_pixels,
          facets_data[facet_idx].data(),
          facets_data[facet_idx].data() + num_facet_pixels);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()