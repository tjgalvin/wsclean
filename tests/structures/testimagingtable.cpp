#include "../../structures/imagingtable.h"

#include <radler/work_table.h>

#include <array>

#include <boost/test/unit_test.hpp>

#include "../../io/cachedimageaccessor.h"
#include "../../io/cachedimageset.h"

#include "../common/smartptr.h"

namespace {

void TestGroupCounts(const ImagingTable& table, size_t indepentGroupCount,
                     size_t facetGroupCount, size_t facetCount,
                     size_t squaredGroupCount) {
  BOOST_TEST(table.IndependentGroupCount() == indepentGroupCount);
  BOOST_TEST(table.FacetGroupCount() == facetGroupCount);
  BOOST_TEST(table.FacetGroups().size() == facetGroupCount);
  BOOST_TEST(table.FacetCount() == facetCount);
  BOOST_TEST(table.Facets().size() == facetCount);
  BOOST_TEST(table.SquaredGroupCount() == squaredGroupCount);
  BOOST_TEST(table.SquaredGroups().size() == squaredGroupCount);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(imagingtable)

BOOST_AUTO_TEST_CASE(constructor) {
  ImagingTable table;
  TestGroupCounts(table, 0, 0, 0, 0);
  BOOST_TEST(table.EntryCount() == 0u);
  BOOST_TEST((table.begin() == table.end()));
}

BOOST_AUTO_TEST_CASE(add_update_clear) {
  ImagingTable table;
  test::UniquePtr<ImagingTableEntry> entry1;
  test::UniquePtr<ImagingTableEntry> entry2;
  entry1->joinedGroupIndex = 42;
  entry2->joinedGroupIndex = 43;
  entry1->facetGroupIndex = 142;
  entry2->facetGroupIndex = 142;
  entry1->squaredDeconvolutionIndex = 242;
  entry2->squaredDeconvolutionIndex = 242;
  table.AddEntry(entry1.take());
  table.AddEntry(entry2.take());

  BOOST_TEST_REQUIRE(table.EntryCount() == 2u);
  BOOST_TEST(&table[0] == entry1.get());
  BOOST_TEST(&table[1] == entry2.get());
  BOOST_TEST(&table.Front() == entry1.get());

  // Test begin() and end() functions, including iterators.
  auto iterator = table.begin();
  BOOST_TEST_REQUIRE((iterator != table.end()));
  BOOST_TEST(&*iterator == entry1.get());
  ++iterator;
  BOOST_TEST_REQUIRE((iterator != table.end()));
  BOOST_TEST(&*iterator == entry2.get());
  ++iterator;
  BOOST_TEST((iterator == table.end()));

  // No Update called -> only EntryCount is correct.
  TestGroupCounts(table, 0, 0, 0, 0);

  table.Update();
  TestGroupCounts(table, 2, 1, 1, 1);

  table.Clear();
  BOOST_TEST(table.EntryCount() == 0u);
  TestGroupCounts(table, 0, 0, 0, 0);
}

BOOST_AUTO_TEST_CASE(independent_groups) {
  ImagingTable table;

  test::UniquePtr<ImagingTableEntry> entry0_0;
  test::UniquePtr<ImagingTableEntry> entry0_1;
  test::UniquePtr<ImagingTableEntry> entry1_0;
  test::UniquePtr<ImagingTableEntry> entry1_1;
  entry0_0->joinedGroupIndex = 42;
  entry0_1->joinedGroupIndex = 42;
  entry1_0->joinedGroupIndex = 43;
  entry1_1->joinedGroupIndex = 43;

  // Deliberately mix the groups when adding the entries:
  table.AddEntry(entry0_0.take());
  table.AddEntry(entry1_0.take());
  table.AddEntry(entry1_1.take());
  table.AddEntry(entry0_1.take());

  BOOST_TEST(table.IndependentGroupCount() == 0u);
  table.Update();
  BOOST_TEST_REQUIRE(table.IndependentGroupCount() == 2u);
  TestGroupCounts(table, 2, 1, 1, 1);

  ImagingTable group0 = table.GetIndependentGroup(0);
  ImagingTable group1 = table.GetIndependentGroup(1);
  BOOST_TEST_REQUIRE(group0.EntryCount() == 2u);
  BOOST_TEST_REQUIRE(group1.EntryCount() == 2u);
  BOOST_TEST(&group0[0] == entry0_0.get());
  BOOST_TEST(&group0[1] == entry0_1.get());
  BOOST_TEST(&group1[0] == entry1_0.get());
  BOOST_TEST(&group1[1] == entry1_1.get());

  // Test that Update() was called on group0 and group1.
  TestGroupCounts(group0, 1, 1, 1, 1);
  TestGroupCounts(group1, 1, 1, 1, 1);
}

BOOST_AUTO_TEST_CASE(facet_groups) {
  ImagingTable table;

  test::UniquePtr<ImagingTableEntry> entry0_0;
  test::UniquePtr<ImagingTableEntry> entry0_1;
  test::UniquePtr<ImagingTableEntry> entry0_2;
  test::UniquePtr<ImagingTableEntry> entry1_0;
  entry0_0->facetGroupIndex = 43;
  entry0_1->facetGroupIndex = 43;
  entry0_2->facetGroupIndex = 43;
  entry1_0->facetGroupIndex = 42;
  table.AddEntry(entry0_0.take());
  table.AddEntry(entry0_1.take());
  table.AddEntry(entry0_2.take());
  table.AddEntry(entry1_0.take());

  BOOST_TEST(table.FacetGroupCount() == 0u);
  table.Update();
  BOOST_TEST_REQUIRE(table.FacetGroupCount() == 2u);
  TestGroupCounts(table, 1, 2, 1, 1);

  // Test GetFacetGroup. The ImagingTable orders the groups by group index.
  ImagingTable group0 = table.GetFacetGroup(0);
  ImagingTable group1 = table.GetFacetGroup(1);
  BOOST_TEST_REQUIRE(group1.EntryCount() == 3u);
  BOOST_TEST_REQUIRE(group0.EntryCount() == 1u);
  BOOST_TEST(&group1[0] == entry0_0.get());
  BOOST_TEST(&group1[1] == entry0_1.get());
  BOOST_TEST(&group1[2] == entry0_2.get());
  BOOST_TEST(&group0[0] == entry1_0.get());
  // Test that Update() was called on group0 and group1.
  TestGroupCounts(group0, 1, 1, 1, 1);
  TestGroupCounts(group1, 1, 1, 1, 1);

  // Test FacetGroups. The ImagingTable orders the groups by group index.
  const ImagingTable::Groups& groups = table.FacetGroups();
  BOOST_TEST_REQUIRE(groups.size() == 2u);
  BOOST_TEST_REQUIRE(groups[1].size() == 3u);
  BOOST_TEST_REQUIRE(groups[0].size() == 1u);
  BOOST_TEST(groups[1][0].get() == entry0_0.get());
  BOOST_TEST(groups[1][1].get() == entry0_1.get());
  BOOST_TEST(groups[1][2].get() == entry0_2.get());
  BOOST_TEST(groups[0][0].get() == entry1_0.get());
}

BOOST_AUTO_TEST_CASE(facets) {
  ImagingTable table;

  test::UniquePtr<ImagingTableEntry> entry0_0;
  test::UniquePtr<ImagingTableEntry> entry0_1;
  test::UniquePtr<ImagingTableEntry> entry0_2;
  test::UniquePtr<ImagingTableEntry> entry1_0;
  entry0_0->facetIndex = 43;
  entry0_1->facetIndex = 43;
  entry0_2->facetIndex = 43;
  entry1_0->facetIndex = 42;
  table.AddEntry(entry0_0.take());
  table.AddEntry(entry0_1.take());
  table.AddEntry(entry0_2.take());
  table.AddEntry(entry1_0.take());

  BOOST_TEST(table.FacetCount() == 0u);
  table.Update();
  BOOST_TEST_REQUIRE(table.FacetCount() == 2u);
  TestGroupCounts(table, 1, 1, 2, 1);

  // Test GetFacet. The ImagingTable orders the facets by facet index.
  ImagingTable group0 = table.GetFacet(0);
  ImagingTable group1 = table.GetFacet(1);
  BOOST_TEST_REQUIRE(group1.EntryCount() == 3u);
  BOOST_TEST_REQUIRE(group0.EntryCount() == 1u);
  BOOST_TEST(&group1[0] == entry0_0.get());
  BOOST_TEST(&group1[1] == entry0_1.get());
  BOOST_TEST(&group1[2] == entry0_2.get());
  BOOST_TEST(&group0[0] == entry1_0.get());
  // Test that Update() was called on group0 and group1.
  TestGroupCounts(group0, 1, 1, 1, 1);
  TestGroupCounts(group1, 1, 1, 1, 1);

  // Test Facets. The ImagingTable orders the facets by facet index.
  const ImagingTable::Groups& groups = table.Facets();
  BOOST_TEST_REQUIRE(groups.size() == 2u);
  BOOST_TEST_REQUIRE(groups[1].size() == 3u);
  BOOST_TEST_REQUIRE(groups[0].size() == 1u);
  BOOST_TEST(groups[1][0].get() == entry0_0.get());
  BOOST_TEST(groups[1][1].get() == entry0_1.get());
  BOOST_TEST(groups[1][2].get() == entry0_2.get());
  BOOST_TEST(groups[0][0].get() == entry1_0.get());
}

BOOST_AUTO_TEST_CASE(squared_groups) {
  ImagingTable table;

  test::UniquePtr<ImagingTableEntry> entry0_0;
  test::UniquePtr<ImagingTableEntry> entry0_1;
  test::UniquePtr<ImagingTableEntry> entry0_2;
  test::UniquePtr<ImagingTableEntry> entry1_0;
  entry0_0->squaredDeconvolutionIndex = 42;
  entry0_1->squaredDeconvolutionIndex = 42;
  entry0_2->squaredDeconvolutionIndex = 42;
  entry1_0->squaredDeconvolutionIndex = 43;

  // Deliberately mix the groups when adding the entries:
  table.AddEntry(entry0_0.take());
  table.AddEntry(entry1_0.take());
  table.AddEntry(entry0_1.take());
  table.AddEntry(entry0_2.take());

  BOOST_TEST(table.SquaredGroups().size() == 0u);
  table.Update();
  BOOST_TEST_REQUIRE(table.SquaredGroups().size() == 2u);
  TestGroupCounts(table, 1, 1, 1, 2);

  // Test GetSquaredGroup
  ImagingTable group0 = table.GetSquaredGroup(0);
  ImagingTable group1 = table.GetSquaredGroup(1);
  BOOST_TEST_REQUIRE(group0.EntryCount() == 3u);
  BOOST_TEST_REQUIRE(group1.EntryCount() == 1u);
  BOOST_TEST(&group0[0] == entry0_0.get());
  BOOST_TEST(&group0[1] == entry0_1.get());
  BOOST_TEST(&group0[2] == entry0_2.get());
  BOOST_TEST(&group1[0] == entry1_0.get());
  // Test that Update() was called on group0 and group1.
  TestGroupCounts(group0, 1, 1, 1, 1);
  TestGroupCounts(group1, 1, 1, 1, 1);

  // Test SquaredGroups
  const ImagingTable::Groups& groups = table.SquaredGroups();
  BOOST_TEST_REQUIRE(groups.size() == 2u);
  BOOST_TEST_REQUIRE(groups[0].size() == 3u);
  BOOST_TEST_REQUIRE(groups[1].size() == 1u);
  BOOST_TEST(groups[0][0].get() == entry0_0.get());
  BOOST_TEST(groups[0][1].get() == entry0_1.get());
  BOOST_TEST(groups[0][2].get() == entry0_2.get());
  BOOST_TEST(groups[1][0].get() == entry1_0.get());
}

BOOST_AUTO_TEST_CASE(create_deconvolution_table_single_entry) {
  ImagingTable table;
  CachedImageSet psf_images;
  CachedImageSet model_images;
  CachedImageSet residual_images;

  // Check for exception on empty table
  BOOST_CHECK_THROW(std::unique_ptr<radler::WorkTable> deconvolution_table =
                        table.CreateDeconvolutionTable(
                            1, psf_images, model_images, residual_images),
                    std::runtime_error);

  test::UniquePtr<ImagingTableEntry> entry;
  entry->polarization = aocommon::PolarizationEnum::StokesV;
  entry->bandStartFrequency = 42.0;
  entry->bandEndFrequency = 43.0;
  entry->outputChannelIndex = 0;
  entry->outputIntervalIndex = 15;
  entry->imageWeight = 44.0;
  entry->imageCount = 2;
  table.AddEntry(entry.take());

  // Check for invalid table
  // Update was not called after adding entry
  BOOST_CHECK_THROW(std::unique_ptr<radler::WorkTable> deconvolution_table =
                        table.CreateDeconvolutionTable(
                            1, psf_images, model_images, residual_images),
                    std::runtime_error);
  table.Update();

  std::unique_ptr<radler::WorkTable> deconvolution_table =
      table.CreateDeconvolutionTable(1, psf_images, model_images,
                                     residual_images);
  BOOST_TEST_REQUIRE(deconvolution_table->OriginalGroups().size() == 1u);
  BOOST_TEST_REQUIRE(deconvolution_table->OriginalGroups()[0].size() == 2u);

  const std::vector<std::vector<size_t>>& deconvolution_groups =
      deconvolution_table->DeconvolutionGroups();
  BOOST_TEST_REQUIRE(deconvolution_groups.size() == 1);
  BOOST_TEST_REQUIRE(deconvolution_groups.front() == std::vector<size_t>(1, 0));

  BOOST_TEST_REQUIRE(deconvolution_table->Size() == entry->imageCount);

  size_t index = 0;
  for (const radler::WorkTableEntry& deconvolution_entry :
       *deconvolution_table) {
    BOOST_TEST(&deconvolution_entry ==
               deconvolution_table->OriginalGroups()[0][index]);
    BOOST_TEST(deconvolution_entry.index == index);
    BOOST_TEST(deconvolution_entry.band_start_frequency ==
               entry->bandStartFrequency);
    BOOST_TEST(deconvolution_entry.band_end_frequency ==
               entry->bandEndFrequency);
    BOOST_TEST(deconvolution_entry.polarization == entry->polarization);
    BOOST_TEST(deconvolution_entry.original_channel_index ==
               entry->outputChannelIndex);
    BOOST_TEST(deconvolution_entry.original_interval_index ==
               entry->outputIntervalIndex);
    BOOST_TEST(deconvolution_entry.image_weight == entry->imageWeight);

    if (index == 0) {
      BOOST_TEST_REQUIRE(deconvolution_entry.psf_accessors.size() == 1);
      auto* psf_accessor = dynamic_cast<CachedImageAccessor*>(
          deconvolution_entry.psf_accessors.front().get());
      BOOST_TEST_REQUIRE(psf_accessor);
      BOOST_TEST(&psf_accessor->GetImageSet() == &psf_images);
      BOOST_TEST(psf_accessor->GetPolarization() == entry->polarization);
      BOOST_TEST(psf_accessor->GetFrequencyIndex() ==
                 entry->outputChannelIndex);
      BOOST_TEST(!psf_accessor->GetIsImaginary());
    } else {
      BOOST_TEST(deconvolution_entry.psf_accessors.empty());
    }

    auto* model_accessor = dynamic_cast<CachedImageAccessor*>(
        deconvolution_entry.model_accessor.get());
    BOOST_TEST_REQUIRE(model_accessor);
    BOOST_TEST(&model_accessor->GetImageSet() == &model_images);
    BOOST_TEST(model_accessor->GetPolarization() == entry->polarization);
    BOOST_TEST(model_accessor->GetFrequencyIndex() ==
               entry->outputChannelIndex);
    BOOST_TEST(model_accessor->GetIsImaginary() == (index == 1));

    auto* residual_accessor = dynamic_cast<CachedImageAccessor*>(
        deconvolution_entry.residual_accessor.get());
    BOOST_TEST_REQUIRE(residual_accessor);
    BOOST_TEST(&residual_accessor->GetImageSet() == &residual_images);
    BOOST_TEST(residual_accessor->GetPolarization() == entry->polarization);
    BOOST_TEST(residual_accessor->GetFrequencyIndex() ==
               entry->outputChannelIndex);
    BOOST_TEST(residual_accessor->GetIsImaginary() == (index == 1));

    ++index;
  }
}

BOOST_AUTO_TEST_CASE(create_deconvolution_table_multiple_groups) {
  ImagingTable table;
  std::array<test::UniquePtr<ImagingTableEntry>, 5> entries;
  for (size_t i = 0; i < entries.size(); ++i) {
    entries[i]->outputChannelIndex = i;
    entries[i]->imageWeight = 42 + i;
    entries[i]->imageCount = 1;
    table.AddEntry(entries[i].take());
  }
  table.Update();

  CachedImageSet images;
  std::unique_ptr<radler::WorkTable> deconvolution_table =
      table.CreateDeconvolutionTable(entries.size(), images, images, images);

  BOOST_TEST_REQUIRE(deconvolution_table->Size() == entries.size());
  BOOST_TEST_REQUIRE(deconvolution_table->OriginalGroups().size() ==
                     entries.size());
  BOOST_TEST_REQUIRE(deconvolution_table->DeconvolutionGroups().size() ==
                     entries.size());

  auto entry_iterator = deconvolution_table->Begin();
  for (int index = 0; index < int(entries.size()); ++index, ++entry_iterator) {
    BOOST_TEST((entry_iterator != deconvolution_table->End()));

    const radler::WorkTableEntry& deconvolution_entry = *entry_iterator;
    BOOST_TEST(deconvolution_entry.original_channel_index ==
               entries[index]->outputChannelIndex);
    // Only use image_weight for checking if the entry is the correct entry.
    // The single entry test above checks if the other fields are copied.
    BOOST_TEST(deconvolution_entry.image_weight == entries[index]->imageWeight);

    const radler::WorkTable::Group& group =
        deconvolution_table->OriginalGroups()[index];
    BOOST_TEST_REQUIRE(group.size() == 1);
    BOOST_TEST(group[0] == &deconvolution_entry);

    BOOST_TEST(deconvolution_table->DeconvolutionGroups()[index] ==
               std::vector<size_t>(1, index));
  }
  BOOST_TEST((entry_iterator == deconvolution_table->End()));
}

BOOST_AUTO_TEST_SUITE_END()
