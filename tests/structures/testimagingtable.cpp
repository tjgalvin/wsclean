#include <boost/test/unit_test.hpp>

#include "../../structures/imagingtable.h"

#include "../common/smartptr.h"

namespace {

void TestGroupCounts(const ImagingTable& table, size_t indepentGroupCount,
                     size_t facetGroupCount, size_t squaredGroupCount) {
  BOOST_TEST(table.IndependentGroupCount() == indepentGroupCount);
  BOOST_TEST(table.FacetGroupCount() == facetGroupCount);
  BOOST_TEST(table.FacetGroups().size() == facetGroupCount);
  BOOST_TEST(table.SquaredGroupCount() == squaredGroupCount);
  BOOST_TEST(table.SquaredGroups().size() == squaredGroupCount);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(imagingtable)

BOOST_AUTO_TEST_CASE(constructor) {
  ImagingTable table;
  TestGroupCounts(table, 0, 0, 0);
  BOOST_TEST(table.EntryCount() == 0u);
  BOOST_TEST((table.begin() == table.end()));
}

BOOST_AUTO_TEST_CASE(add_update_clear) {
  ImagingTable table;
  test::UniquePtr<ImagingTableEntry> entry1;
  test::UniquePtr<ImagingTableEntry> entry2;
  entry1->joinedGroupIndex = 42;
  entry2->joinedGroupIndex = 43;
  entry1->facetIndex = 142;
  entry2->facetIndex = 142;
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
  TestGroupCounts(table, 0, 0, 0);

  table.Update();
  TestGroupCounts(table, 2, 1, 1);

  table.Clear();
  BOOST_TEST(table.EntryCount() == 0u);
  TestGroupCounts(table, 0, 0, 0);
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
  TestGroupCounts(table, 2, 1, 1);

  ImagingTable group0 = table.GetIndependentGroup(0);
  ImagingTable group1 = table.GetIndependentGroup(1);
  BOOST_TEST_REQUIRE(group0.EntryCount() == 2u);
  BOOST_TEST_REQUIRE(group1.EntryCount() == 2u);
  BOOST_TEST(&group0[0] == entry0_0.get());
  BOOST_TEST(&group0[1] == entry0_1.get());
  BOOST_TEST(&group1[0] == entry1_0.get());
  BOOST_TEST(&group1[1] == entry1_1.get());

  // Test that Update() was called on group0 and group1.
  TestGroupCounts(group0, 1, 1, 1);
  TestGroupCounts(group1, 1, 1, 1);
}

BOOST_AUTO_TEST_CASE(facet_groups) {
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

  BOOST_TEST(table.FacetGroupCount() == 0u);
  table.Update();
  BOOST_TEST_REQUIRE(table.FacetGroupCount() == 2u);
  TestGroupCounts(table, 1, 2, 1);

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
  TestGroupCounts(group0, 1, 1, 1);
  TestGroupCounts(group1, 1, 1, 1);

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
  TestGroupCounts(table, 1, 1, 2);

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
  TestGroupCounts(group0, 1, 1, 1);
  TestGroupCounts(group1, 1, 1, 1);

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

BOOST_AUTO_TEST_SUITE_END()
