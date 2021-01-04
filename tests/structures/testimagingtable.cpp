#include <boost/test/unit_test.hpp>

#include "../../structures/imagingtable.h"

#include "../common/smartptr.h"

namespace {

void TestGroupCounts(const ImagingTable& table, size_t indepentGroupCount,
                     size_t facetGroupCount, size_t squaredGroupCount) {
  BOOST_TEST(table.IndependentGroupCount() == indepentGroupCount);
  BOOST_TEST(table.FacetGroupCount() == facetGroupCount);
  BOOST_TEST(table.SquaredGroupCount() == squaredGroupCount);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(imagingtable)

BOOST_AUTO_TEST_CASE(constructor) {
  ImagingTable table;
  TestGroupCounts(table, 0, 0, 0);
  BOOST_TEST(table.EntryCount() == 0u);
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
  TestGroupCounts(table, 1, 2, 1);

  ImagingTable group0 = table.GetFacetGroup(0);
  ImagingTable group1 = table.GetFacetGroup(1);
  BOOST_TEST_REQUIRE(group0.EntryCount() == 3u);
  BOOST_TEST_REQUIRE(group1.EntryCount() == 1u);
  BOOST_TEST(&group0[0] == entry0_0.get());
  BOOST_TEST(&group0[1] == entry0_1.get());
  BOOST_TEST(&group0[2] == entry0_2.get());
  BOOST_TEST(&group1[0] == entry1_0.get());

  // Test that Update() was called on group0 and group1.
  TestGroupCounts(group0, 1, 1, 1);
  TestGroupCounts(group1, 1, 1, 1);
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

  BOOST_TEST(table.SquaredGroupCount() == 0u);
  table.Update();
  BOOST_TEST_REQUIRE(table.SquaredGroupCount() == 2u);
  TestGroupCounts(table, 1, 1, 2);

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
}

BOOST_AUTO_TEST_CASE(higher_lower_frequency) {
  ImagingTable table;
  test::UniquePtr<ImagingTableEntry> entry1;
  test::UniquePtr<ImagingTableEntry> entry2;
  test::UniquePtr<ImagingTableEntry> entry3;
  entry1->bandStartFrequency = 50.0;
  entry1->bandEndFrequency = 150.0;
  entry2->bandStartFrequency = 150.0;
  entry2->bandEndFrequency = 250.0;
  entry3->bandStartFrequency = 250.0;
  entry3->bandEndFrequency = 350.0;
  BOOST_TEST_REQUIRE(entry1->CentralFrequency() == 100.0);
  BOOST_TEST_REQUIRE(entry2->CentralFrequency() == 200.0);
  BOOST_TEST_REQUIRE(entry3->CentralFrequency() == 300.0);

  // Add the entries out-of-order in this test.
  table.AddEntry(entry1.take());
  table.AddEntry(entry3.take());
  table.AddEntry(entry2.take());

  BOOST_TEST(table.FirstWithHigherFrequency(0.0) == entry1.get());
  BOOST_TEST(table.FirstWithHigherFrequency(199.0) == entry2.get());
  BOOST_TEST(table.FirstWithHigherFrequency(201.0) == entry3.get());
  BOOST_TEST(table.FirstWithHigherFrequency(342.0) == nullptr);

  BOOST_TEST(table.FirstWithLowerFrequency(0.0) == nullptr);
  BOOST_TEST(table.FirstWithLowerFrequency(199.0) == entry1.get());
  BOOST_TEST(table.FirstWithLowerFrequency(201.0) == entry2.get());
  BOOST_TEST(table.FirstWithLowerFrequency(342.0) == entry3.get());
}

BOOST_AUTO_TEST_SUITE_END()
