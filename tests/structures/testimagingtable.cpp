#include <boost/test/unit_test.hpp>

#include "../../structures/imagingtable.h"

BOOST_AUTO_TEST_SUITE(imagingtable)

BOOST_AUTO_TEST_CASE(constructor) {
  ImagingTable table;
  BOOST_TEST(table.IndependentGroupCount() == 0u);
  BOOST_TEST(table.SquaredGroupCount() == 0u);
  BOOST_TEST(table.EntryCount() == 0u);
}

BOOST_AUTO_TEST_CASE(add_update_clear) {
  ImagingTable table;
  ImagingTableEntry& entry1 = table.AddEntry();
  ImagingTableEntry& entry2 = table.AddEntry();
  entry1.joinedGroupIndex = 42;
  entry2.joinedGroupIndex = 43;
  entry1.squaredDeconvolutionIndex = 142;
  entry2.squaredDeconvolutionIndex = 142;

  BOOST_TEST_REQUIRE(table.EntryCount() == 2u);
  BOOST_TEST(&table[0] == &entry1);
  BOOST_TEST(&table[1] == &entry2);
  BOOST_TEST(&table.Front() == &entry1);

  // No Update called -> only EntryCount is correct.
  BOOST_TEST(table.IndependentGroupCount() == 0u);
  BOOST_TEST(table.SquaredGroupCount() == 0u);

  table.Update();
  BOOST_TEST(table.IndependentGroupCount() == 2u);
  BOOST_TEST(table.SquaredGroupCount() == 1u);

  table.Clear();
  BOOST_TEST(table.EntryCount() == 0u);
  BOOST_TEST(table.IndependentGroupCount() == 0u);
  BOOST_TEST(table.SquaredGroupCount() == 0u);
}

BOOST_AUTO_TEST_CASE(independent_groups) {
  ImagingTable table;
  // Deliberately mix the groups when adding the entries:
  ImagingTableEntry& entry0_0 = table.AddEntry();
  ImagingTableEntry& entry1_0 = table.AddEntry();
  ImagingTableEntry& entry1_1 = table.AddEntry();
  ImagingTableEntry& entry0_1 = table.AddEntry();
  entry0_0.joinedGroupIndex = 42;
  entry0_1.joinedGroupIndex = 42;
  entry1_0.joinedGroupIndex = 43;
  entry1_1.joinedGroupIndex = 43;
  BOOST_TEST(table.IndependentGroupCount() == 0u);
  table.Update();
  BOOST_TEST_REQUIRE(table.IndependentGroupCount() == 2u);
  BOOST_TEST(table.SquaredGroupCount() == 1u);

  ImagingTable group0 = table.GetIndependentGroup(0);
  ImagingTable group1 = table.GetIndependentGroup(1);
  BOOST_TEST_REQUIRE(group0.EntryCount() == 2u);
  BOOST_TEST_REQUIRE(group1.EntryCount() == 2u);
  BOOST_TEST(&group0[0] == &entry0_0);
  BOOST_TEST(&group0[1] == &entry0_1);
  BOOST_TEST(&group1[0] == &entry1_0);
  BOOST_TEST(&group1[1] == &entry1_1);

  // Test that Update() was called on group0 and group1.
  BOOST_TEST(group0.IndependentGroupCount() == 1u);
  BOOST_TEST(group1.IndependentGroupCount() == 1u);
  BOOST_TEST(group0.SquaredGroupCount() == 1u);
  BOOST_TEST(group1.SquaredGroupCount() == 1u);
}

BOOST_AUTO_TEST_CASE(squared_groups) {
  ImagingTable table;
  // Deliberately mix the groups when adding the entries:
  ImagingTableEntry& entry0_0 = table.AddEntry();
  ImagingTableEntry& entry1_0 = table.AddEntry();
  ImagingTableEntry& entry0_1 = table.AddEntry();
  ImagingTableEntry& entry0_2 = table.AddEntry();
  entry0_0.squaredDeconvolutionIndex = 42;
  entry0_1.squaredDeconvolutionIndex = 42;
  entry0_2.squaredDeconvolutionIndex = 42;
  entry1_0.squaredDeconvolutionIndex = 43;
  BOOST_TEST(table.SquaredGroupCount() == 0u);
  table.Update();
  BOOST_TEST_REQUIRE(table.SquaredGroupCount() == 2u);
  BOOST_TEST(table.IndependentGroupCount() == 1u);

  ImagingTable group0 = table.GetSquaredGroup(0);
  ImagingTable group1 = table.GetSquaredGroup(1);
  BOOST_TEST_REQUIRE(group0.EntryCount() == 3u);
  BOOST_TEST_REQUIRE(group1.EntryCount() == 1u);
  BOOST_TEST(&group0[0] == &entry0_0);
  BOOST_TEST(&group0[1] == &entry0_1);
  BOOST_TEST(&group0[2] == &entry0_2);
  BOOST_TEST(&group1[0] == &entry1_0);

  // Test that Update() was called on group0 and group1.
  BOOST_TEST(group0.IndependentGroupCount() == 1u);
  BOOST_TEST(group1.IndependentGroupCount() == 1u);
  BOOST_TEST(group0.SquaredGroupCount() == 1u);
  BOOST_TEST(group1.SquaredGroupCount() == 1u);
}

BOOST_AUTO_TEST_CASE(higher_lower_frequency) {
  ImagingTable table;
  // Add the entries out-of-order in this test.
  ImagingTableEntry& entry1 = table.AddEntry();
  ImagingTableEntry& entry3 = table.AddEntry();
  ImagingTableEntry& entry2 = table.AddEntry();

  entry1.bandStartFrequency = 50.0;
  entry1.bandEndFrequency = 150.0;
  entry2.bandStartFrequency = 150.0;
  entry2.bandEndFrequency = 250.0;
  entry3.bandStartFrequency = 250.0;
  entry3.bandEndFrequency = 350.0;
  BOOST_TEST_REQUIRE(entry1.CentralFrequency() == 100.0);
  BOOST_TEST_REQUIRE(entry2.CentralFrequency() == 200.0);
  BOOST_TEST_REQUIRE(entry3.CentralFrequency() == 300.0);

  BOOST_TEST(table.FirstWithHigherFrequency(0.0) == &entry1);
  BOOST_TEST(table.FirstWithHigherFrequency(199.0) == &entry2);
  BOOST_TEST(table.FirstWithHigherFrequency(201.0) == &entry3);
  BOOST_TEST(table.FirstWithHigherFrequency(342.0) == nullptr);

  BOOST_TEST(table.FirstWithLowerFrequency(0.0) == nullptr);
  BOOST_TEST(table.FirstWithLowerFrequency(199.0) == &entry1);
  BOOST_TEST(table.FirstWithLowerFrequency(201.0) == &entry2);
  BOOST_TEST(table.FirstWithLowerFrequency(342.0) == &entry3);
}

BOOST_AUTO_TEST_SUITE_END()
