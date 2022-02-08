#include <boost/test/unit_test.hpp>

#include "../deconvolution/componentlist.h"

BOOST_AUTO_TEST_SUITE(component_list)

struct ComponentListFixture {
  ComponentListFixture() : list(512, 512, 4, 3) {
    aocommon::UVector<float> values = {1.0, 2.0, 3.0};
    list.Add(256, 256, 1, values.data());
    values = {5.0, 6.0, 7.0};
    list.Add(256, 256, 1, values.data());
    values = {8.0, 9.0, 10.0};
    list.Add(511, 511, 0, values.data());
    values = {11.0, 12.0, 13.0};
    list.Add(13, 42, 3, values.data());

    list.MergeDuplicates();
  }

  ComponentList list;
  aocommon::UVector<float> values;
};

BOOST_FIXTURE_TEST_CASE(adding_values, ComponentListFixture) {
  BOOST_CHECK_EQUAL(list.ComponentCount(0), 1u);
  BOOST_CHECK_EQUAL(list.ComponentCount(1), 1u);
  BOOST_CHECK_EQUAL(list.ComponentCount(2), 0u);
  BOOST_CHECK_EQUAL(list.ComponentCount(3), 1u);

  size_t x;
  size_t y;
  aocommon::UVector<float> values(3);

  list.GetComponent(0, 0, x, y, values.data());
  BOOST_CHECK_EQUAL(x, 511u);
  BOOST_CHECK_EQUAL(y, 511u);
  BOOST_CHECK_CLOSE_FRACTION(values[0], 8.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(values[1], 9.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(values[2], 10.0, 1e-5);

  list.GetComponent(1, 0, x, y, values.data());
  BOOST_CHECK_EQUAL(x, 256u);
  BOOST_CHECK_EQUAL(y, 256u);
  BOOST_CHECK_CLOSE_FRACTION(values[0], 6.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(values[1], 8.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(values[2], 10.0, 1e-5);

  list.GetComponent(3, 0, x, y, values.data());
  BOOST_CHECK_EQUAL(x, 13u);
  BOOST_CHECK_EQUAL(y, 42u);
  BOOST_CHECK_CLOSE_FRACTION(values[0], 11.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(values[1], 12.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(values[2], 13.0, 1e-5);
}

BOOST_FIXTURE_TEST_CASE(get_position, ComponentListFixture) {
  BOOST_CHECK_EQUAL(list.GetPositions(0).size(), 1u);
  BOOST_CHECK_EQUAL(list.GetPositions(1).size(), 1u);
  BOOST_CHECK_EQUAL(list.GetPositions(2).size(), 0u);
  BOOST_CHECK_EQUAL(list.GetPositions(3).size(), 1u);

  BOOST_CHECK_EQUAL(list.GetPositions(0)[0].x, 511u);
  BOOST_CHECK_EQUAL(list.GetPositions(0)[0].y, 511u);

  BOOST_CHECK_EQUAL(list.GetPositions(1)[0].x, 256u);
  BOOST_CHECK_EQUAL(list.GetPositions(1)[0].y, 256u);

  BOOST_CHECK_EQUAL(list.GetPositions(3)[0].x, 13u);
  BOOST_CHECK_EQUAL(list.GetPositions(3)[0].y, 42u);
}

BOOST_FIXTURE_TEST_CASE(multiply_scale_component, ComponentListFixture) {
  for (size_t i = 0; i != list.NScales(); ++i) {
    // Second scale is empty (captured in DEBUG mode via an assert)
    if (i != 2u) {
      for (size_t j = 0; j != list.NFrequencies(); ++j) {
        list.MultiplyScaleComponent(i, 0, j, static_cast<double>(j + 1u));
      }
    }
  }

  size_t x;
  size_t y;
  aocommon::UVector<float> values(3);

  list.GetComponent(0, 0, x, y, values.data());
  BOOST_CHECK_CLOSE_FRACTION(values[0], 8.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(values[1], 18.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(values[2], 30.0, 1e-5);

  list.GetComponent(1, 0, x, y, values.data());
  BOOST_CHECK_CLOSE_FRACTION(values[0], 6.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(values[1], 16.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(values[2], 30.0, 1e-5);

  list.GetComponent(3, 0, x, y, values.data());
  BOOST_CHECK_CLOSE_FRACTION(values[0], 11.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(values[1], 24.0, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(values[2], 39.0, 1e-5);
}

BOOST_AUTO_TEST_SUITE_END()
