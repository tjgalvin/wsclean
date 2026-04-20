#include "../../structures/preallocatedvector.h"

#include <boost/test/unit_test.hpp>

namespace wsclean {

BOOST_AUTO_TEST_SUITE(preallocated_vector)

BOOST_AUTO_TEST_CASE(default_construction) {
  PreallocatedVector<int> test_vector;
  BOOST_CHECK_EQUAL(test_vector.Size(), 0);
  BOOST_CHECK_EQUAL(test_vector.Capacity(), 0);
  BOOST_CHECK(test_vector.Empty());
}

BOOST_AUTO_TEST_CASE(reserve_and_push_back) {
  PreallocatedVector<int> test_vector;
  test_vector.Reserve(10);

  BOOST_CHECK_EQUAL(test_vector.Capacity(), 10);
  BOOST_CHECK_EQUAL(test_vector.Size(), 0);

  test_vector.PushBack(42);
  BOOST_CHECK_EQUAL(test_vector.Size(), 1);
  BOOST_CHECK_EQUAL(test_vector[0], 42);

  test_vector.PushBack(100);
  BOOST_CHECK_EQUAL(test_vector.Size(), 2);
  BOOST_CHECK_EQUAL(test_vector[1], 100);
}

BOOST_AUTO_TEST_CASE(element_access) {
  PreallocatedVector<std::string> test_vector;
  test_vector.Reserve(3);

  test_vector.PushBack("first");
  test_vector.PushBack("second");
  test_vector.PushBack("third");

  BOOST_CHECK_EQUAL(test_vector[0], "first");
  BOOST_CHECK_EQUAL(test_vector[1], "second");
  BOOST_CHECK_EQUAL(test_vector[2], "third");
}

BOOST_AUTO_TEST_CASE(const_access) {
  PreallocatedVector<int> test_vector;
  test_vector.Reserve(2);
  test_vector.PushBack(10);
  test_vector.PushBack(20);

  const PreallocatedVector<int>& const_vector = test_vector;
  BOOST_CHECK_EQUAL(const_vector[0], 10);
  BOOST_CHECK_EQUAL(const_vector[1], 20);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace wsclean
