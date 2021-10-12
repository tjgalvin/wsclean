#include "../system/throwruntimeerror.h"

#include <boost/test/unit_test.hpp>

#include <string>

BOOST_AUTO_TEST_SUITE(throw_runtime_error)

BOOST_AUTO_TEST_CASE(throw_runtime_error) {
  BOOST_CHECK_EXCEPTION(
      ThrowRuntimeError(), std::runtime_error,
      [](const std::runtime_error& e) { return e.what() == std::string(); });

  BOOST_CHECK_EXCEPTION(ThrowRuntimeError('a', "bc", std::string("def"), 123),
                        std::runtime_error, [](const std::runtime_error& e) {
                          return e.what() == std::string("abcdef123");
                        });
}

BOOST_AUTO_TEST_SUITE_END()
