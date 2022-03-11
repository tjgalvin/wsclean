#include "../../system/mappedfile.h"

#include <boost/filesystem/operations.hpp>
#include <boost/test/unit_test.hpp>

#include <fstream>

BOOST_AUTO_TEST_SUITE(mapped_file)

BOOST_AUTO_TEST_CASE(map) {
  constexpr size_t kSize = 16;
  constexpr const char* kFilename = "tmappedfile-test.tmp";

  std::ofstream file(kFilename);
  char data[kSize];
  std::fill_n(data, kSize, 1);
  file.write(data, kSize);
  file.close();

  MappedFile file_a;
  MappedFile file_b(kFilename, kSize);
  for (size_t i = 0; i != kSize; ++i) {
    file_b.Data()[i] = i + 7;
  }
  file_a = std::move(file_b);
  for (size_t i = 0; i != kSize; ++i) {
    BOOST_CHECK_EQUAL(file_a.Data()[i], i + 7);
  }

  boost::filesystem::remove(kFilename);
}

BOOST_AUTO_TEST_SUITE_END()
