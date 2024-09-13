#include "../../msproviders/reordering.h"
#include "../../msproviders/reorderedfilewriter.h"
#include "aocommon/polarization.h"

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/filesystem.hpp>

#include <stdexcept>
#include <vector>
#include <set>
#include <complex>
#include <cstddef>

using aocommon::Polarization;
using aocommon::PolarizationEnum;
using wsclean::reordering::ChannelRange;
using wsclean::reordering::MetaHeader;
using wsclean::reordering::MetaRecord;
using wsclean::reordering::PartHeader;
using wsclean::reordering::ReorderedFileWriter;
using wsclean::reordering::ReorderedHandleData;

class FixtureDirectory {
 public:
  /// Create the temporary directory and set it as working directory
  FixtureDirectory() {
    boost::filesystem::create_directories(kPath);
    boost::filesystem::current_path(kPath);
  }

  FixtureDirectory(const FixtureDirectory&) = delete;
  FixtureDirectory& operator=(const FixtureDirectory&) = delete;

  /// Remove the temporary diectory
  /// Will always run
  ~FixtureDirectory() {
    boost::filesystem::current_path(kWorkDir);
    boost::filesystem::remove_all(kPath);
  }

 private:
  const boost::filesystem::path kPath = boost::filesystem::unique_path();
  const boost::filesystem::path kWorkDir = boost::filesystem::current_path();
};

using ComplexVector = std::vector<std::complex<float>>;

const ComplexVector kTestData{
    10.0f, 11.0f, 12.0f, 13.0f, 20.0f, 21.0f, 22.0f, 23.0f,
    30.0f, 31.0f, 32.0f, 33.0f, 40.0f, 41.0f, 42.0f, 43.0f,
};
const std::vector<float> kTestWeights{0.1f, 0.1f, 0.1f, 0.1f, 0.2f, 0.2f,
                                      0.2f, 0.2f, 0.3f, 0.3f, 0.3f, 0.3f,
                                      0.4f, 0.4f, 0.4f, 0.4f};
// std::vector<bool> does not have a .data() method
const bool kTestFlags[] = {false, false, false, false, false, false,
                           false, false, true,  true,  true,  true,
                           false, false, false, false};

const std::vector<ChannelRange> kChannelRanges{{0, 0, 2}, {0, 6, 8}, {1, 0, 4}};
const schaapcommon::reordering::MSSelection kSelection;
const std::set<PolarizationEnum> kPolsOut{Polarization::StokesI,
                                          Polarization::StokesQ};
const aocommon::MultiBandData kBands;
const std::string kTemporaryDirectory = "tmp";
const double kStartTime = 1.0;

const ReorderedHandleData kData("test.ms", "DATA", kTemporaryDirectory,
                                kChannelRanges, true, false, kPolsOut,
                                kSelection, kBands, 6, true,
                                [](ReorderedHandleData) {});

const std::map<size_t, std::set<aocommon::PolarizationEnum>>
    kMsPolarizationsPerDataDescId{{0, {Polarization::XX, Polarization::YY}},
                                  {1, {Polarization::XX, Polarization::YY}}};

BOOST_AUTO_TEST_SUITE(reordering_writer)

BOOST_FIXTURE_TEST_CASE(reordering_writer_assert_file_creation,
                        FixtureDirectory) {
  {
    boost::filesystem::create_directory(kTemporaryDirectory);
    ReorderedFileWriter reordering_writer(kData, kMsPolarizationsPerDataDescId,
                                          kStartTime);
  }

  const std::vector<std::string> expected_reorder_files{
      "test.ms-spw0-parted-meta.tmp", "test.ms-spw1-parted-meta.tmp",
      "test.ms-part0000-I-b0.tmp",    "test.ms-part0000-I-b0-w.tmp",
      "test.ms-part0000-I-b0-m.tmp",  "test.ms-part0001-I-b0.tmp",
      "test.ms-part0001-I-b0-w.tmp",  "test.ms-part0001-I-b0-m.tmp",
      "test.ms-part0002-I-b1.tmp",    "test.ms-part0002-I-b1-w.tmp",
      "test.ms-part0002-I-b1-m.tmp",  "test.ms-part0000-Q-b0-m.tmp",
      "test.ms-part0000-Q-b0.tmp",    "test.ms-part0000-Q-b0-w.tmp",
      "test.ms-part0001-Q-b0-m.tmp",  "test.ms-part0001-Q-b0.tmp",
      "test.ms-part0001-Q-b0-w.tmp",  "test.ms-part0002-Q-b1-m.tmp",
      "test.ms-part0002-Q-b1.tmp",    "test.ms-part0002-Q-b1-w.tmp"};

  for (const std::string& reorder_file : expected_reorder_files) {
    std::string file = kTemporaryDirectory + "/" + reorder_file;
    BOOST_CHECK(boost::filesystem::exists(file));
  }
}

BOOST_FIXTURE_TEST_CASE(reordering_writer_assert_file_content,
                        FixtureDirectory) {
  {
    // Write content for spw0 and assert the contents of each file
    boost::filesystem::create_directory(kTemporaryDirectory);
    ReorderedFileWriter reordering_writer(kData, kMsPolarizationsPerDataDescId,
                                          kStartTime);

    // spw 0
    reordering_writer.WriteMetaRow(0.1, 0.2, 0.3, 1.1, 0, 0, 1, 1);

    reordering_writer.WriteDataRow(kTestData.data(), kTestData.data(),
                                   kTestWeights.data(), kTestFlags, 0);

    reordering_writer.UpdateMetaHeaders();
    reordering_writer.UpdatePartHeaders(true);
  }

  // Assert the header of the meta file before updating the headers
  std::ifstream meta_file("tmp/test.ms-spw0-parted-meta.tmp");

  MetaHeader meta_header;
  meta_header.Read(meta_file);
  BOOST_CHECK_CLOSE_FRACTION(meta_header.start_time, 1.0, 1e-5);
  BOOST_CHECK_EQUAL(meta_header.selected_row_count, 1);
  BOOST_CHECK_EQUAL(meta_header.filename_length, 7);

  std::vector<char> ms_path(meta_header.filename_length);
  meta_file.read(ms_path.data(), meta_header.filename_length);
  std::string ms_path_str(ms_path.begin(), ms_path.end());

  BOOST_CHECK_EQUAL(ms_path_str, "test.ms");

  // Assert meta row
  MetaRecord meta_record;
  meta_record.Read(meta_file);

  BOOST_CHECK_CLOSE_FRACTION(meta_record.u, 0.1, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(meta_record.v, 0.2, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(meta_record.w, 0.3, 1e-5);
  BOOST_CHECK_CLOSE_FRACTION(meta_record.time, 1.1, 1e-5);
  BOOST_CHECK_EQUAL(meta_record.antenna1, 0);
  BOOST_CHECK_EQUAL(meta_record.antenna2, 1);
  BOOST_CHECK_EQUAL(meta_record.field_id, 1);

  // Assert the contents of data
  std::ifstream data_file("tmp/test.ms-part0000-I-b0.tmp");
  PartHeader part_header;
  part_header.Read(data_file);
  BOOST_CHECK_EQUAL(part_header.channel_count, 2);
  BOOST_CHECK_EQUAL(part_header.channel_start, 0);
  BOOST_CHECK_EQUAL(part_header.data_desc_id, 0);
  BOOST_CHECK(part_header.has_model);

  constexpr size_t kMaxChannels = 2;
  std::vector<std::complex<float>> data_buffer(kMaxChannels);
  std::vector<std::complex<float>> expected_data_buffer{{10.5f, 0.0f},
                                                        {12.5f, 0.0f}};

  data_file.read(reinterpret_cast<char*>(data_buffer.data()),
                 2 * sizeof(std::complex<float>));
  BOOST_CHECK_EQUAL_COLLECTIONS(data_buffer.begin(), data_buffer.end(),
                                expected_data_buffer.begin(),
                                expected_data_buffer.end());

  std::ifstream model_file("tmp/test.ms-part0000-I-b0-m.tmp");
  std::vector<std::complex<float>> model_buffer(kMaxChannels);
  std::vector<std::complex<float>> expected_model_buffer{{10.5f, 0.0f},
                                                         {12.5f, 0.0f}};

  model_file.read(reinterpret_cast<char*>(model_buffer.data()),
                  2 * sizeof(std::complex<float>));
  BOOST_CHECK_EQUAL_COLLECTIONS(model_buffer.begin(), model_buffer.end(),
                                expected_model_buffer.begin(),
                                expected_model_buffer.end());

  std::ifstream weight_file("tmp/test.ms-part0000-I-b0-w.tmp");
  std::vector<float> weight_buffer(kMaxChannels);
  std::vector<float> expected_weight_buffer{0.4f, 0.4f};

  weight_file.read(reinterpret_cast<char*>(weight_buffer.data()),
                   2 * sizeof(std::complex<float>));
  BOOST_CHECK_EQUAL_COLLECTIONS(weight_buffer.begin(), weight_buffer.end(),
                                expected_weight_buffer.begin(),
                                expected_weight_buffer.end());
}

BOOST_FIXTURE_TEST_CASE(reordering_writer_assert_zeros_model_creation,
                        FixtureDirectory) {
  {
    boost::filesystem::create_directory(kTemporaryDirectory);
    const ReorderedHandleData handle_no_model_column(
        "test.ms", "DATA", kTemporaryDirectory, kChannelRanges, false, false,
        kPolsOut, kSelection, kBands, 6, true, [](ReorderedHandleData) {});

    ReorderedFileWriter reordering_writer(
        handle_no_model_column, kMsPolarizationsPerDataDescId, kStartTime);

    reordering_writer.WriteMetaRow(0.1, 0.2, 0.3, 1.1, 0, 0, 1, 1);

    reordering_writer.WriteDataRow(kTestData.data(), kTestData.data(),
                                   kTestWeights.data(), kTestFlags, 0);

    reordering_writer.UpdateMetaHeaders();
    reordering_writer.UpdatePartHeaders(true);
    reordering_writer.PopulateModelWithZeros([](size_t, size_t) {});
  }

  constexpr size_t kMaxChannels = 2;
  std::ifstream model_file("tmp/test.ms-part0000-I-b0-m.tmp");
  std::vector<std::complex<float>> model_buffer(kMaxChannels);
  std::vector<std::complex<float>> expected_model_buffer(kMaxChannels,
                                                         {0.0f, 0.0f});

  model_file.read(reinterpret_cast<char*>(model_buffer.data()),
                  2 * sizeof(std::complex<float>));
  BOOST_CHECK_EQUAL_COLLECTIONS(model_buffer.begin(), model_buffer.end(),
                                expected_model_buffer.begin(),
                                expected_model_buffer.end());
}

BOOST_AUTO_TEST_SUITE_END()
