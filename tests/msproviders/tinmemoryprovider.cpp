#include "../../msproviders/inmemoryprovider.h"
#include "../../msproviders/msreaders/msreader.h"
#include "../../structures/inmemorypart.h"

#include <boost/test/unit_test.hpp>

namespace wsclean {
namespace {
const std::string kDescription = "In memory data part";
const std::vector<std::string> kAntennas{"Antenna1", "Antenna2", "Antenna3"};
const double kInterval = 1.0;
const double kRa = 2.0;
const double kDec = -1.0;
const std::string kTelescopeName = "Large super telescope";
const double kStartTime = 1000.0;

ObservationInfo GetObservationInfo() {
  ObservationInfo observation_info;
  observation_info.fieldName = "Field";
  observation_info.observer = "Observer";
  observation_info.phaseCentreDec = kDec;
  observation_info.phaseCentreRA = kRa;
  observation_info.telescopeName = kTelescopeName;
  return observation_info;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(in_memory_provider)

BOOST_AUTO_TEST_CASE(read_write) {
  InMemoryPart data;
  data.description = kDescription;
  data.info = std::make_shared<InMemoryInfo>();
  data.info->antenna_names = kAntennas;
  data.info->has_frequency_bda = true;
  data.info->interval = kInterval;
  data.info->observation_info = GetObservationInfo();
  data.info->start_time = kStartTime;
  data.is_regular = false;
  data.polarization = aocommon::PolarizationEnum::DiagonalInstrumental;

  const std::vector<aocommon::ChannelInfo> channels3{
      {100e6, 1e6}, {110e6, 1e6}, {120e6, 1e6}};
  const aocommon::BandData band3(channels3, channels3[0].Frequency());
  data.selected_bands.SetBand(3, band3);

  const std::vector<aocommon::ChannelInfo> channels4{{105e6, 1e6},
                                                     {115e6, 1e6}};
  const aocommon::BandData band4(channels4, channels4[0].Frequency());
  data.selected_bands.SetBand(4, band4);

  data.meta_data = std::make_shared<std::vector<MSProvider::MetaData>>();

  MSProvider::MetaData meta0;
  meta0.antenna1 = 0;
  meta0.antenna2 = 1;
  meta0.data_desc_id = 3;
  meta0.field_id = 4;
  meta0.time = kStartTime + 1.0;
  meta0.u_in_m = -300.0;
  meta0.v_in_m = 400.0;
  meta0.w_in_m = -500.0;
  data.meta_data->emplace_back(meta0);

  MSProvider::MetaData meta2 = meta0;
  meta2.antenna2 = 2;
  data.meta_data->emplace_back(meta2);

  MSProvider::MetaData meta3 = meta2;
  meta3.antenna1 = 1;
  meta3.data_desc_id = 4;
  data.meta_data->emplace_back(meta3);
  data.meta_data->emplace_back(meta0);

  std::vector<std::complex<float>> data0{42.0, 43.0, 44.0, 45.0, 46.0, 47.0};
  InMemoryPartRow row0;
  row0.data = data0;
  row0.model_data.assign(data0.size(), 0.0);
  row0.weights.assign(data0.size(), 13.0f);
  data.rows.emplace_back(row0);

  std::vector<std::complex<float>> data1{48.0, 49.0, 50.0, 51.0, 52.0, 53.0};
  InMemoryPartRow row1;
  row1.data = data1;
  row1.model_data.assign(data1.size(), 0.0);
  row1.weights.assign(data1.size(), 14.0f);
  data.rows.emplace_back(row1);

  std::vector<std::complex<float>> data2{54.0, 55.0, 56.0, 57.0};
  InMemoryPartRow row2;
  row2.data = data2;
  row2.model_data.assign(data2.size(), 0.0);
  row2.weights.assign(data2.size(), 15.0f);
  data.rows.emplace_back(row2);

  std::vector<std::complex<float>> data3{58.0, 59.0, 60.0, 61.0, 62.0, 63.0};
  InMemoryPartRow row3;
  row3.data = data3;
  row3.model_data.assign(data3.size(), 0.0);
  row3.weights.assign(data3.size(), 16.0f);
  data.rows.emplace_back(row3);

  InMemoryProvider provider(data);
  BOOST_CHECK(provider.GetAntennaNames() == kAntennas);
  BOOST_CHECK(provider.GetObservationInfo() == GetObservationInfo());
  BOOST_CHECK_EQUAL(provider.Interval(), kInterval);
  BOOST_CHECK_EQUAL(provider.NAntennas(), kAntennas.size());
  BOOST_CHECK_EQUAL(provider.NMaxChannels(), 3);
  BOOST_CHECK_EQUAL(provider.NPolarizations(), 2);
  BOOST_CHECK_EQUAL(provider.NRows(), 4);
  BOOST_CHECK_EQUAL(provider.Polarization(),
                    aocommon::PolarizationEnum::DiagonalInstrumental);
  BOOST_CHECK_EQUAL(provider.SelectedBands().BandCount(), 2);
  BOOST_CHECK_EQUAL(provider.SelectedBands()[3].ChannelCount(), 3);
  BOOST_CHECK_EQUAL(provider.SelectedBands()[4].ChannelCount(), 2);
  BOOST_CHECK_EQUAL(provider.StartTime(), kStartTime);
  BOOST_CHECK(provider.HasFrequencyBda());
  BOOST_CHECK(!provider.IsRegular());
  BOOST_CHECK_EQUAL(provider.PartDescription(), kDescription);
  BOOST_CHECK(!provider.MsIfAvailable());

  std::vector<std::complex<float>> ones(6, 1.0f);
  // The last row is on purpose not written to also test an early reset.
  for (size_t i = 0; i != 3; ++i) {
    provider.WriteModel(ones.data(), false);
    provider.NextOutputRow();
  }

  provider.ResetWritePosition();

  for (size_t i = 0; i != provider.NRows(); ++i) {
    provider.WriteModel(ones.data(), true);
    provider.NextOutputRow();
  }

  std::unique_ptr<wsclean::MSReader> reader = provider.MakeReader();

  std::complex<float> buffer[6];
  float weights[6];
  MSProvider::MetaData metadata;

  BOOST_REQUIRE(reader->CurrentRowAvailable());
  BOOST_CHECK_EQUAL(reader->RowId(), 0);
  reader->ReadData(buffer);
  BOOST_CHECK_EQUAL(buffer[5], 47.0f);
  reader->ReadModel(buffer);
  BOOST_CHECK_EQUAL(buffer[5], 2.0f);
  reader->ReadWeights(weights);
  BOOST_CHECK_EQUAL(weights[5], 13.0f);
  reader->ReadMeta(metadata);
  BOOST_CHECK_EQUAL(metadata.antenna1, 0);
  BOOST_CHECK_EQUAL(metadata.antenna2, 1);
  reader->NextInputRow();

  BOOST_REQUIRE(reader->CurrentRowAvailable());
  BOOST_CHECK_EQUAL(reader->RowId(), 1);
  reader->ReadData(buffer);
  BOOST_CHECK_EQUAL(buffer[5], 53.0f);
  reader->ReadModel(buffer);
  BOOST_CHECK_EQUAL(buffer[5], 2.0f);
  reader->ReadWeights(weights);
  BOOST_CHECK_EQUAL(weights[5], 14.0f);
  reader->ReadMeta(metadata);
  BOOST_CHECK_EQUAL(metadata.antenna1, 0);
  BOOST_CHECK_EQUAL(metadata.antenna2, 2);
  reader->NextInputRow();

  BOOST_REQUIRE(reader->CurrentRowAvailable());
  BOOST_CHECK_EQUAL(reader->RowId(), 2);
  reader->ReadData(buffer);
  BOOST_CHECK_EQUAL(buffer[3], 57.0f);
  // This row only has 4 vis, check if it didn't overwrite the old value:
  BOOST_CHECK_EQUAL(buffer[5], 2.0f);
  buffer[5] = -1.0f;
  reader->ReadModel(buffer);
  BOOST_CHECK_EQUAL(buffer[3], 2.0f);
  BOOST_CHECK_EQUAL(buffer[5], -1.0f);
  reader->ReadWeights(weights);
  BOOST_CHECK_EQUAL(weights[3], 15.0f);
  BOOST_CHECK_EQUAL(weights[5], 14.0f);
  reader->ReadMeta(metadata);
  BOOST_CHECK_EQUAL(metadata.antenna1, 1);
  BOOST_CHECK_EQUAL(metadata.antenna2, 2);
  reader->NextInputRow();

  BOOST_REQUIRE(reader->CurrentRowAvailable());
  BOOST_CHECK_EQUAL(reader->RowId(), 3);
  reader->ReadData(buffer);
  BOOST_CHECK_EQUAL(buffer[5], 63.0f);
  reader->ReadModel(buffer);
  BOOST_CHECK_EQUAL(buffer[5], 1.0f);
  reader->ReadWeights(weights);
  BOOST_CHECK_EQUAL(weights[5], 16.0f);
  reader->ReadMeta(metadata);
  BOOST_CHECK_EQUAL(metadata.antenna1, 0);
  BOOST_CHECK_EQUAL(metadata.antenna2, 1);
  reader->NextInputRow();

  BOOST_CHECK(!reader->CurrentRowAvailable());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace wsclean
