#include "../../msproviders/bdamsrowprovider.h"

#include "tbdamsrowproviderdata.h"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(bda_ms_row_provider)

BOOST_AUTO_TEST_CASE(bda_ms_row_provider_constructor_no_bda_tables) {
  BOOST_CHECK_EXCEPTION(
      (BdaMsRowProvider{{"test_data/MWA_MOCK.ms"},
                        MSSelection{},
                        std::map<size_t, size_t>{{0, 0}},
                        "DATA",
                        false}),
      std::runtime_error, [](const std::runtime_error& e) {
        return e.what() ==
               std::string(
                   "A BDA measurement set requires a BDA_FACTORS table.");
      });
}

static void create_bda_ms_row_provider_with_selection(
    const MSSelection& selection) {
  BdaMsRowProvider{{"test_data/MWA_BDA_MOCK.ms"},
                   selection,
                   std::map<size_t, size_t>{{0, 0}},
                   "DATA",
                   false};
}

static void create_bda_ms_row_provider_with_selection_interval() {
  MSSelection selection;
  selection.SetInterval(0, 1);
  create_bda_ms_row_provider_with_selection(selection);
}

static void create_bda_ms_row_provider_with_selection_even_timesteps() {
  MSSelection selection;
  selection.SetEvenOrOddTimesteps(MSSelection::EvenTimesteps);
  create_bda_ms_row_provider_with_selection(selection);
}

static void create_bda_ms_row_provider_with_selection_odd_timesteps() {
  MSSelection selection;
  selection.SetEvenOrOddTimesteps(MSSelection::OddTimesteps);
  create_bda_ms_row_provider_with_selection(selection);
}

BOOST_AUTO_TEST_CASE(bda_ms_row_provider_constructor_invalid_selection) {
  BOOST_CHECK_EXCEPTION(create_bda_ms_row_provider_with_selection_interval(),
                        std::runtime_error, [](const std::runtime_error& e) {
                          return e.what() ==
                                 std::string(
                                     "An interval selection isn't supported "
                                     "for a BDA measurement set.");
                        });
  BOOST_CHECK_EXCEPTION(
      create_bda_ms_row_provider_with_selection_even_timesteps(),
      std::runtime_error, [](const std::runtime_error& e) {
        return e.what() == std::string(
                               "An interval selection isn't supported "
                               "for a BDA measurement set.");
      });
  BOOST_CHECK_EXCEPTION(
      create_bda_ms_row_provider_with_selection_odd_timesteps(),
      std::runtime_error, [](const std::runtime_error& e) {
        return e.what() == std::string(
                               "An interval selection isn't supported "
                               "for a BDA measurement set.");
      });
}

BOOST_AUTO_TEST_CASE(bda_ms_row_provider) {
  BdaMsRowProvider provider{{"test_data/MWA_BDA_MOCK.ms"},
                            MSSelection{},
                            std::map<size_t, size_t>{{0, 0}},
                            "DATA",
                            false};

  BOOST_REQUIRE_EQUAL(provider.BeginRow(), 0);
  BOOST_REQUIRE_EQUAL(provider.EndRow(), 21);
  BOOST_CHECK_EQUAL(provider.CurrentProgress(), 0);

  for (const auto& row : MwaBdaMockMs::kMs) {
    BdaMsRowProvider::DataArray data;
    BdaMsRowProvider::FlagArray flag;
    BdaMsRowProvider::WeightArray weight;
    std::array<double, 3> uvw;
    uint32_t data_description_id;
    std::array<uint32_t, 2> antenna;
    uint32_t field_id;
    double time;
    provider.ReadData(data, flag, weight, uvw[0], uvw[1], uvw[2],
                      data_description_id, antenna[0], antenna[1], field_id,
                      time);

    BOOST_CHECK_EQUAL(uint64_t(time), row[MwaBdaMockMs::kTime]);
    BOOST_CHECK_EQUAL(antenna[0], row[MwaBdaMockMs::kAntenna1]);
    BOOST_CHECK_EQUAL(antenna[1], row[MwaBdaMockMs::kAntenna2]);

    // TODO Add tests for the model.

    // The last row in the real MWA_BDA_MOCK.ms is an autocorrelation so
    // always go to the next row to make sure we skip it.
    provider.NextRow();
  }
  BOOST_CHECK(provider.AtEnd());
  BOOST_CHECK_EQUAL(provider.CurrentProgress(), 21);
}

BOOST_AUTO_TEST_SUITE_END()
