#include "../../msproviders/rowproviders/bdamsrowprovider.h"

#include "tbdamsrowproviderdata.h"

#include <boost/test/unit_test.hpp>

using schaapcommon::reordering::MSSelection;

namespace wsclean {

BOOST_AUTO_TEST_SUITE(bda_ms_row_provider)

BOOST_AUTO_TEST_CASE(bda_ms_row_provider_constructor_no_bda_tables) {
  BOOST_CHECK_THROW((BdaMsRowProvider({"test_data/MWA_MOCK.ms"}, MSSelection{},
                                      std::map<size_t, size_t>{{0, 0}}, "DATA",
                                      "MODEL_DATA", false)),
                    std::runtime_error);
}

static void CreateBdaMsRowProviderWithSelection(const MSSelection& selection) {
  BdaMsRowProvider({"test_data/MWA_BDA_MOCK.ms"}, selection,
                   std::map<size_t, size_t>{{0, 0}}, "DATA", "MODEL_DATA",
                   false);
}

static void CreateBdaMsRowProviderWithSelectionInterval() {
  MSSelection selection;
  selection.SetInterval(0, 1);
  CreateBdaMsRowProviderWithSelection(selection);
}

static void CreateBdaMsRowProviderWithSelectionEvenTimesteps() {
  MSSelection selection;
  selection.SetEvenOrOddTimesteps(MSSelection::kEvenTimesteps);
  CreateBdaMsRowProviderWithSelection(selection);
}

static void CreateBdaMsRowProviderWithSelectionOddTimesteps() {
  MSSelection selection;
  selection.SetEvenOrOddTimesteps(MSSelection::kOddTimesteps);
  CreateBdaMsRowProviderWithSelection(selection);
}

BOOST_AUTO_TEST_CASE(bda_ms_row_provider_constructor_invalid_selection) {
  BOOST_CHECK_EXCEPTION(CreateBdaMsRowProviderWithSelectionInterval(),
                        std::runtime_error, [](const std::runtime_error& e) {
                          return e.what() ==
                                 std::string(
                                     "An interval selection isn't supported "
                                     "for a BDA measurement set.");
                        });
  BOOST_CHECK_EXCEPTION(CreateBdaMsRowProviderWithSelectionEvenTimesteps(),
                        std::runtime_error, [](const std::runtime_error& e) {
                          return e.what() ==
                                 std::string(
                                     "An interval selection isn't supported "
                                     "for a BDA measurement set.");
                        });
  BOOST_CHECK_EXCEPTION(CreateBdaMsRowProviderWithSelectionOddTimesteps(),
                        std::runtime_error, [](const std::runtime_error& e) {
                          return e.what() ==
                                 std::string(
                                     "An interval selection isn't supported "
                                     "for a BDA measurement set.");
                        });
}

struct RowData {
  bool is_found;
  size_t time;  // our test set has integer values
  size_t antenna1;
  size_t antenna2;
};

BOOST_AUTO_TEST_CASE(bda_ms_row_provider) {
  BdaMsRowProvider provider({"test_data/MWA_BDA_MOCK.ms"}, MSSelection{},
                            std::map<size_t, size_t>{{0, 0}}, "DATA",
                            "MODEL_DATA", false);

  BOOST_REQUIRE_EQUAL(provider.BeginRow(), 0);
  BOOST_REQUIRE_EQUAL(provider.EndRow(), 21);
  BOOST_CHECK_EQUAL(provider.CurrentProgress(), 0);

  std::vector<RowData> rows;
  for (const auto& row : MwaBdaMockMs::kMs)
    rows.emplace_back(RowData{false, row[0], row[1], row[2]});

  double previous_time = 0;
  for (size_t i = 0; i != MwaBdaMockMs::kMs.size(); ++i) {
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
    BOOST_CHECK_LE(previous_time, time);
    previous_time = time;

    // We don't know exactly in which order the rows come out: they should be
    // time sorted, but there are multiple rows with the same time, and they may
    // be returned in any order. Therefore, this checks just if all rows are
    // read and no rows are duplicated.
    for (RowData& data : rows) {
      if (time == data.time && antenna[0] == data.antenna1 &&
          antenna[1] == data.antenna2) {
        // Each row should be provided only once
        BOOST_CHECK(!data.is_found);
        data.is_found = true;
      }
    }

    // TODO Add tests for the model.

    // The last row in the real MWA_BDA_MOCK.ms is an autocorrelation so
    // always go to the next row to make sure we skip it.
    provider.NextRow();
  }
  BOOST_CHECK(provider.AtEnd());
  BOOST_CHECK_EQUAL(provider.CurrentProgress(), 21);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace wsclean
