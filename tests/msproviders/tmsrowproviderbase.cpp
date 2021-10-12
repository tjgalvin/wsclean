#include "../../msproviders/msrowproviderbase.h"

#include "../../msproviders/directmsrowprovider.h"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(ms_row_provider_base)

BOOST_AUTO_TEST_CASE(create_ms_row_proverider_invalid_ms) {
  BOOST_CHECK_THROW(
      ((void)MakeMsRowProvider("invalid_measurement_set", MSSelection{},
                               std::map<size_t, size_t>{}, "DATA", false)),
      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(create_ms_row_proverider_direct_ms) {
  std::unique_ptr<MsRowProviderBase> provider =
      MakeMsRowProvider("test_data/MWA_MOCK.ms", MSSelection{},
                        std::map<size_t, size_t>{}, "DATA", false);

  BOOST_REQUIRE(provider);
  BOOST_CHECK_NO_THROW(dynamic_cast<DirectMSRowProvider&>(*provider));
}

// TODO AST-630 Add test for BDA row provider

BOOST_AUTO_TEST_SUITE_END()
