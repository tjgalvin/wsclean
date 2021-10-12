#include "msrowproviderbase.h"

#include "directmsrowprovider.h"
#include "msprovider.h"
#include "../system/throwruntimeerror.h"

#include <boost/make_unique.hpp>

MsRowProviderBase::MsRowProviderBase(const casacore::MeasurementSet& ms,
                                     const MSSelection& selection,
                                     const std::string& data_column_name)
    : ms_(ms), selection_(selection), columns_(ms_, data_column_name) {
  MSProvider::GetRowRange(ms_, selection_, begin_row_, end_row_);
}

std::unique_ptr<MsRowProviderBase> MakeMsRowProvider(
    const std::string& ms_name, const MSSelection& selection,
    const std::map<size_t, size_t>& selected_data_description_ids,
    const std::string& data_column_name, bool require_model) {
  if (!casacore::Table::isReadable(ms_name)) {
    ThrowRuntimeError("The measurement set ", ms_name,
                      " can't be opened for reading.");
  }

  casacore::MeasurementSet ms(ms_name);

  // TODO AST-630 Return BdaMsRowProvider when BDA tables are present.

  return boost::make_unique<DirectMSRowProvider>(
      ms, selection, selected_data_description_ids, data_column_name,
      require_model);
}
