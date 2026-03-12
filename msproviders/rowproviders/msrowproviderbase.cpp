#include "msrowproviderbase.h"

#include "bdamsrowprovider.h"
#include "directmsrowprovider.h"

#include "../operations.h"

#include <aocommon/logger.h>
#include <aocommon/throwruntimeerror.h>

#include <casacore/tables/Tables/TableRecord.h>

#include <memory>

namespace wsclean {

void ThrowTimeOrderException() {
  throw std::runtime_error(
      "The data in this measurement set are not ordered in time. This is "
      "required for interval selection. To make wsclean order the data in "
      "time while reading, add -sort-ms to the command line. This is not the "
      "default because it is slightly slower.");
}

MsRowProviderBase::MsRowProviderBase(const casacore::MeasurementSet& ms,
                                     const schaapcommon::reordering::MSSelection& selection,
                                     const std::string& data_column_name,
                                     const std::string& model_column_name)
    : ms_(ms),
      selection_(selection),
      columns_(ms_, data_column_name),
      model_column_name_(model_column_name) {
  GetRowRange(ms_, selection_, begin_row_, end_row_);
}

std::unique_ptr<MsRowProviderBase> MakeMsRowProvider(
    const std::string& ms_name, const schaapcommon::reordering::MSSelection& selection,
    const std::map<size_t, size_t>& selected_data_description_ids,
    const std::string& data_column_name, const std::string& model_column_name, bool require_model,
    bool sort_in_time) {
  if (!casacore::Table::isReadable(ms_name)) {
    aocommon::ThrowRuntimeError("The measurement set ", ms_name, " can't be opened for reading.");
  }

  casacore::MeasurementSet ms(ms_name);
  if (sort_in_time) {
    aocommon::Logger::Debug << "Sorting measurement set...\n";
    ms = ms.sort(casacore::MS::columnName(casacore::MSMainEnums::TIME));
  }
  if (MsHasBdaData(ms))
    return std::make_unique<BdaMsRowProvider>(ms, selection, selected_data_description_ids,
                                              data_column_name, model_column_name, require_model);

  return std::make_unique<DirectMSRowProvider>(ms, selection, selected_data_description_ids,
                                               data_column_name, model_column_name, require_model);
}

bool MsHasBdaData(const casacore::MeasurementSet& ms) {
  return ms.keywordSet().isDefined(BdaMsRowProvider::kBdaFactorsTable) &&
         ms.keywordSet().asTable(BdaMsRowProvider::kBdaFactorsTable).nrow() != 0;
}

double GetBdaMaxTimeInterval(const casacore::MeasurementSet& ms) {
  if (!ms.keywordSet().isDefined(BdaMsRowProvider::kBdaFactorsTable))
    throw std::runtime_error("To process BDA data with WSClean, the " +
                             BdaMsRowProvider::kBdaFactorsTable +
                             " table needs to be present in the measurement set");
  casacore::Table table(ms.keywordSet().asTable(BdaMsRowProvider::kBdaTimeAxisTable));
  if (table.nrow() != 1)
    throw std::runtime_error(
        "Measurement set with BDA data and multiple fields or spectral "
        "windows is not supported");
  casacore::ScalarColumn<double> max_time_interval_column(
      table, BdaMsRowProvider::kBdaMaxTimeIntervalColumn);
  return max_time_interval_column(0);
}

}  // namespace wsclean
