#ifndef MSPROVIDERS_MSCOLUMNS_H
#define MSPROVIDERS_MSCOLUMNS_H

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>

#include <string>

namespace casacore {
class MeasurementSet;
}  // namespace casacore

/** Helper class to provide the columns in a measurement set. */
struct MsColumns {
  explicit MsColumns(casacore::MeasurementSet& ms,
                     const std::string& data_column_name);

  casacore::ScalarColumn<int> antenna_1;
  casacore::ScalarColumn<int> antenna_2;
  casacore::ScalarColumn<int> field_id;
  casacore::ScalarColumn<double> time;
  casacore::MEpoch::ScalarColumn epoch_as_time;
  casacore::ArrayColumn<double> uvw;
  casacore::ArrayColumn<casacore::Complex> data;
  casacore::ArrayColumn<bool> flag;
  casacore::ScalarColumn<int> data_description_id;
};

#endif  //  MSPROVIDERS_MSCOLUMNS_H
