#include "mscolumns.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

MsColumns::MsColumns(casacore::MeasurementSet& ms,
                     const std::string& data_column_name)
    : antenna_1(ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA1)),
      antenna_2(ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA2)),
      field_id(ms, casacore::MS::columnName(casacore::MSMainEnums::FIELD_ID)),
      time(ms, casacore::MS::columnName(casacore::MSMainEnums::TIME)),
      epoch_as_time(ms, casacore::MS::columnName(casacore::MSMainEnums::TIME)),
      uvw(ms, casacore::MS::columnName(casacore::MSMainEnums::UVW)),
      data(ms, data_column_name),
      flag(ms, casacore::MS::columnName(casacore::MSMainEnums::FLAG)),
      data_description_id(
          ms, casacore::MS::columnName(casacore::MSMainEnums::DATA_DESC_ID)) {}
