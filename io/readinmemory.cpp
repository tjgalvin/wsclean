#include "readinmemory.h"

#include "../msproviders/operations.h"

using schaapcommon::reordering::MSSelection;

namespace wsclean {
namespace {

void SetFlagsToZero(std::span<float> values, std::span<const bool> flags) {
  assert(values.size() == flags.size());
  std::span<float>::iterator value_iter = values.begin();
  for (bool flag : flags) {
    if (flag) *value_iter = 0.0f;
    ++value_iter;
  }
}

void ReadMsMetaData(casacore::MeasurementSet& ms, InMemoryMs& result,
                    const MSSelection& selection) {
  result.antenna_names = wsclean::GetAntennaNames(ms.antenna());
  result.bands =
      aocommon::MultiBandData(ms.spectralWindow(), ms.dataDescription());

  const ObservationInfo info =
      ReadObservationInfo(ms, selection.FieldIds().front());
  result.field_name = info.fieldName;
  result.observer = info.observer;
  result.phase_centre_dec = info.phaseCentreDec;
  result.phase_centre_ra = info.phaseCentreRA;
  const std::set<aocommon::PolarizationEnum> polarizations =
      GetMSPolarizations(*result.bands.DataDescIds().begin(), ms);
  result.polarizations.assign(polarizations.begin(), polarizations.end());
  result.telescope_name = info.telescopeName;

  result.has_frequency_bda = HasFrequencyBda(ms);
  result.has_time_bda = false;

  casacore::ScalarColumn<double> time_column(ms, ms.columnName(ms.TIME));
  casacore::ScalarColumn<double> interval_column(ms,
                                                 ms.columnName(ms.INTERVAL));
  if (ms.nrow() == 0) {
    result.start_time = 0.0;
    result.interval = 0.0;
  } else {
    result.start_time = time_column(0);
    result.interval = interval_column(0);
  }
}

std::vector<InMemoryRow> ReadMsRows(
    casacore::MeasurementSet& ms, const std::string& data_column_name,
    const schaapcommon::reordering::MSSelection& selection) {
  casacore::ArrayColumn<std::complex<float>> data_column(ms, data_column_name);
  casacore::ArrayColumn<bool> flag_column(ms, ms.columnName(ms.FLAG));
  casacore::ArrayColumn<float> weight_column(ms,
                                             ms.columnName(ms.WEIGHT_SPECTRUM));

  casacore::ScalarColumn<int> antenna1_column(ms, ms.columnName(ms.ANTENNA1));
  casacore::ScalarColumn<int> antenna2_column(ms, ms.columnName(ms.ANTENNA1));
  casacore::ScalarColumn<int> data_desc_id_column(
      ms, ms.columnName(ms.DATA_DESC_ID));
  casacore::ScalarColumn<int> field_id_column(ms, ms.columnName(ms.FIELD_ID));
  casacore::ScalarColumn<double> time_column(ms, ms.columnName(ms.TIME));
  casacore::ArrayColumn<double> uvw_column(ms, ms.columnName(ms.UVW));

  casacore::Array<std::complex<float>> data_array;
  casacore::Array<float> weight_array;
  casacore::Array<bool> flag_array;
  casacore::Array<double> uvw(casacore::IPosition{3});

  std::vector<InMemoryRow> in_memory_rows;
  in_memory_rows.reserve(ms.nrow());
  for (size_t row_index = 0; row_index != ms.nrow(); ++row_index) {
    data_column.get(row_index, data_array, true);
    weight_column.get(row_index, weight_array, true);
    flag_column.get(row_index, flag_array, true);
    uvw_column.get(row_index, uvw, false);

    InMemoryRow& row = in_memory_rows.emplace_back();
    row.data.assign(data_array.cbegin(), data_array.cend());
    row.weights.assign(weight_array.cbegin(), weight_array.cend());
    SetFlagsToZero(row.weights,
                   std::span(flag_array.cbegin(), flag_array.cend()));

    row.antenna1 = antenna1_column(row_index);
    row.antenna2 = antenna2_column(row_index);
    row.data_desc_id = data_desc_id_column(row_index);
    row.field_id = field_id_column(row_index);
    row.time = time_column(row_index);
    std::copy_n(uvw.cbegin(), 3, row.uvw.begin());
  }

  return in_memory_rows;
}

}  // namespace

InMemoryMs ReadMsInMemory(
    const std::string& measurement_set, const std::string& data_column_name,
    const schaapcommon::reordering::MSSelection& selection) {
  casacore::MeasurementSet ms(measurement_set);

  InMemoryMs result;

  ReadMsMetaData(ms, result, selection);

  result.rows = ReadMsRows(ms, data_column_name, selection);

  return result;
}

}  // namespace wsclean
