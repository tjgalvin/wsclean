#include "readinmemory.h"

#include <aocommon/logger.h>

#include "../msproviders/operations.h"
#include "../msproviders/rowproviders/msrowprovider.h"

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
  result.bands = aocommon::MultiBandData(ms.spectralWindow(), ms.dataDescription());

  const ObservationInfo info = ReadObservationInfo(ms, selection.FieldIds().front());
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
  casacore::ScalarColumn<double> interval_column(ms, ms.columnName(ms.INTERVAL));
  if (ms.nrow() == 0) {
    result.start_time = 0.0;
    result.interval = 0.0;
  } else {
    result.start_time = time_column(0);
    result.interval = interval_column(0);
  }
}

std::vector<InMemoryRow> ReadMsRows(casacore::MeasurementSet& ms,
                                    const std::string& data_column_name,
                                    const schaapcommon::reordering::MSSelection& selection) {
  casacore::ArrayColumn<std::complex<float>> data_column(ms, data_column_name);
  casacore::ArrayColumn<bool> flag_column(ms, ms.columnName(ms.FLAG));
  casacore::ArrayColumn<float> weight_column(ms, ms.columnName(ms.WEIGHT_SPECTRUM));

  casacore::ScalarColumn<int> antenna1_column(ms, ms.columnName(ms.ANTENNA1));
  casacore::ScalarColumn<int> antenna2_column(ms, ms.columnName(ms.ANTENNA1));
  casacore::ScalarColumn<int> data_desc_id_column(ms, ms.columnName(ms.DATA_DESC_ID));
  casacore::ScalarColumn<int> field_id_column(ms, ms.columnName(ms.FIELD_ID));
  casacore::ScalarColumn<double> time_column(ms, ms.columnName(ms.TIME));
  casacore::ArrayColumn<double> uvw_column(ms, ms.columnName(ms.UVW));

  casacore::Array<std::complex<float>> data_array;
  casacore::Array<float> weight_array;
  casacore::Array<bool> flag_array;
  casacore::Array<double> uvw(casacore::IPosition{3});

  size_t start_row;
  size_t end_row;
  GetRowRange(ms, selection, start_row, end_row);
  const bool require_time_ordering = selection.HasInterval();
  double time = time_column(start_row);
  size_t timestep = selection.HasInterval() ? selection.IntervalStart() : 0;

  std::vector<InMemoryRow> in_memory_rows;
  in_memory_rows.reserve(ms.nrow());
  for (size_t row_index = 0; row_index != ms.nrow(); ++row_index) {
    const int antenna1 = antenna1_column(row_index);
    const int antenna2 = antenna2_column(row_index);
    const int field_id = field_id_column(row_index);
    uvw_column.get(row_index, uvw, false);

    const double next_time = time_column(row_index);
    if (time != next_time) {
      ++timestep;
      if (require_time_ordering && time > next_time) ThrowTimeOrderException();
      time = next_time;
    }
    if (selection.IsSelected(field_id, timestep, antenna1, antenna2, uvw.data())) {
      data_column.get(row_index, data_array, true);
      weight_column.get(row_index, weight_array, true);
      flag_column.get(row_index, flag_array, true);

      InMemoryRow& row = in_memory_rows.emplace_back();
      row.data.assign(data_array.cbegin(), data_array.cend());
      row.weights.assign(weight_array.cbegin(), weight_array.cend());
      SetFlagsToZero(row.weights, std::span(flag_array.cbegin(), flag_array.cend()));

      row.antenna1 = antenna1;
      row.antenna2 = antenna2;
      row.data_desc_id = data_desc_id_column(row_index);
      row.field_id = field_id;
      row.time = time;
      std::copy_n(uvw.cbegin(), 3, row.uvw.begin());
    }
  }

  return in_memory_rows;
}

}  // namespace

InMemoryMs ReadMsInMemory(const std::string& measurement_set, const std::string& data_column_name,
                          const schaapcommon::reordering::MSSelection& selection,
                          bool sort_in_time) {
  casacore::MeasurementSet ms(measurement_set);
  if (sort_in_time) {
    aocommon::Logger::Debug << "Sorting measurement set...\n";
    ms = ms.sort(casacore::MS::columnName(casacore::MSMainEnums::TIME));
  }

  InMemoryMs result;

  ReadMsMetaData(ms, result, selection);

  result.rows = ReadMsRows(ms, data_column_name, selection);

  return result;
}

void ApplyRowSelectionToMs(InMemoryMs& ms, const schaapcommon::reordering::MSSelection& selection) {
  if (ms.rows.empty() || selection.IsUnset()) return;

  const bool require_time_ordering = selection.HasInterval();
  // By setting the previous time and timestep to -1, the first row will be
  // considered to be a "new" timestep, set time and increase timestep, such
  // that the first row is timestep 0.
  double time = -1.0;
  size_t timestep = static_cast<size_t>(-1);
  std::vector<InMemoryRow> selected_rows;

  for (InMemoryRow& row : ms.rows) {
    if (time != row.time) {
      ++timestep;
      if (require_time_ordering && time > row.time) ThrowTimeOrderException();
      time = row.time;
    }
    if (selection.IsSelected(row.field_id, timestep, row.antenna1, row.antenna2, row.uvw.data())) {
      selected_rows.emplace_back(std::move(row));
    }
  }
  ms.rows = std::move(selected_rows);
}

}  // namespace wsclean
