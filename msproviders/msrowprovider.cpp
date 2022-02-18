#include "directmsrowprovider.h"
#include "msprovider.h"

#include <aocommon/throwruntimeerror.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

MSRowProvider::MSRowProvider(
    const string& msPath, const MSSelection& selection,
    const std::map<size_t, size_t>& selectedDataDescIds,
    const std::string& dataColumnName, bool requireModel)
    : MsRowProviderBase(casacore::MeasurementSet(msPath), selection,
                        dataColumnName),
      _selectedDataDescIds(selectedDataDescIds),
      _requireModel(requireModel) {
  Initialize();
}

MSRowProvider::MSRowProvider(
    const casacore::MeasurementSet& ms, const MSSelection& selection,
    const std::map<size_t, size_t>& selected_data_description_ids,
    const std::string& data_column_name, bool require_model)
    : MsRowProviderBase(ms, selection, data_column_name),
      _selectedDataDescIds(selected_data_description_ids),
      _requireModel(require_model) {
  Initialize();
}

void MSRowProvider::Initialize() {
  if (MsHasBdaData(Ms()))
    aocommon::ThrowRuntimeError(
        "Measurement set contains BDA data, but isn't opened for BDA "
        "processing.");

  if (_requireModel)
    _modelColumn.reset(new casacore::ArrayColumn<casacore::Complex>(
        Ms(), casacore::MS::columnName(casacore::MSMainEnums::MODEL_DATA)));

  MsColumns& columns = Columns();
  _msHasWeights =
      MSProvider::OpenWeightSpectrumColumn(Ms(), _weightSpectrumColumn);
  if (!_msHasWeights) {
    const size_t nCorrelations = columns.data.shape(0)[0];
    const casacore::IPosition scalarShape(1, nCorrelations);
    _scratchWeightScalarArray = casacore::Array<float>(scalarShape);
    _weightScalarColumn.reset(new casacore::ArrayColumn<float>(
        Ms(), casacore::MS::columnName(casacore::MSMainEnums::WEIGHT)));
  }

  // Determine last timestep
  _startTimestep = Selection().HasInterval() ? Selection().IntervalStart() : 0;
  size_t timestep = _startTimestep;
  double time = columns.time(BeginRow());
  for (size_t row = BeginRow(); row != EndRow(); ++row) {
    if (time != columns.time(row)) {
      ++timestep;
      time = columns.time(row);
    }
  }
  _endTimestep = timestep;

  _currentRow = BeginRow();
  _currentTimestep = _startTimestep;
  _currentTime = columns.time(BeginRow());
  _currentUVWArray = columns.uvw(BeginRow());
  _currentDataDescId = columns.data_description_id(BeginRow());

  // If this row is not selected, it is necessary to continue to the first
  // selected row.
  const int a1 = columns.antenna_1(_currentRow);
  const int a2 = columns.antenna_2(_currentRow);
  const int fieldId = columns.field_id(_currentRow);
  if (!isCurrentRowSelected(fieldId, a1, a2)) MSRowProvider::NextRow();
}

void MSRowProvider::NextRow() {
  bool isRowSelected;
  do {
    ++_currentRow;
    if (_currentRow == EndRow()) {
      break;
    } else {
      MsColumns& columns = Columns();
      const int a1 = columns.antenna_1(_currentRow);
      const int a2 = columns.antenna_2(_currentRow);
      const int fieldId = columns.field_id(_currentRow);
      _currentDataDescId = columns.data_description_id(_currentRow);

      if (_currentTime != columns.time(_currentRow)) {
        ++_currentTimestep;
        _currentTime = columns.time(_currentRow);
      }
      _currentUVWArray = columns.uvw(_currentRow);
      isRowSelected = isCurrentRowSelected(fieldId, a1, a2);
    }
  } while (!isRowSelected);
}

bool MSRowProvider::isCurrentRowSelected(int fieldId, int a1, int a2) const {
  std::map<size_t, size_t>::const_iterator dataDescIdIter =
      _selectedDataDescIds.find(_currentDataDescId);
  bool isDataDescIdSelected = dataDescIdIter != _selectedDataDescIds.end();
  return _selection.IsSelected(fieldId, _currentTimestep, a1, a2,
                               _currentUVWArray) &&
         isDataDescIdSelected;
}

void MSRowProvider::getCurrentWeights(WeightArray& weights,
                                      const casacore::IPosition& shape) {
  if (_msHasWeights)
    _weightSpectrumColumn->get(_currentRow, weights, true);
  else {
    _weightScalarColumn->get(_currentRow, _scratchWeightScalarArray);
    weights.resize(shape);
    MSProvider::ExpandScalarWeights(_scratchWeightScalarArray, weights);
  }
}
