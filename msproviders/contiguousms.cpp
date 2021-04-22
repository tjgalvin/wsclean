#include "contiguousms.h"
#include "../io/logger.h"
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>

ContiguousMS::ContiguousMS(const string& msPath,
                           const std::string& dataColumnName,
                           const MSSelection& selection,
                           aocommon::PolarizationEnum polOut, size_t dataDescId)
    : _currentInputRow(0),
      _currentInputTimestep(0),
      _currentInputTime(0.0),
      _currentOutputRow(0),
      _currentOutputTimestep(0),
      _currentOutputTime(0.0),
      _currentRowId(0),
      _dataDescId(dataDescId),
      _nAntenna(0),
      _isModelColumnPrepared(false),
      _selection(selection),
      _polOut(polOut),
      _msPath(msPath),
      _dataColumnName(dataColumnName) {
  open();
}

void ContiguousMS::open() {
  Logger::Info << "Opening " << _msPath << ", spw " << _dataDescId
               << " with contiguous MS reader.\n";

  _ms = SynchronizedMS(_msPath);
  _antenna1Column = casacore::ScalarColumn<int>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA1));
  _antenna2Column = casacore::ScalarColumn<int>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA2));
  _fieldIdColumn = casacore::ScalarColumn<int>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::FIELD_ID));
  _dataDescIdColumn = casacore::ScalarColumn<int>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::DATA_DESC_ID));
  _timeColumn = casacore::ScalarColumn<double>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::TIME));
  _uvwColumn = casacore::ArrayColumn<double>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::UVW));
  _dataColumn = casacore::ArrayColumn<casacore::Complex>(*_ms, _dataColumnName);
  _flagColumn = casacore::ArrayColumn<bool>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::FLAG));

  _inputPolarizations = GetMSPolarizations(*_ms);

  const casacore::IPosition shape(_dataColumn.shape(0));
  _dataArray = casacore::Array<std::complex<float>>(shape);
  _weightSpectrumArray = casacore::Array<float>(shape);
  _imagingWeightSpectrumArray = casacore::Array<float>(shape);
  _flagArray = casacore::Array<bool>(shape);
  _bandData = MultiBandData(_ms->spectralWindow(), _ms->dataDescription());

  if (_bandData.BandCount() > 1) {
    throw std::runtime_error(
        "This set contains multiple spws, and can therefore not be opened "
        "directly due to possible synchronization issues between spws. You can "
        "force reordering of the measurement by adding -reorder to the command "
        "line.");
  }
  _nAntenna = _ms->antenna().nrow();

  _msHasWeightSpectrum =
      openWeightSpectrumColumn(*_ms, _weightSpectrumColumn, shape);
  if (!_msHasWeightSpectrum) {
    casacore::IPosition scalarShape(1, shape[0]);
    _weightScalarArray = casacore::Array<float>(scalarShape);
    _weightScalarColumn.reset(new casacore::ArrayColumn<float>(
        *_ms, casacore::MS::columnName(casacore::MSMainEnums::WEIGHT)));
  }

  getRowRangeAndIDMap(*_ms, _selection, _startRow, _endRow,
                      std::set<size_t>{size_t(_dataDescId)}, _idToMSRow);
  Reset();
}

void ContiguousMS::Reset() {
  _currentInputRow = _startRow - 1;
  _currentInputTime = 0.0;
  if (_selection.HasInterval())
    _currentInputTimestep = _selection.IntervalStart() - 1;
  else
    _currentInputTimestep = -1;
  _currentOutputRow = _currentInputRow;
  _currentOutputTime = _currentInputTime;
  _currentOutputTimestep = _currentInputTimestep;
  _currentRowId = size_t(-1);
  NextInputRow();
  NextOutputRow();
}

bool ContiguousMS::CurrentRowAvailable() {
  if (_currentInputRow >= _endRow) return false;

  int fieldId = _fieldIdColumn(_currentInputRow);
  int a1 = _antenna1Column(_currentInputRow);
  int a2 = _antenna2Column(_currentInputRow);
  int dataDescId = _dataDescIdColumn(_currentInputRow);
  casacore::Vector<double> uvw = _uvwColumn(_currentInputRow);

  while (!_selection.IsSelected(fieldId, _currentInputTimestep, a1, a2, uvw) ||
         dataDescId != _dataDescId) {
    ++_currentInputRow;
    if (_currentInputRow >= _endRow) return false;

    fieldId = _fieldIdColumn(_currentInputRow);
    a1 = _antenna1Column(_currentInputRow);
    a2 = _antenna2Column(_currentInputRow);
    uvw = _uvwColumn(_currentInputRow);
    dataDescId = _dataDescIdColumn(_currentInputRow);
    if (_currentInputTime != _timeColumn(_currentInputRow)) {
      ++_currentInputTimestep;
      _currentInputTime = _timeColumn(_currentInputRow);
    }

    _isDataRead = false;
    _isWeightRead = false;
    _isModelRead = false;
  }

  return true;
}

void ContiguousMS::NextInputRow() {
  _isDataRead = false;
  _isWeightRead = false;
  _isModelRead = false;

  ++_currentRowId;
  int fieldId, a1, a2, dataDescId;
  casacore::Vector<double> uvw;
  do {
    ++_currentInputRow;
    if (_currentInputRow >= _endRow) return;

    fieldId = _fieldIdColumn(_currentInputRow);
    a1 = _antenna1Column(_currentInputRow);
    a2 = _antenna2Column(_currentInputRow);
    uvw = _uvwColumn(_currentInputRow);
    dataDescId = _dataDescIdColumn(_currentInputRow);
    if (_currentInputTime != _timeColumn(_currentInputRow)) {
      ++_currentInputTimestep;
      _currentInputTime = _timeColumn(_currentInputRow);
    }
  } while (
      !_selection.IsSelected(fieldId, _currentInputTimestep, a1, a2, uvw) ||
      (dataDescId != _dataDescId));
}

void ContiguousMS::NextOutputRow() {
  int fieldId, a1, a2, dataDescId;
  casacore::Vector<double> uvw;
  do {
    ++_currentOutputRow;
    if (_currentOutputRow >= _endRow) return;

    fieldId = _fieldIdColumn(_currentOutputRow);
    a1 = _antenna1Column(_currentOutputRow);
    a2 = _antenna2Column(_currentOutputRow);
    uvw = _uvwColumn(_currentOutputRow);
    dataDescId = _dataDescIdColumn(_currentOutputRow);
    if (_currentOutputTime != _timeColumn(_currentOutputRow)) {
      ++_currentOutputTimestep;
      _currentOutputTime = _timeColumn(_currentOutputRow);
    }
  } while (
      !_selection.IsSelected(fieldId, _currentOutputTimestep, a1, a2, uvw) ||
      (dataDescId != _dataDescId));
}

double ContiguousMS::StartTime() {
  return casacore::MEpoch::ScalarColumn(
             *_ms, casacore::MS::columnName(casacore::MS::TIME))(_startRow)
      .getValue()
      .get();
}

void ContiguousMS::ReadMeta(double& u, double& v, double& w,
                            size_t& dataDescId) {
  casacore::Vector<double> uvwArray = _uvwColumn(_currentInputRow);
  u = uvwArray(0);
  v = uvwArray(1);
  w = uvwArray(2);
  dataDescId = _dataDescId;
}

void ContiguousMS::ReadMeta(MetaData& metaData) {
  casacore::Vector<double> uvwArray = _uvwColumn(_currentInputRow);
  metaData.uInM = uvwArray(0);
  metaData.vInM = uvwArray(1);
  metaData.wInM = uvwArray(2);
  metaData.dataDescId = _dataDescId;
  metaData.fieldId = _fieldIdColumn(_currentInputRow);
  metaData.antenna1 = _antenna1Column(_currentInputRow);
  metaData.antenna2 = _antenna2Column(_currentInputRow);
  metaData.time = _timeColumn(_currentInputRow);
}

void ContiguousMS::ReadData(std::complex<float>* buffer) {
  readData();
  readWeights();
  size_t startChannel, endChannel;
  if (_selection.HasChannelRange()) {
    startChannel = _selection.ChannelRangeStart();
    endChannel = _selection.ChannelRangeEnd();
  } else {
    startChannel = 0;
    endChannel = _bandData[_dataDescId].ChannelCount();
  }
  copyData(buffer, startChannel, endChannel, _inputPolarizations, _dataArray,
           _polOut);
}

size_t ContiguousMS::NChannels() {
  if (_selection.HasChannelRange())
    return _selection.ChannelRangeEnd() - _selection.ChannelRangeStart();
  else
    return _bandData[_dataDescId].ChannelCount();
}

size_t ContiguousMS::NPolarizations() { return _inputPolarizations.size(); }

void ContiguousMS::prepareModelColumn() {
  initializeModelColumn(*_ms);

  _modelColumn.reset(new casacore::ArrayColumn<casacore::Complex>(
      *_ms, casacore::MS::columnName(casacore::MSMainEnums::MODEL_DATA)));
  const casacore::IPosition shape(_modelColumn->shape(0));
  _modelArray = casacore::Array<std::complex<float>>(shape);
  _isModelColumnPrepared = true;
}

void ContiguousMS::ReadModel(std::complex<float>* buffer) {
  if (!_isModelColumnPrepared) prepareModelColumn();

  readModel();
  readWeights();
  size_t startChannel, endChannel;
  if (_selection.HasChannelRange()) {
    startChannel = _selection.ChannelRangeStart();
    endChannel = _selection.ChannelRangeEnd();
  } else {
    startChannel = 0;
    endChannel = _bandData[_dataDescId].ChannelCount();
  }
  copyData(buffer, startChannel, endChannel, _inputPolarizations, _modelArray,
           _polOut);
}

void ContiguousMS::WriteModel(const std::complex<float>* buffer, bool addToMS) {
  if (!_isModelColumnPrepared) prepareModelColumn();

  size_t startChannel, endChannel;
  if (_selection.HasChannelRange()) {
    startChannel = _selection.ChannelRangeStart();
    endChannel = _selection.ChannelRangeEnd();
  } else {
    startChannel = 0;
    endChannel = _bandData[_dataDescId].ChannelCount();
  }

  _modelColumn->get(_currentOutputRow, _modelArray);
  if (addToMS) {
    reverseCopyData<true>(_modelArray, startChannel, endChannel,
                          _inputPolarizations, buffer, _polOut);
  } else {
    reverseCopyData<false>(_modelArray, startChannel, endChannel,
                           _inputPolarizations, buffer, _polOut);
  }
  _modelColumn->put(_currentOutputRow, _modelArray);
}

void ContiguousMS::ReadWeights(std::complex<float>* buffer) {
  readData();
  readWeights();
  size_t startChannel, endChannel;
  if (_selection.HasChannelRange()) {
    startChannel = _selection.ChannelRangeStart();
    endChannel = _selection.ChannelRangeEnd();
  } else {
    startChannel = 0;
    endChannel = _bandData[_dataDescId].ChannelCount();
  }
  copyWeights(buffer, startChannel, endChannel, _inputPolarizations, _dataArray,
              _weightSpectrumArray, _flagArray, _polOut);
}

void ContiguousMS::ReadWeights(float* buffer) {
  readData();
  readWeights();
  size_t startChannel, endChannel;
  if (_selection.HasChannelRange()) {
    startChannel = _selection.ChannelRangeStart();
    endChannel = _selection.ChannelRangeEnd();
  } else {
    startChannel = 0;
    endChannel = _bandData[_dataDescId].ChannelCount();
  }
  copyWeights(buffer, startChannel, endChannel, _inputPolarizations, _dataArray,
              _weightSpectrumArray, _flagArray, _polOut);
}

void ContiguousMS::WriteImagingWeights(const float* buffer) {
  if (_imagingWeightsColumn == nullptr) {
    _imagingWeightsColumn.reset(
        new casacore::ArrayColumn<float>(initializeImagingWeightColumn(*_ms)));
  }
  size_t msRowId = _idToMSRow[_currentInputRow];
  size_t dataDescId = _dataDescIdColumn(msRowId);
  size_t startChannel, endChannel;
  if (_selection.HasChannelRange()) {
    startChannel = _selection.ChannelRangeStart();
    endChannel = _selection.ChannelRangeEnd();
  } else {
    startChannel = 0;
    endChannel = _bandData[dataDescId].ChannelCount();
  }

  _imagingWeightsColumn->get(msRowId, _imagingWeightSpectrumArray);
  reverseCopyWeights(_imagingWeightSpectrumArray, startChannel, endChannel,
                     _inputPolarizations, buffer, _polOut);
  _imagingWeightsColumn->put(msRowId, _imagingWeightSpectrumArray);
}

void ContiguousMS::MakeIdToMSRowMapping(std::vector<size_t>& idToMSRow) {
  idToMSRow = _idToMSRow;
}
