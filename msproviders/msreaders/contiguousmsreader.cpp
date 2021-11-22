#include "contiguousmsreader.h"
#include "../contiguousms.h"

ContiguousMSReader::ContiguousMSReader(ContiguousMS* contiguousms)
    : MSReader(contiguousms),
      _currentInputRow(contiguousms->_startRow - 1),
      _currentInputTimestep(size_t(-1)),
      _currentInputTime(0.0),
      _currentRowId(size_t(-1)),
      _isDataRead(false),
      _isModelRead(false),
      _isWeightRead(false),
      _imagingWeightsColumn() {
  if (contiguousms->_selection.HasInterval())
    _currentInputTimestep = contiguousms->_selection.IntervalStart() - 1;
  NextInputRow();
}

bool ContiguousMSReader::CurrentRowAvailable() {
  const ContiguousMS& contiguousms =
      static_cast<const ContiguousMS&>(*_msProvider);

  if (_currentInputRow >= contiguousms._endRow) return false;

  int fieldId = contiguousms._fieldIdColumn(_currentInputRow);
  int a1 = contiguousms._antenna1Column(_currentInputRow);
  int a2 = contiguousms._antenna2Column(_currentInputRow);
  int dataDescId = contiguousms._dataDescIdColumn(_currentInputRow);
  casacore::Vector<double> uvw = contiguousms._uvwColumn(_currentInputRow);

  while (!contiguousms._selection.IsSelected(fieldId, _currentInputTimestep, a1,
                                             a2, uvw) ||
         dataDescId != contiguousms._dataDescId) {
    ++_currentInputRow;
    if (_currentInputRow >= contiguousms._endRow) return false;

    fieldId = contiguousms._fieldIdColumn(_currentInputRow);
    a1 = contiguousms._antenna1Column(_currentInputRow);
    a2 = contiguousms._antenna2Column(_currentInputRow);
    uvw = contiguousms._uvwColumn(_currentInputRow);
    dataDescId = contiguousms._dataDescIdColumn(_currentInputRow);
    if (_currentInputTime != contiguousms._timeColumn(_currentInputRow)) {
      ++_currentInputTimestep;
      _currentInputTime = contiguousms._timeColumn(_currentInputRow);
    }

    _isDataRead = false;
    _isWeightRead = false;
    _isModelRead = false;
  }

  return true;
}

void ContiguousMSReader::NextInputRow() {
  const ContiguousMS& contiguousms =
      static_cast<const ContiguousMS&>(*_msProvider);

  _isDataRead = false;
  _isWeightRead = false;
  _isModelRead = false;

  ++_currentRowId;
  int fieldId, a1, a2, dataDescId;
  casacore::Vector<double> uvw;

  do {
    ++_currentInputRow;
    if (_currentInputRow >= contiguousms._endRow) return;

    fieldId = contiguousms._fieldIdColumn(_currentInputRow);
    a1 = contiguousms._antenna1Column(_currentInputRow);
    a2 = contiguousms._antenna2Column(_currentInputRow);
    uvw = contiguousms._uvwColumn(_currentInputRow);
    dataDescId = contiguousms._dataDescIdColumn(_currentInputRow);
    if (_currentInputTime != contiguousms._timeColumn(_currentInputRow)) {
      ++_currentInputTimestep;
      _currentInputTime = contiguousms._timeColumn(_currentInputRow);
    }
  } while (!contiguousms._selection.IsSelected(fieldId, _currentInputTimestep,
                                               a1, a2, uvw) ||
           (dataDescId != contiguousms._dataDescId));
}

void ContiguousMSReader::ReadMeta(double& u, double& v, double& w) {
  const ContiguousMS& contiguousms =
      static_cast<const ContiguousMS&>(*_msProvider);

  casacore::Vector<double> uvwArray = contiguousms._uvwColumn(_currentInputRow);
  u = uvwArray(0);
  v = uvwArray(1);
  w = uvwArray(2);
}

void ContiguousMSReader::ReadMeta(MSProvider::MetaData& metaData) {
  const ContiguousMS& contiguousms =
      static_cast<const ContiguousMS&>(*_msProvider);

  casacore::Vector<double> uvwArray = contiguousms._uvwColumn(_currentInputRow);
  metaData.uInM = uvwArray(0);
  metaData.vInM = uvwArray(1);
  metaData.wInM = uvwArray(2);
  metaData.fieldId = contiguousms._fieldIdColumn(_currentInputRow);
  metaData.antenna1 = contiguousms._antenna1Column(_currentInputRow);
  metaData.antenna2 = contiguousms._antenna2Column(_currentInputRow);
  metaData.time = contiguousms._timeColumn(_currentInputRow);
}

void ContiguousMSReader::ReadData(std::complex<float>* buffer) {
  ContiguousMS& contiguousms = static_cast<ContiguousMS&>(*_msProvider);

  readData();
  readWeights();
  size_t startChannel, endChannel;
  if (contiguousms._selection.HasChannelRange()) {
    startChannel = contiguousms._selection.ChannelRangeStart();
    endChannel = contiguousms._selection.ChannelRangeEnd();
  } else {
    startChannel = 0;
    endChannel =
        contiguousms._bandData[contiguousms._dataDescId].ChannelCount();
  }
  copyData(buffer, startChannel, endChannel, contiguousms._inputPolarizations,
           contiguousms._dataArray, contiguousms._polOut);
}

void ContiguousMSReader::ReadModel(std::complex<float>* buffer) {
  ContiguousMS& contiguousms = static_cast<ContiguousMS&>(*_msProvider);

  if (!contiguousms._isModelColumnPrepared) contiguousms.prepareModelColumn();

  readModel();
  readWeights();
  size_t startChannel, endChannel;
  if (contiguousms._selection.HasChannelRange()) {
    startChannel = contiguousms._selection.ChannelRangeStart();
    endChannel = contiguousms._selection.ChannelRangeEnd();
  } else {
    startChannel = 0;
    endChannel =
        contiguousms._bandData[contiguousms._dataDescId].ChannelCount();
  }
  copyData(buffer, startChannel, endChannel, contiguousms._inputPolarizations,
           contiguousms._modelArray, contiguousms._polOut);
}

void ContiguousMSReader::ReadWeights(std::complex<float>* buffer) {
  const ContiguousMS& contiguousms =
      static_cast<const ContiguousMS&>(*_msProvider);

  readData();
  readWeights();
  size_t startChannel, endChannel;
  if (contiguousms._selection.HasChannelRange()) {
    startChannel = contiguousms._selection.ChannelRangeStart();
    endChannel = contiguousms._selection.ChannelRangeEnd();
  } else {
    startChannel = 0;
    endChannel =
        contiguousms._bandData[contiguousms._dataDescId].ChannelCount();
  }
  MSProvider::CopyWeights(
      buffer, startChannel, endChannel, contiguousms._inputPolarizations,
      contiguousms._dataArray, contiguousms._weightSpectrumArray,
      contiguousms._flagArray, contiguousms._polOut);
}

void ContiguousMSReader::ReadWeights(float* buffer) {
  const ContiguousMS& contiguousms =
      static_cast<const ContiguousMS&>(*_msProvider);

  readData();
  readWeights();
  size_t startChannel, endChannel;
  if (contiguousms._selection.HasChannelRange()) {
    startChannel = contiguousms._selection.ChannelRangeStart();
    endChannel = contiguousms._selection.ChannelRangeEnd();
  } else {
    startChannel = 0;
    endChannel =
        contiguousms._bandData[contiguousms._dataDescId].ChannelCount();
  }
  MSProvider::CopyWeights(
      buffer, startChannel, endChannel, contiguousms._inputPolarizations,
      contiguousms._dataArray, contiguousms._weightSpectrumArray,
      contiguousms._flagArray, contiguousms._polOut);
}

void ContiguousMSReader::WriteImagingWeights(const float* buffer) {
  ContiguousMS& contiguousms = static_cast<ContiguousMS&>(*_msProvider);

  if (_imagingWeightsColumn == nullptr) {
    _imagingWeightsColumn.reset(new casacore::ArrayColumn<float>(
        MSProvider::InitializeImagingWeightColumn(*(contiguousms._ms))));
  }
  size_t dataDescId = contiguousms._dataDescIdColumn(_currentInputRow);
  size_t startChannel, endChannel;
  if (contiguousms._selection.HasChannelRange()) {
    startChannel = contiguousms._selection.ChannelRangeStart();
    endChannel = contiguousms._selection.ChannelRangeEnd();
  } else {
    startChannel = 0;
    endChannel = contiguousms._bandData[dataDescId].ChannelCount();
  }

  _imagingWeightsColumn->get(_currentInputRow,
                             contiguousms._imagingWeightSpectrumArray);
  MSProvider::ReverseCopyWeights(
      contiguousms._imagingWeightSpectrumArray, startChannel, endChannel,
      contiguousms._inputPolarizations, buffer, contiguousms._polOut);
  _imagingWeightsColumn->put(_currentInputRow,
                             contiguousms._imagingWeightSpectrumArray);
}

void ContiguousMSReader::readData() {
  ContiguousMS& contiguousms = static_cast<ContiguousMS&>(*_msProvider);
  if (!_isDataRead) {
    contiguousms._dataColumn.get(_currentInputRow, contiguousms._dataArray);
    _isDataRead = true;
  }
}

void ContiguousMSReader::readWeights() {
  ContiguousMS& contiguousms = static_cast<ContiguousMS&>(*_msProvider);

  if (!_isWeightRead) {
    contiguousms._flagColumn.get(_currentInputRow, contiguousms._flagArray);
    if (contiguousms._msHasWeightSpectrum)
      contiguousms._weightSpectrumColumn->get(
          _currentInputRow, contiguousms._weightSpectrumArray);
    else {
      contiguousms._weightScalarColumn->get(_currentInputRow,
                                            contiguousms._weightScalarArray);
      contiguousms.ExpandScalarWeights(contiguousms._weightScalarArray,
                                       contiguousms._weightSpectrumArray);
    }
    _isWeightRead = true;
  }
}

void ContiguousMSReader::readModel() {
  ContiguousMS& contiguousms = static_cast<ContiguousMS&>(*_msProvider);
  if (!_isModelRead) {
    contiguousms._modelColumn.get(_currentInputRow, contiguousms._modelArray);
    _isModelRead = true;
  }
}
