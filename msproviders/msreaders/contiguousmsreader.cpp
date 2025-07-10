#include "contiguousmsreader.h"
#include "../contiguousms.h"

#include <schaapcommon/reordering/reordering.h>

namespace wsclean {

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
      static_cast<const ContiguousMS&>(*ms_provider_);

  if (_currentInputRow >= contiguousms._endRow) return false;

  int fieldId = contiguousms._fieldIdColumn(_currentInputRow);
  int a1 = contiguousms._antenna1Column(_currentInputRow);
  int a2 = contiguousms._antenna2Column(_currentInputRow);
  int dataDescId = contiguousms._dataDescIdColumn(_currentInputRow);
  casacore::Vector<double> uvw = contiguousms._uvwColumn(_currentInputRow);

  while (!contiguousms._selection.IsSelected(fieldId, _currentInputTimestep, a1,
                                             a2, uvw.data()) ||
         dataDescId != contiguousms._dataDescId) {
    ++_currentInputRow;
    if (_currentInputRow >= contiguousms._endRow) return false;

    fieldId = contiguousms._fieldIdColumn(_currentInputRow);
    a1 = contiguousms._antenna1Column(_currentInputRow);
    a2 = contiguousms._antenna2Column(_currentInputRow);
    contiguousms._uvwColumn.get(_currentInputRow, uvw, true);
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
      static_cast<const ContiguousMS&>(*ms_provider_);

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
                                               a1, a2, uvw.data()) ||
           (dataDescId != contiguousms._dataDescId));
}

void ContiguousMSReader::ReadMeta(MSProvider::MetaData& metadata) {
  const ContiguousMS& contiguousms =
      static_cast<const ContiguousMS&>(*ms_provider_);

  casacore::Vector<double> uvwArray = contiguousms._uvwColumn(_currentInputRow);
  metadata.u_in_m = uvwArray(0);
  metadata.v_in_m = uvwArray(1);
  metadata.w_in_m = uvwArray(2);
  metadata.time = contiguousms._timeColumn(_currentInputRow);
  metadata.data_desc_id = contiguousms._dataDescId;
  metadata.field_id = contiguousms._fieldIdColumn(_currentInputRow);
  metadata.antenna1 = contiguousms._antenna1Column(_currentInputRow);
  metadata.antenna2 = contiguousms._antenna2Column(_currentInputRow);
}

void ContiguousMSReader::ReadData(std::complex<float>* buffer) {
  ContiguousMS& contiguousms = static_cast<ContiguousMS&>(*ms_provider_);

  readData();
  readWeights();
  size_t startChannel, endChannel;
  if (contiguousms._selection.HasChannelRange()) {
    startChannel = contiguousms._selection.ChannelRangeStart();
    endChannel = contiguousms._selection.ChannelRangeEnd();
  } else {
    startChannel = 0;
    endChannel =
        contiguousms.original_bands_[contiguousms._dataDescId].ChannelCount();
  }
  schaapcommon::reordering::ExtractData(
      buffer, startChannel, endChannel, contiguousms._inputPolarizations,
      contiguousms._dataArray.data(), contiguousms._outputPolarization);
}

void ContiguousMSReader::ReadModel(std::complex<float>* buffer) {
  ContiguousMS& contiguousms = static_cast<ContiguousMS&>(*ms_provider_);

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
        contiguousms.original_bands_[contiguousms._dataDescId].ChannelCount();
  }
  schaapcommon::reordering::ExtractData(
      buffer, startChannel, endChannel, contiguousms._inputPolarizations,
      contiguousms._modelArray.data(), contiguousms._outputPolarization);
}

void ContiguousMSReader::ReadWeights(float* buffer) {
  const ContiguousMS& contiguousms =
      static_cast<const ContiguousMS&>(*ms_provider_);

  readData();
  readWeights();
  size_t startChannel, endChannel;
  if (contiguousms._selection.HasChannelRange()) {
    startChannel = contiguousms._selection.ChannelRangeStart();
    endChannel = contiguousms._selection.ChannelRangeEnd();
  } else {
    startChannel = 0;
    endChannel =
        contiguousms.original_bands_[contiguousms._dataDescId].ChannelCount();
  }
  schaapcommon::reordering::ExtractWeights(
      buffer, startChannel, endChannel, contiguousms._inputPolarizations,
      contiguousms._dataArray.data(), contiguousms._weightSpectrumArray.data(),
      contiguousms._flagArray.data(), contiguousms._outputPolarization);
}

void ContiguousMSReader::WriteImagingWeights(const float* buffer) {
  ContiguousMS& contiguousms = static_cast<ContiguousMS&>(*ms_provider_);

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
    endChannel = contiguousms.original_bands_[dataDescId].ChannelCount();
  }

  _imagingWeightsColumn->get(_currentInputRow,
                             contiguousms._imagingWeightSpectrumArray);
  schaapcommon::reordering::StoreWeights(
      contiguousms._imagingWeightSpectrumArray.data(), startChannel, endChannel,
      contiguousms._inputPolarizations, buffer,
      contiguousms._outputPolarization);
  _imagingWeightsColumn->put(_currentInputRow,
                             contiguousms._imagingWeightSpectrumArray);
}

void ContiguousMSReader::readData() {
  ContiguousMS& contiguousms = static_cast<ContiguousMS&>(*ms_provider_);
  if (!_isDataRead) {
    contiguousms._dataColumn.get(_currentInputRow, contiguousms._dataArray);
    _isDataRead = true;
  }
}

void ContiguousMSReader::readWeights() {
  ContiguousMS& contiguousms = static_cast<ContiguousMS&>(*ms_provider_);

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
  ContiguousMS& contiguousms = static_cast<ContiguousMS&>(*ms_provider_);
  if (!_isModelRead) {
    contiguousms._modelColumn.get(_currentInputRow, contiguousms._modelArray);
    _isModelRead = true;
  }
}

}  // namespace wsclean
