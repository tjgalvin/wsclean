#include "contiguousmsreader.h"

#include "../contiguousms.h"
#include "../operations.h"

#include <schaapcommon/reordering/reordering.h>

using schaapcommon::reordering::ContainsDataDescId;

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
  if (contiguousms->selection_.HasInterval())
    _currentInputTimestep = contiguousms->selection_.IntervalStart() - 1;
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

  while (!contiguousms.selection_.IsSelected(fieldId, _currentInputTimestep, a1,
                                             a2, uvw.data()) ||
         !ContainsDataDescId(contiguousms.channel_ranges_, dataDescId)) {
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
  } while (!contiguousms.selection_.IsSelected(fieldId, _currentInputTimestep,
                                               a1, a2, uvw.data()) ||
           !ContainsDataDescId(contiguousms.channel_ranges_, dataDescId));
}

void ContiguousMSReader::ReadMeta(MSProvider::MetaData& metadata) {
  const ContiguousMS& contiguousms =
      static_cast<const ContiguousMS&>(*ms_provider_);

  casacore::Vector<double> uvwArray = contiguousms._uvwColumn(_currentInputRow);
  metadata.u_in_m = uvwArray(0);
  metadata.v_in_m = uvwArray(1);
  metadata.w_in_m = uvwArray(2);
  metadata.time = contiguousms._timeColumn(_currentInputRow);
  metadata.data_desc_id = contiguousms._dataDescIdColumn(_currentInputRow);
  metadata.field_id = contiguousms._fieldIdColumn(_currentInputRow);
  metadata.antenna1 = contiguousms._antenna1Column(_currentInputRow);
  metadata.antenna2 = contiguousms._antenna2Column(_currentInputRow);
}

void ContiguousMSReader::ReadData(std::complex<float>* buffer) {
  ContiguousMS& contiguousms = static_cast<ContiguousMS&>(*ms_provider_);

  readData();
  readWeights();
  const size_t data_desc_id = contiguousms._dataDescIdColumn(_currentInputRow);
  const size_t start_channel = contiguousms.channel_ranges_[data_desc_id].start;
  const size_t end_channel = contiguousms.channel_ranges_[data_desc_id].end;
  schaapcommon::reordering::ExtractData(
      buffer, start_channel, end_channel,
      contiguousms.input_polarizations_[data_desc_id],
      contiguousms._dataArray.data(), contiguousms._outputPolarization);
}

void ContiguousMSReader::ReadModel(std::complex<float>* buffer) {
  ContiguousMS& contiguousms = static_cast<ContiguousMS&>(*ms_provider_);

  if (!contiguousms._isModelColumnPrepared) contiguousms.prepareModelColumn();

  readModel();
  readWeights();
  const size_t data_desc_id = contiguousms._dataDescIdColumn(_currentInputRow);
  const size_t start_channel = contiguousms.channel_ranges_[data_desc_id].start;
  const size_t end_channel = contiguousms.channel_ranges_[data_desc_id].end;
  schaapcommon::reordering::ExtractData(
      buffer, start_channel, end_channel,
      contiguousms.input_polarizations_[data_desc_id],
      contiguousms._modelArray.data(), contiguousms._outputPolarization);
}

void ContiguousMSReader::ReadWeights(float* buffer) {
  const ContiguousMS& contiguousms =
      static_cast<const ContiguousMS&>(*ms_provider_);

  readData();
  readWeights();
  const size_t data_desc_id = contiguousms._dataDescIdColumn(_currentInputRow);
  const size_t start_channel = contiguousms.channel_ranges_[data_desc_id].start;
  const size_t end_channel = contiguousms.channel_ranges_[data_desc_id].end;
  schaapcommon::reordering::ExtractWeights(
      buffer, start_channel, end_channel,
      contiguousms.input_polarizations_[data_desc_id],
      contiguousms._dataArray.data(), contiguousms._weightSpectrumArray.data(),
      contiguousms._flagArray.data(), contiguousms._outputPolarization);
}

void ContiguousMSReader::WriteImagingWeights(const float* buffer) {
  ContiguousMS& contiguousms = static_cast<ContiguousMS&>(*ms_provider_);

  if (_imagingWeightsColumn == nullptr) {
    _imagingWeightsColumn.reset(new casacore::ArrayColumn<float>(
        InitializeImagingWeightColumn(*(contiguousms._ms))));
  }
  const size_t data_desc_id = contiguousms._dataDescIdColumn(_currentInputRow);
  const size_t start_channel = contiguousms.channel_ranges_[data_desc_id].start;
  const size_t end_channel = contiguousms.channel_ranges_[data_desc_id].end;

  _imagingWeightsColumn->get(_currentInputRow,
                             contiguousms._imagingWeightSpectrumArray);
  schaapcommon::reordering::StoreWeights(
      contiguousms._imagingWeightSpectrumArray.data(), start_channel,
      end_channel, contiguousms.input_polarizations_[data_desc_id], buffer,
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
      ExpandScalarWeights(contiguousms._weightScalarArray,
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
