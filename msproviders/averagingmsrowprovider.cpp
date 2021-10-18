#include "averagingmsrowprovider.h"

#include "../structures/multibanddata.h"

#include "../io/logger.h"

#include <casacore/tables/Tables/ArrayColumn.h>

AveragingMSRowProvider::AveragingMSRowProvider(
    double nWavelengthsAveraging, const string& msPath,
    const MSSelection& selection,
    const std::map<size_t, size_t>& selectedDataDescIds, size_t fieldId,
    const string& dataColumnName, bool requireModel)
    : MSRowProvider(msPath, selection, selectedDataDescIds, dataColumnName,
                    requireModel),
      _fieldId(fieldId) {
  casacore::MSAntenna antennaTable(_ms.antenna());
  _nAntennae = antennaTable.nrow();

  casacore::ArrayColumn<double> positionColumn(
      antennaTable,
      casacore::MSAntenna::columnName(casacore::MSAntennaEnums::POSITION));
  std::vector<Pos> positions(_nAntennae);

  casacore::Array<double> posArr(casacore::IPosition(1, 3));
  for (size_t i = 0; i != _nAntennae; ++i) {
    positionColumn.get(i, posArr);
    positions[i] = Pos(posArr.data()[0], posArr.data()[1], posArr.data()[2]);
  }

  // dataDescId x ant x ant
  _nElements = selectedDataDescIds.size() * _nAntennae * _nAntennae;
  _averagingFactors.assign(_nElements, 0.0);
  _buffers.resize(_nElements);
  MultiBandData bands(_ms.spectralWindow(), _ms.dataDescription());

  double dt = (EndTime() - StartTime()) / (EndTimestep() - StartTimestep());
  Logger::Debug << "Assuming integration time of " << dt * (24.0 * 60.0 * 60.0)
                << " seconds.\n";

  size_t element = 0;
  size_t minAvgFactor = std::numeric_limits<size_t>::max(), maxAvgFactor = 0;
  for (size_t a1 = 0; a1 != _nAntennae; ++a1) {
    Pos pos1 = positions[a1];
    for (size_t a2 = 0; a2 != _nAntennae; ++a2) {
      Pos pos2 = positions[a2];
      double dx = std::get<0>(pos1) - std::get<0>(pos2);
      double dy = std::get<1>(pos1) - std::get<1>(pos2);
      double dz = std::get<2>(pos1) - std::get<2>(pos2);
      double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
      for (std::map<size_t, size_t>::const_iterator spwIter =
               selectedDataDescIds.begin();
           spwIter != selectedDataDescIds.end(); ++spwIter) {
        BandData band = bands[spwIter->first];
        double lambda = band.SmallestWavelength();
        double nWavelengthsPerIntegration = 2.0 * M_PI * dist / lambda * dt;
        _averagingFactors[element] = std::max<size_t>(
            size_t(floor(nWavelengthsAveraging / nWavelengthsPerIntegration)),
            1);
        if (a1 != a2) {
          minAvgFactor =
              std::min<size_t>(minAvgFactor, _averagingFactors[element]);
          maxAvgFactor =
              std::max<size_t>(maxAvgFactor, _averagingFactors[element]);
        }
        ++element;
      }
    }
  }
  Logger::Info << "Averaging factor for longest baseline: " << minAvgFactor
               << " x . For the shortest: " << maxAvgFactor << " x \n";

  _spwIndexToDataDescId.resize(selectedDataDescIds.size());
  for (std::map<size_t, size_t>::const_iterator spwIter =
           selectedDataDescIds.begin();
       spwIter != selectedDataDescIds.end(); ++spwIter) {
    _spwIndexToDataDescId[spwIter->second] = spwIter->first;
  }

  _averageFactorSum = 0.0;
  _rowCount = 0;
  _averagedRowCount = 0;

  _averagedDataDescId = _currentDataDescId;
  _flushPosition = 0;

  if (!MSRowProvider::AtEnd()) {
    bool timestepAvailable = processCurrentTimestep();
    if (!timestepAvailable) NextRow();
  }
}

bool AveragingMSRowProvider::processCurrentTimestep() {
  size_t a1 = _antenna1Column(_currentRow);
  size_t a2 = _antenna2Column(_currentRow);
  _averagedDataDescId = _currentDataDescId;
  _averagedAntenna1Index = a1;
  _averagedAntenna2Index = a2;

  size_t spwCount = selectedDataDescIds().size();
  size_t elementIndex = spwCount * (a2 + a1 * _nAntennae) +
                        selectedDataDescIds().find(_averagedDataDescId)->second;
  size_t avgFactor = _averagingFactors[elementIndex];
  _averageFactorSum += avgFactor;
  ++_rowCount;

  _dataColumn.get(_currentRow, _currentData, true);
  _flagColumn.get(_currentRow, _currentFlags, true);
  getCurrentWeights(_currentWeights, _currentData.shape());
  if (requireModel()) _modelColumn->get(_currentRow, _currentModel, true);

  if (avgFactor == 1)
    return true;
  else {
    size_t bufferSize = _currentData.shape()[0] * _currentData.shape()[1];
    AveragingBuffer& buffer = _buffers[elementIndex];
    if (!buffer.IsInitialized()) buffer.Initialize(bufferSize, requireModel());

    if (requireModel())
      buffer.AddDataAndModel(bufferSize, _currentData.data(),
                             _currentModel.data(), _currentFlags.data(),
                             _currentWeights.data(), _currentUVWArray.data(),
                             _currentTime);
    else
      buffer.AddData(bufferSize, _currentData.data(), _currentFlags.data(),
                     _currentWeights.data(), _currentUVWArray.data(),
                     _currentTime);

    bool foundFullBuffer = (buffer.AveragedDataCount() == avgFactor);
    if (foundFullBuffer) {
      if (requireModel())
        buffer.Get(bufferSize, _currentData.data(), _currentModel.data(),
                   _currentFlags.data(), _currentWeights.data(),
                   _currentUVWArray.data(), _currentTime);
      else
        buffer.Get(bufferSize, _currentData.data(), _currentFlags.data(),
                   _currentWeights.data(), _currentUVWArray.data(),
                   _currentTime);
      buffer.Reset(bufferSize);
    }
    return foundFullBuffer;
  }
}

void AveragingMSRowProvider::NextRow() {
  ++_averagedRowCount;
  if (!MSRowProvider::AtEnd()) {
    bool foundFullBuffer = false;
    while (!foundFullBuffer) {
      MSRowProvider::NextRow();
      if (MSRowProvider::AtEnd()) break;

      foundFullBuffer = processCurrentTimestep();
    }
    // TODO I think this should now say "if(foundFullBuffer) return; "
  }

  if (MSRowProvider::AtEnd()) {
    // There might be residual data in the buffers which have to be read out
    AveragingBuffer* buffer = 0;
    do {
      buffer = &_buffers[_flushPosition];
      ++_flushPosition;
    } while (buffer->AveragedDataCount() == 0 && _flushPosition < _nElements);

    if (buffer != 0 && buffer->AveragedDataCount() != 0) {
      size_t bufferSize = _currentData.shape()[0] * _currentData.shape()[1];
      if (requireModel())
        buffer->Get(bufferSize, _currentData.data(), _currentModel.data(),
                    _currentFlags.data(), _currentWeights.data(),
                    _currentUVWArray.data(), _currentTime);
      else
        buffer->Get(bufferSize, _currentData.data(), _currentFlags.data(),
                    _currentWeights.data(), _currentUVWArray.data(),
                    _currentTime);

      size_t elementIndex = _flushPosition - 1;
      size_t spwCount = selectedDataDescIds().size();
      size_t spwIndex = elementIndex % spwCount;
      size_t subIndex = elementIndex / spwCount;
      _averagedAntenna1Index = subIndex / _nAntennae;
      _averagedAntenna2Index = subIndex % _nAntennae;
      _averagedDataDescId = _spwIndexToDataDescId[spwIndex];
    }
  }
}

void AveragingMSRowProvider::ReadData(MSRowProvider::DataArray& data,
                                      MSRowProvider::FlagArray& flags,
                                      MSRowProvider::WeightArray& weights,
                                      double& u, double& v, double& w,
                                      uint32_t& dataDescId, uint32_t& antenna1,
                                      uint32_t& antenna2, uint32_t& fieldId,
                                      double& time) {
  const size_t bufferSize = data.shape()[0] * data.shape()[1];
  std::copy_n(_currentData.data(), bufferSize, data.data());
  std::copy_n(_currentFlags.data(), bufferSize, flags.data());
  std::copy_n(_currentWeights.data(), bufferSize, weights.data());
  u = _currentUVWArray.data()[0];
  v = _currentUVWArray.data()[1];
  w = _currentUVWArray.data()[2];
  dataDescId = _averagedDataDescId;
  antenna1 = _averagedAntenna1Index;
  antenna2 = _averagedAntenna2Index;
  time = _currentTime;
  fieldId = _fieldId;  // TODO multi-fields are not supported
}

void AveragingMSRowProvider::ReadModel(MSRowProvider::DataArray& model) {
  size_t bufferSize = _currentModel.shape()[0] * _currentModel.shape()[1];
  std::copy_n(_currentModel.data(), bufferSize, model.data());
}

void AveragingMSRowProvider::OutputStatistics() const {
  // Logger::Info << "Average selected integration-time averaging factor after
  // row selection: " << double(_averageFactorSum)/_rowCount << '\n';
  Logger::Info << "Baseline averaging reduced the number of rows to "
               << 0.1 * round(_averagedRowCount * 1000.0 / _rowCount) << "%.\n";
}

void AveragingMSRowProvider::AveragingBuffer::Initialize(size_t bufferSize,
                                                         bool includeModel) {
  _modelData.resize(includeModel ? bufferSize : 0);
  Reset(bufferSize);
}

void AveragingMSRowProvider::AveragingBuffer::AddData(
    size_t n, const std::complex<float>* data, const bool* flags,
    const float* weights, const double* uvw, double time) {
  double weightSum = 0.0;
  for (size_t i = 0; i != n; ++i) {
    if (!flags[i] && std::isfinite(data[i].real()) &&
        std::isfinite(data[i].imag())) {
      _data[i] += data[i] * weights[i];
      _weights[i] += weights[i];
      weightSum += weights[i];
    }
  }
  // This division (and in the next function) seems not strictly required and
  // could be removed, because all rows will have the same n, and it just scales
  // both weightSum and _summedWeight down by n, and they're divided by each
  // other later on, so it cancels out. TODO verify, fix and test this.
  weightSum /= n;
  _uvw[0] += uvw[0] * weightSum;
  _uvw[1] += uvw[1] * weightSum;
  _uvw[2] += uvw[2] * weightSum;
  _time += time * weightSum;
  _unweightedTime = time;
  _summedWeight += weightSum;
  ++_averagedDataCount;
}

void AveragingMSRowProvider::AveragingBuffer::AddDataAndModel(
    size_t n, const std::complex<float>* data,
    const std::complex<float>* modelData, const bool* flags,
    const float* weights, const double* uvw, double time) {
  double weightSum = 0.0;
  for (size_t i = 0; i != n; ++i) {
    if (!flags[i] && std::isfinite(data[i].real()) &&
        std::isfinite(data[i].imag()) && std::isfinite(modelData[i].real()) &&
        std::isfinite(modelData[i].imag())) {
      _data[i] += data[i] * weights[i];
      _modelData[i] += modelData[i] * weights[i];
      _weights[i] += weights[i];
      weightSum += weights[i];
    }
  }
  weightSum /= n;
  _uvw[0] += uvw[0] * weightSum;
  _uvw[1] += uvw[1] * weightSum;
  _uvw[2] += uvw[2] * weightSum;
  _time += time * weightSum;
  _unweightedTime = time;
  _summedWeight += weightSum;
  ++_averagedDataCount;
}

void AveragingMSRowProvider::AveragingBuffer::Get(size_t n,
                                                  std::complex<float>* data,
                                                  bool* flags, float* weights,
                                                  double* uvw, double& time) {
  for (size_t i = 0; i != n; ++i) {
    flags[i] = (weights[i] == 0.0);
    data[i] = flags[i] ? 0.0 : _data[i] / _weights[i];
    weights[i] = _weights[i];
  }
  if (_summedWeight == 0.0) {
    for (size_t i = 0; i != 3; ++i) uvw[i] = 0.0;
    time = _unweightedTime;
  } else {
    for (size_t i = 0; i != 3; ++i) uvw[i] = _uvw[i] / _summedWeight;
    time = _time / _summedWeight;
  }
}

void AveragingMSRowProvider::AveragingBuffer::Get(
    size_t n, std::complex<float>* data, std::complex<float>* modelData,
    bool* flags, float* weights, double* uvw, double& time) {
  for (size_t i = 0; i != n; ++i) {
    flags[i] = (weights[i] == 0.0);
    data[i] = flags[i] ? 0.0 : _data[i] / _weights[i];
    modelData[i] = flags[i] ? 0.0 : _modelData[i] / weights[i];
    weights[i] = _weights[i];
  }
  if (_summedWeight == 0.0) {
    for (size_t i = 0; i != 3; ++i) uvw[i] = 0.0;
    time = _unweightedTime;
  } else {
    for (size_t i = 0; i != 3; ++i) uvw[i] = _uvw[i] / _summedWeight;
    time = _time / _summedWeight;
  }
}

void AveragingMSRowProvider::AveragingBuffer::Reset(size_t n) {
  _data.assign(n, 0.0);
  _weights.assign(n, 0.0);
  if (!_modelData.empty()) _modelData.assign(n, 0.0);
  _uvw[0] = 0.0;
  _uvw[1] = 0.0;
  _uvw[2] = 0.0;
  _time = 0.0;
  _unweightedTime = 0.0;
  _averagedDataCount = 0;
  _summedWeight = 0.0;
}
