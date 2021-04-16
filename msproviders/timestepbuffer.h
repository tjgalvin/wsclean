#ifndef MSPROVIDERS_TIMESTEP_BUFFER_H
#define MSPROVIDERS_TIMESTEP_BUFFER_H

#include "msprovider.h"

/**
 * This class wraps any MSProvider to make it read whole blocks of rows
 * at once that correspond to the same timestep.
 *
 * This is used in IDGMSGridder to be able to get the UVWs for calculating the
 * a-terms.
 */
class TimestepBuffer final : public MSProvider {
 public:
  TimestepBuffer(MSProvider* msProvider, bool readModel)
      : _msProvider(msProvider), _bufferPosition(0), _readModel(readModel) {
    _msProvider->Reset();
    readTimeblock();
  }

  SynchronizedMS MS() override { return _msProvider->MS(); }

  const std::string& DataColumnName() override {
    return _msProvider->DataColumnName();
  }

  bool CurrentRowAvailable() override {
    return !_buffer.empty() || _msProvider->CurrentRowAvailable();
  }

  void NextRow() override {
    ++_bufferPosition;
    if (_bufferPosition == _buffer.size()) readTimeblock();
  }

  void Reset() override {
    _msProvider->Reset();
    readTimeblock();
  }

  size_t RowId() const override { return _buffer[_bufferPosition].rowId; }

  void ReadMeta(double& u, double& v, double& w, size_t& dataDescId) override {
    MSProvider::MetaData& m = _buffer[_bufferPosition].metaData;
    u = m.uInM;
    v = m.vInM;
    w = m.wInM;
    dataDescId = m.dataDescId;
  }

  void ReadMeta(MSProvider::MetaData& metaData) override {
    metaData = _buffer[_bufferPosition].metaData;
  }

  void ReadData(std::complex<float>* buffer) override {
    std::copy(_buffer[_bufferPosition].data.begin(),
              _buffer[_bufferPosition].data.end(), buffer);
  }

  void ReadModel(std::complex<float>* buffer) override {
    std::copy(_buffer[_bufferPosition].model.begin(),
              _buffer[_bufferPosition].model.end(), buffer);
  }

  virtual void WriteModel(size_t rowId, const std::complex<float>* buffer,
                          bool addToMS) override {
    _msProvider->WriteModel(rowId, buffer, addToMS);
  }

  virtual void WriteImagingWeights(size_t rowId, const float* buffer) override {
    _msProvider->WriteImagingWeights(rowId, buffer);
  }

  virtual void ReadWeights(float* buffer) override {
    std::copy(_buffer[_bufferPosition].weights.begin(),
              _buffer[_bufferPosition].weights.end(), buffer);
  }

  virtual void ReadWeights(std::complex<float>* buffer) override {
    std::copy(_buffer[_bufferPosition].weights.begin(),
              _buffer[_bufferPosition].weights.end(), buffer);
  }

  /**
   * Returns an Array containing the uvws for baselines (antenna1, antenna2)
   * that have antenna1=0, sorted by antenna2.
   * @param uvws should have the correct size on input (nantenna * 3)
   */
  void GetUVWsForTimestep(aocommon::UVector<double>& uvws) {
    for (size_t i = 0; i != _buffer.size(); ++i) {
      if (_buffer[i].metaData.antenna1 == 0) {
        size_t index = _buffer[i].metaData.antenna2 * 3;
        if (index >= _buffer.size()) _buffer.resize(index + 3);
        uvws[index + 0] = _buffer[i].metaData.uInM;
        uvws[index + 1] = _buffer[i].metaData.vInM;
        uvws[index + 2] = _buffer[i].metaData.wInM;
      }
    }
    uvws[0] = 0.0;
    uvws[1] = 0.0;
    uvws[2] = 0.0;
  }

  void ReopenRW() override { _msProvider->ReopenRW(); }

  double StartTime() override { return _msProvider->StartTime(); }

  void MakeIdToMSRowMapping(std::vector<size_t>& idToMSRow) override {
    _msProvider->MakeIdToMSRowMapping(idToMSRow);
  }

  aocommon::PolarizationEnum Polarization() override {
    return _msProvider->Polarization();
  }

  size_t NChannels() override { return _msProvider->NChannels(); }

  size_t NAntennas() override { return _msProvider->NAntennas(); }

  size_t NPolarizations() override { return _msProvider->NPolarizations(); }

 private:
  void readTimeblock() {
    _bufferPosition = 0;
    _buffer.clear();
    MSProvider::MetaData metaData;
    size_t dataSize = _msProvider->NPolarizations() * _msProvider->NChannels();

    if (_msProvider->CurrentRowAvailable()) {
      _msProvider->ReadMeta(metaData);
      double blockTime = metaData.time, curTime = blockTime;
      size_t writePos = 0;
      do {
        if (_buffer.size() <= writePos) {
          _buffer.emplace_back();
          RowData& newRow = _buffer.back();
          newRow.data.resize(dataSize);
          if (_readModel) newRow.model.resize(dataSize);
          newRow.weights.resize(dataSize);
        }
        RowData& row = _buffer[writePos];
        row.metaData = metaData;
        _msProvider->ReadData(row.data.data());
        if (_readModel) _msProvider->ReadModel(row.model.data());
        _msProvider->ReadWeights(row.weights.data());
        row.rowId = _msProvider->RowId();

        _msProvider->NextRow();
        ++writePos;
        if (_msProvider->CurrentRowAvailable()) {
          _msProvider->ReadMeta(metaData);
          curTime = metaData.time;
        }
      } while (_msProvider->CurrentRowAvailable() && blockTime == curTime);
      _buffer.resize(writePos);
    }
  }

  struct RowData {
    std::vector<std::complex<float>> data, model;
    std::vector<float> weights;
    MSProvider::MetaData metaData;
    size_t rowId;
  };

  MSProvider* _msProvider;

  size_t _bufferPosition;
  std::vector<RowData> _buffer;
  bool _readModel;
};

#endif
