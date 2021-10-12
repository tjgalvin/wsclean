#include "timestepbufferreader.h"

TimestepBufferReader::TimestepBufferReader(TimestepBuffer* timestepBuffer)
    : MSReader(timestepBuffer),
      _msReader(timestepBuffer->_msProvider->MakeReader()),
      _bufferPosition(0) {
  readTimeblock();
}

bool TimestepBufferReader::CurrentRowAvailable() {
  return !_buffer.empty() || _msReader->CurrentRowAvailable();
}

void TimestepBufferReader::NextInputRow() {
  ++_bufferPosition;
  if (_bufferPosition == _buffer.size()) {
    readTimeblock();
  }
}

void TimestepBufferReader::ReadMeta(double& u, double& v, double& w,
                                    size_t& dataDescId) {
  MSProvider::MetaData& m = _buffer[_bufferPosition].metaData;
  u = m.uInM;
  v = m.vInM;
  w = m.wInM;
  dataDescId = m.dataDescId;
}

void TimestepBufferReader::ReadMeta(MSProvider::MetaData& metaData) {
  metaData = _buffer[_bufferPosition].metaData;
}

void TimestepBufferReader::ReadData(std::complex<float>* buffer) {
  std::copy(_buffer[_bufferPosition].data.begin(),
            _buffer[_bufferPosition].data.end(), buffer);
}

void TimestepBufferReader::ReadModel(std::complex<float>* buffer) {
  std::copy(_buffer[_bufferPosition].model.begin(),
            _buffer[_bufferPosition].model.end(), buffer);
}

void TimestepBufferReader::ReadWeights(float* buffer) {
  std::copy(_buffer[_bufferPosition].weights.begin(),
            _buffer[_bufferPosition].weights.end(), buffer);
}

void TimestepBufferReader::ReadWeights(std::complex<float>* buffer) {
  std::copy(_buffer[_bufferPosition].weights.begin(),
            _buffer[_bufferPosition].weights.end(), buffer);
}

void TimestepBufferReader::WriteImagingWeights(const float* buffer) {
  _msReader->WriteImagingWeights(buffer);
}

void TimestepBufferReader::readTimeblock() {
  // Beware that the _msProvider data member is a TimestepBuffer,
  // which in turn has its own _msProvider
  TimestepBuffer& tstepbuffer = static_cast<TimestepBuffer&>(*_msProvider);

  _bufferPosition = 0;
  _buffer.clear();
  MSProvider::MetaData metaData;
  size_t dataSize = tstepbuffer._msProvider->NPolarizations() *
                    tstepbuffer._msProvider->NChannels();

  if (_msReader->CurrentRowAvailable()) {
    _msReader->ReadMeta(metaData);
    double blockTime = metaData.time, curTime = blockTime;
    size_t writePos = 0;
    do {
      if (_buffer.size() <= writePos) {
        _buffer.emplace_back();
        TimestepBuffer::RowData& newRow = _buffer.back();
        newRow.data.resize(dataSize);
        if (tstepbuffer._readModel) newRow.model.resize(dataSize);
        newRow.weights.resize(dataSize);
      }
      TimestepBuffer::RowData& row = _buffer[writePos];
      row.metaData = metaData;
      _msReader->ReadData(row.data.data());
      if (tstepbuffer._readModel) _msReader->ReadModel(row.model.data());
      _msReader->ReadWeights(row.weights.data());
      row.rowId = _msReader->RowId();

      _msReader->NextInputRow();
      ++writePos;
      if (_msReader->CurrentRowAvailable()) {
        _msReader->ReadMeta(metaData);
        curTime = metaData.time;
      }
    } while (_msReader->CurrentRowAvailable() && blockTime == curTime);
    _buffer.resize(writePos);
  }
}
