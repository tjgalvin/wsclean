#ifndef MSPROVIDERS_MSREADERS_TIMESTEPBUFFERREADER_H
#define MSPROVIDERS_MSREADERS_TIMESTEPBUFFERREADER_H

#include "msreader.h"
#include "../timestepbuffer.h"

class TimestepBufferReader final : public MSReader {
 public:
  TimestepBufferReader(TimestepBuffer* timestepBuffer);

  virtual ~TimestepBufferReader(){};

  size_t RowId() const final override { return _buffer[_bufferPosition].rowId; }

  bool CurrentRowAvailable() final override;

  void NextInputRow() final override;

  void ReadMeta(double& u, double& v, double& w,
                size_t& dataDescId) final override;

  void ReadMeta(MSProvider::MetaData& metaData) final override;

  void ReadData(std::complex<float>* buffer) final override;

  void ReadModel(std::complex<float>* buffer) final override;

  void ReadWeights(std::complex<float>* buffer) final override;

  void ReadWeights(float* buffer) final override;

  void WriteImagingWeights(const float* buffer) final override;

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

 private:
  void readTimeblock();

  std::unique_ptr<MSReader> _msReader;

  size_t _bufferPosition;
  std::vector<TimestepBuffer::RowData> _buffer;
};

#endif