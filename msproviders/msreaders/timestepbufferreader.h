#ifndef MSPROVIDERS_MSREADERS_TIMESTEPBUFFERREADER_H
#define MSPROVIDERS_MSREADERS_TIMESTEPBUFFERREADER_H

#include "msreader.h"
#include "../timestepbuffer.h"

namespace wsclean {

class TimestepBufferReader final : public MSReader {
 public:
  TimestepBufferReader(TimestepBuffer* timestep_buffer);

  virtual ~TimestepBufferReader(){};

  size_t RowId() const final override {
    return buffer_[buffer_position_].row_id;
  }

  bool CurrentRowAvailable() final override;

  void NextInputRow() final override;

  void ReadMeta(MSProvider::MetaData& metadata) final override;

  void ReadData(std::complex<float>* buffer) final override;

  void ReadModel(std::complex<float>* buffer) final override;

  void ReadWeights(float* buffer) final override;

  void WriteImagingWeights(const float* buffer) final override;

  /**
   * Returns an Array containing the uvws for baselines (antenna1, antenna2)
   * that have antenna1=0, sorted by antenna2.
   * @param uvws should have the correct size on input (nantenna * 3)
   */
  void GetUVWsForTimestep(aocommon::UVector<double>& uvws) {
    for (size_t i = 0; i != buffer_.size(); ++i) {
      if (buffer_[i].metadata.antenna1 == 0) {
        size_t index = buffer_[i].metadata.antenna2 * 3;
        if (index >= buffer_.size()) buffer_.resize(index + 3);
        uvws[index + 0] = buffer_[i].metadata.u_in_m;
        uvws[index + 1] = buffer_[i].metadata.v_in_m;
        uvws[index + 2] = buffer_[i].metadata.w_in_m;
      }
    }
    uvws[0] = 0.0;
    uvws[1] = 0.0;
    uvws[2] = 0.0;
  }

 private:
  void readTimeblock();

  std::unique_ptr<MSReader> ms_reader_;

  size_t buffer_position_;
  std::vector<TimestepBuffer::RowData> buffer_;
};

}  // namespace wsclean

#endif
