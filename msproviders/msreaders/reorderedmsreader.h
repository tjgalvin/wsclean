#ifndef MSPROVIDERS_MSREADERS_REORDEREDMSREADER_H_
#define MSPROVIDERS_MSREADERS_REORDEREDMSREADER_H_

#include "msreader.h"

#include <aocommon/uvector.h>

#include <cstdint>
#include <fstream>

namespace wsclean {

class ReorderedMsProvider;

class ReorderedMsReader final : public MSReader {
 public:
  ReorderedMsReader(ReorderedMsProvider* reordered_ms);
  virtual ~ReorderedMsReader(){};

  size_t RowId() const override { return current_input_row_; }

  bool CurrentRowAvailable() override;

  void NextInputRow() override;

  void ReadMeta(MSProvider::MetaData& metadata) override;

  void ReadData(std::complex<float>* buffer) override;

  void ReadModel(std::complex<float>* buffer) override;

  void ReadWeights(float* buffer) override;

  void WriteImagingWeights(const float* buffer) override;

 private:
  uint64_t NVisibilitiesPerRow() const;

  size_t current_input_row_;

  // Chunkoffset counts the amount of data rows we are ahead or behind.
  // Positive values mean we are lagging, whereas negative values mean we are
  //  ahead of the current time step.
  long read_ptr_row_offset_, meta_ptr_row_offset_, weight_ptr_row_offset_;

  std::ifstream meta_file_, weight_file_, data_file_;

  aocommon::UVector<float> imaging_weight_buffer_;
  std::unique_ptr<std::fstream> imaging_weights_file_;
};

}  // namespace wsclean

#endif
