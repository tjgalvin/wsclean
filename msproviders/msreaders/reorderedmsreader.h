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

  size_t CurrentDataDescId() const { return metadata_.data_desc_id; }

 private:
  uint64_t NVisibilitiesPerRow(size_t data_desc_id) const;
  void CacheMeta();

  size_t current_input_row_ = 0;
  size_t current_value_position_ = 0;
  MSProvider::MetaData metadata_;

  // Chunkoffset counts the amount of data rows we are ahead or behind.
  // Positive values mean we are lagging, whereas negative values mean we are
  //  ahead of the current time step.
  bool has_data_ = false;
  bool has_weights_ = false;

  std::ifstream meta_file_;
  std::ifstream weight_file_;
  std::ifstream data_file_;

  aocommon::UVector<float> imaging_weight_buffer_;
  aocommon::UVector<std::complex<float>> data_buffer_;
  aocommon::UVector<float> weight_buffer_;
  std::unique_ptr<std::fstream> imaging_weights_file_;
};

}  // namespace wsclean

#endif
