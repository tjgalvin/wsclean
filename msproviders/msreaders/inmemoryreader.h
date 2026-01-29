#ifndef MSPROVIDERS_MSREADERS_INMEMORYREADER_H_
#define MSPROVIDERS_MSREADERS_INMEMORYREADER_H_

#include <span>

#include "msreader.h"

#include "../../structures/inmemorypart.h"

namespace wsclean {

class InMemoryReader final : public MSReader {
 public:
  InMemoryReader(
      MSProvider* ms_provider, std::span<const InMemoryPartRow> rows,
      std::vector<MSProvider::MetaData>::const_iterator meta_data_iterator)
      : MSReader(ms_provider),
        rows_(rows),
        data_iterator_(rows.begin()),
        meta_data_iterator_(meta_data_iterator){};

  size_t RowId() const final { return data_iterator_ - rows_.begin(); }

  bool CurrentRowAvailable() final { return data_iterator_ != rows_.end(); }

  void NextInputRow() final {
    ++data_iterator_;
    ++meta_data_iterator_;
  }

  void ReadMeta(MSProvider::MetaData& metadata) final {
    metadata = *meta_data_iterator_;
  }

  void ReadData(std::complex<float>* buffer) final {
    std::copy(data_iterator_->data.begin(), data_iterator_->data.end(), buffer);
  }

  void ReadModel(std::complex<float>* buffer) final {
    std::copy(data_iterator_->model_data.begin(),
              data_iterator_->model_data.end(), buffer);
  }

  void ReadWeights(float* buffer) final {
    std::copy(data_iterator_->weights.begin(), data_iterator_->weights.end(),
              buffer);
  }

  void WriteImagingWeights(const float* buffer) final {
    // In memory reader does not hold imaging weights, so do nothing.
  }

 private:
  std::span<const InMemoryPartRow> rows_;
  std::span<const InMemoryPartRow>::iterator data_iterator_;
  std::vector<MSProvider::MetaData>::const_iterator meta_data_iterator_;
};

}  // namespace wsclean

#endif
