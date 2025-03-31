#ifndef MSPROVIDERS_TIMESTEP_BUFFER_H
#define MSPROVIDERS_TIMESTEP_BUFFER_H

#include "msprovider.h"

#include <aocommon/uvector.h>

namespace wsclean {

class TimestepBufferReader;

/**
 * This class wraps any MSProvider to make it read whole blocks of rows
 * at once that correspond to the same timestep.
 *
 * This is used in IDGMSGridder to be able to get the UVWs for calculating the
 * a-terms.
 */
class TimestepBuffer final : public MSProvider {
  friend class TimestepBufferReader;

 public:
  TimestepBuffer(MSProvider* ms_provider, bool read_model)
      : ms_provider_(ms_provider), read_model_(read_model) {
    ms_provider_->ResetWritePosition();
  }

  virtual ~TimestepBuffer(){};

  SynchronizedMS MS() override { return ms_provider_->MS(); }

  std::unique_ptr<MSReader> MakeReader() final override;

  const std::string& DataColumnName() override {
    return ms_provider_->DataColumnName();
  }

  void NextOutputRow() override { ms_provider_->NextOutputRow(); }

  void ResetWritePosition() override { ms_provider_->ResetWritePosition(); }

  virtual void WriteModel(const std::complex<float>* buffer,
                          bool addToMS) override {
    ms_provider_->WriteModel(buffer, addToMS);
  }

  void ReopenRW() override { ms_provider_->ReopenRW(); }

  double StartTime() override { return ms_provider_->StartTime(); }

  size_t DataDescId() override { return ms_provider_->DataDescId(); }

  void MakeIdToMSRowMapping(std::vector<size_t>& idToMSRow) override {
    ms_provider_->MakeIdToMSRowMapping(idToMSRow);
  }

  aocommon::PolarizationEnum Polarization() override {
    return ms_provider_->Polarization();
  }

  size_t NChannels() override { return ms_provider_->NChannels(); }

  size_t NAntennas() override { return ms_provider_->NAntennas(); }

  size_t NPolarizations() override { return ms_provider_->NPolarizations(); }

  const aocommon::BandData& Band() override { return ms_provider_->Band(); }

 private:
  struct RowData {
    std::vector<std::complex<float>> data;
    std::vector<std::complex<float>> model;
    std::vector<float> weights;
    MSProvider::MetaData metadata;
    size_t row_id;
  };

  MSProvider* ms_provider_;

  bool read_model_;
};

}  // namespace wsclean

#endif
