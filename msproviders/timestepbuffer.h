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

  std::unique_ptr<MSReader> MakeReader() final;

  void NextOutputRow() final { ms_provider_->NextOutputRow(); }

  void ResetWritePosition() final { ms_provider_->ResetWritePosition(); }

  virtual void WriteModel(const std::complex<float>* buffer,
                          bool addToMS) final {
    ms_provider_->WriteModel(buffer, addToMS);
  }

  void ReopenRW() final { ms_provider_->ReopenRW(); }

  double StartTime() final { return ms_provider_->StartTime(); }

  aocommon::PolarizationEnum Polarization() final {
    return ms_provider_->Polarization();
  }

  bool IsRegular() const final { return ms_provider_->IsRegular(); }

  bool HasFrequencyBda() const override {
    return ms_provider_->HasFrequencyBda();
  }

  size_t NMaxChannels() final { return ms_provider_->NMaxChannels(); }

  size_t NAntennas() final { return ms_provider_->NAntennas(); }

  size_t NPolarizations() final { return ms_provider_->NPolarizations(); }

  size_t NRows() final { return ms_provider_->NRows(); }

  double Interval() final { return ms_provider_->Interval(); }

  ObservationInfo GetObservationInfo() final {
    return ms_provider_->GetObservationInfo();
  }

  std::vector<std::string> GetAntennaNames() final {
    return ms_provider_->GetAntennaNames();
  }

  const aocommon::MultiBandData& SelectedBands() final {
    return ms_provider_->SelectedBands();
  }

  std::optional<SynchronizedMS> MsIfAvailable() final {
    return ms_provider_->MsIfAvailable();
  }

  std::string PartDescription() const final {
    return ms_provider_->PartDescription();
  }

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
