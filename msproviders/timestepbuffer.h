#ifndef MSPROVIDERS_TIMESTEP_BUFFER_H
#define MSPROVIDERS_TIMESTEP_BUFFER_H

#include "msprovider.h"

#include <aocommon/uvector.h>

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
  TimestepBuffer(MSProvider* msProvider, bool readModel)
      : _msProvider(msProvider), _readModel(readModel) {
    _msProvider->ResetWritePosition();
  }

  virtual ~TimestepBuffer(){};

  SynchronizedMS MS() override { return _msProvider->MS(); }

  std::unique_ptr<MSReader> MakeReader() final override;

  const std::string& DataColumnName() override {
    return _msProvider->DataColumnName();
  }

  void NextOutputRow() override { _msProvider->NextOutputRow(); }

  void ResetWritePosition() override { _msProvider->ResetWritePosition(); }

  virtual void WriteModel(const std::complex<float>* buffer,
                          bool addToMS) override {
    _msProvider->WriteModel(buffer, addToMS);
  }

  void ReopenRW() override { _msProvider->ReopenRW(); }

  double StartTime() override { return _msProvider->StartTime(); }

  size_t DataDescId() override { return _msProvider->DataDescId(); }

  void MakeIdToMSRowMapping(std::vector<size_t>& idToMSRow) override {
    _msProvider->MakeIdToMSRowMapping(idToMSRow);
  }

  aocommon::PolarizationEnum Polarization() override {
    return _msProvider->Polarization();
  }

  size_t NChannels() override { return _msProvider->NChannels(); }

  size_t NAntennas() override { return _msProvider->NAntennas(); }

  size_t NPolarizations() override { return _msProvider->NPolarizations(); }

  const aocommon::BandData& Band() override { return _msProvider->Band(); }

 private:
  struct RowData {
    std::vector<std::complex<float>> data, model;
    std::vector<float> weights;
    MSProvider::MetaData metaData;
    size_t rowId;
  };

  MSProvider* _msProvider;

  bool _readModel;
};

#endif
