#ifndef MSPROVIDERS_REORDERED_MS_PROVIDER_
#define MSPROVIDERS_REORDERED_MS_PROVIDER_

#include "msprovider.h"

#include <complex>
#include <functional>
#include <fstream>
#include <optional>
#include <string>

#include <aocommon/multibanddata.h>
#include <aocommon/polarization.h>
#include <aocommon/uvector.h>

#include <schaapcommon/reordering/reordering.h>

#include "synchronizedms.h"

#include "../structures/observationinfo.h"

#include "../system/mappedfile.h"

#include "reordering/reorderedhandle.h"

#include "msreaders/reorderedmsreader.h"

namespace wsclean {

class ReorderedMsReader;

class ReorderedMsProvider final : public MSProvider {
  friend class ReorderedMsReader;

 public:
  ReorderedMsProvider(const ReorderedHandle& handle, size_t part_index,
                      aocommon::PolarizationEnum polarization);

  ~ReorderedMsProvider() final;

  ReorderedMsProvider(const ReorderedMsProvider&) = delete;
  ReorderedMsProvider& operator=(const ReorderedMsProvider&) = delete;

  std::unique_ptr<MSReader> MakeReader() override;

  void NextOutputRow() override;

  void ResetWritePosition() override;

  void WriteModel(const std::complex<float>* buffer, bool add_to_MS) override;

  void ResetModelColumn() override;

  void ReopenRW() override {}

  double StartTime() override { return meta_header_.start_time; }

  aocommon::PolarizationEnum Polarization() override { return polarization_; }

  std::string PartDescription() const final;

  size_t NRows() final { return meta_header_.selected_row_count; }

  size_t NMaxChannels() override { return part_header_.max_channel_count; }
  bool IsRegular() const override { return meta_header_.data_desc_id.HasValue(); }
  bool HasFrequencyBda() const override { return handle_.data_->has_frequency_bda_; }
  size_t NPolarizations() override { return polarization_count_in_file_; }
  size_t NAntennas() override { return handle_.data_->n_antennas_; }

  const aocommon::MultiBandData& SelectedBands() final {
    return handle_.data_->bands_per_part_[part_index_];
  }

  ObservationInfo GetObservationInfo() final;

  std::vector<std::string> GetAntennaNames() final;

  std::optional<SynchronizedMS> MsIfAvailable() final { return MS(); }

  bool IsThreadSafe() const final { return true; }

 private:
  SynchronizedMS MS() { return SynchronizedMS(handle_.data_->ms_path_.data()); }

  const ReorderedHandle handle_;
  const size_t part_index_;
  MappedFile model_file_;
  std::string model_file_path_;
  size_t current_output_row_ = 0;
  size_t current_output_position_ = 0;
  std::unique_ptr<std::ofstream> model_data_file_;
  const aocommon::PolarizationEnum polarization_;
  size_t polarization_count_in_file_;
  double interval_ = 0.0;
  std::optional<ReorderedMsReader> reader_;

  schaapcommon::reordering::MetaHeader meta_header_;
  schaapcommon::reordering::PartHeader part_header_;
};

}  // namespace wsclean

#endif
