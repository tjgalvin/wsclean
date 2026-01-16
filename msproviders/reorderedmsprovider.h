#ifndef MSPROVIDERS_REORDERED_MS_PROVIDER_
#define MSPROVIDERS_REORDERED_MS_PROVIDER_

#include "msprovider.h"

#include "../structures/msselection.h"
#include "../system/mappedfile.h"

#include "reorderedhandle.h"
#include "msreaders/reorderedmsreader.h"

#include <schaapcommon/reordering/handledata.h>
#include <schaapcommon/reordering/filewriter.h>
#include <schaapcommon/reordering/reordering.h>

#include <aocommon/polarization.h>
#include <aocommon/uvector.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <functional>
#include <fstream>
#include <string>

namespace wsclean {

class ReorderedMsReader;
class Settings;

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
  bool IsRegular() const override {
    return meta_header_.data_desc_id.HasValue();
  }
  bool HasFrequencyBda() const override {
    return handle_.data_->has_frequency_bda_;
  }
  size_t NPolarizations() override { return polarization_count_in_file_; }
  size_t NAntennas() override { return handle_.data_->n_antennas_; }

  const aocommon::MultiBandData& SelectedBands() final {
    return handle_.data_->bands_per_part_[part_index_];
  }

  double Interval() final;

  ObservationInfo GetObservationInfo() final;

  std::vector<std::string> GetAntennaNames() final;

  static void StoreReorderedInMS(
      const schaapcommon::reordering::HandleData& handle);

  std::optional<SynchronizedMS> MsIfAvailable() final { return MS(); }

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
  std::optional<ReorderedMsReader> reader_;

  schaapcommon::reordering::MetaHeader meta_header_;
  schaapcommon::reordering::PartHeader part_header_;
};

ReorderedHandle ReorderMS(
    const std::string& ms_path,
    const std::vector<
        aocommon::VectorMap<schaapcommon::reordering::ChannelRange>>& channels,
    const schaapcommon::reordering::MSSelection& selection,
    const std::string& data_column_name, const std::string& model_column_name,
    schaapcommon::reordering::StorageManagerType model_storage_manager,
    bool include_model, bool initial_model_required, const Settings& settings);

}  // namespace wsclean

#endif
