#ifndef MSPROVIDERS_REORDERED_MS_PROVIDER_
#define MSPROVIDERS_REORDERED_MS_PROVIDER_

#include "msprovider.h"

#include "../structures/msselection.h"
#include "../system/mappedfile.h"

#include <schaapcommon/reordering/reorderedhandle.h>
#include <schaapcommon/reordering/reorderedfilewriter.h>
#include <schaapcommon/reordering/reordering.h>

#include <aocommon/io/serialstreamfwd.h>
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

class ReorderedMsProvider final : public MSProvider {
  friend class ReorderedMsReader;

 public:
  class ReorderedHandle;

  ReorderedMsProvider(const ReorderedHandle& handle, size_t part_index,
                      aocommon::PolarizationEnum polarization,
                      size_t data_desc_id);

  virtual ~ReorderedMsProvider();

  ReorderedMsProvider(const ReorderedMsProvider&) = delete;
  ReorderedMsProvider& operator=(const ReorderedMsProvider&) = delete;

  std::unique_ptr<MSReader> MakeReader() override;

  SynchronizedMS MS() override {
    return SynchronizedMS(handle_.data_->ms_path_.data());
  }

  const std::string& DataColumnName() override {
    return handle_.data_->data_column_name_;
  }

  void NextOutputRow() override;

  void ResetWritePosition() override { current_output_row_ = 0; };

  void WriteModel(const std::complex<float>* buffer, bool add_to_MS) override;

  void ReopenRW() override {}

  double StartTime() override { return meta_header_.start_time; }

  void MakeIdToMSRowMapping(std::vector<size_t>& id_to_MS_row) override;

  aocommon::PolarizationEnum Polarization() override { return polarization_; }

  size_t NChannels() override { return part_header_.channel_count; }
  size_t NPolarizations() override { return polarization_count_in_file_; }
  size_t NAntennas() override { return handle_.data_->n_antennas_; }

  size_t DataDescId() override { return part_header_.data_desc_id; }

  const aocommon::BandData& Band() override {
    return handle_.data_->bands_[data_desc_id_];
  }

  class ReorderedHandle {
    // ReorderedMsReader is a friend of ReorderedHandle
    // in order to access the data_ member.
    friend class ReorderedMsReader;
    friend class ReorderedMsProvider;

   public:
    ReorderedHandle() = default;

    ReorderedHandle(
        const std::string& ms_path, const string& data_column_name,
        const std::string& model_column_name,
        schaapcommon::reordering::StorageManagerType model_storage_manager,
        const std::string& temporary_directory,
        const std::vector<schaapcommon::reordering::ChannelRange>& channels,
        bool initial_model_required, bool model_update_required,
        const std::set<aocommon::PolarizationEnum>& polarizations,
        const schaapcommon::reordering::MSSelection& selection,
        const aocommon::MultiBandData& bands, size_t n_antennas,
        bool keep_temporary_files,
        std::function<
            void(schaapcommon::reordering::ReorderedHandleData& handle)>
            cleanup_callback)
        : data_(std::make_shared<schaapcommon::reordering::ReorderedHandleData>(
              ms_path, data_column_name, model_column_name,
              model_storage_manager, temporary_directory, channels,
              initial_model_required, model_update_required, polarizations,
              selection, bands, n_antennas, keep_temporary_files,
              std::move(cleanup_callback))) {}

    void Serialize(aocommon::SerialOStream& stream) const;
    void Unserialize(aocommon::SerialIStream& stream);

   private:
    std::shared_ptr<schaapcommon::reordering::ReorderedHandleData> data_;
  };

  static void StoreReorderedInMS(
      const schaapcommon::reordering::ReorderedHandleData& handle);

 private:
  const ReorderedHandle handle_;
  const size_t part_index_;
  const size_t data_desc_id_;
  MappedFile model_file_;
  size_t current_output_row_;
  std::unique_ptr<std::ofstream> model_data_file_;
  const aocommon::PolarizationEnum polarization_;
  size_t polarization_count_in_file_;

  schaapcommon::reordering::MetaHeader meta_header_;
  schaapcommon::reordering::PartHeader part_header_;
};

ReorderedMsProvider::ReorderedHandle ReorderMS(
    const std::string& ms_path,
    const std::vector<schaapcommon::reordering::ChannelRange>& channels,
    const schaapcommon::reordering::MSSelection& selection,
    const std::string& data_column_name, const std::string& model_column_name,
    schaapcommon::reordering::StorageManagerType model_storage_manager,
    bool include_model, bool initial_model_required,
    const class Settings& settings);

}  // namespace wsclean

#endif
