#ifndef MSPROVIDERS_REORDERED_HANDLE_
#define MSPROVIDERS_REORDERED_HANDLE_

#include "../structures/msselection.h"

#include <functional>

#include <aocommon/io/serialstreamfwd.h>
#include <aocommon/polarization.h>

using schaapcommon::reordering::MSSelection;

namespace wsclean {
namespace reordering {

struct ChannelRange {
  size_t data_desc_id;
  size_t start, end;
  bool operator<(const ChannelRange& rhs) const {
    if (data_desc_id < rhs.data_desc_id) return true;
    if (data_desc_id > rhs.data_desc_id) return false;
    if (start < rhs.start) return true;
    if (start > rhs.start) return false;
    return end < rhs.end;
  }
};

struct ReorderedHandleData {
  ReorderedHandleData() = default;

  ReorderedHandleData(
      const std::string& ms_path, const string& data_column_name,
      const std::string& temporary_directory,
      const std::vector<reordering::ChannelRange>& channels,
      bool initial_model_required, bool model_update_required,
      const std::set<aocommon::PolarizationEnum>& polarizations,
      const MSSelection& selection, const aocommon::MultiBandData& bands,
      size_t n_antennas, bool keep_temporary_files,
      std::function<void(ReorderedHandleData& handle)> cleanup_callback)
      : ms_path_(ms_path),
        data_column_name_(data_column_name),
        temporary_directory_(temporary_directory),
        channels_(channels),
        initial_model_required_(initial_model_required),
        model_update_required_(model_update_required),
        polarizations_(polarizations),
        selection_(selection),
        bands_(bands),
        n_antennas_(n_antennas),
        keep_temporary_files_(keep_temporary_files),
        cleanup_callback_(std::move(cleanup_callback)) {}

  ~ReorderedHandleData();

  std::string ms_path_;
  std::string data_column_name_;
  std::string temporary_directory_;
  std::vector<ChannelRange> channels_;
  bool initial_model_required_;
  bool model_update_required_;
  std::set<aocommon::PolarizationEnum> polarizations_;
  MSSelection selection_;
  aocommon::MultiBandData bands_;
  size_t n_antennas_;
  bool is_copy_ = false;
  bool keep_temporary_files_;
  std::function<void(ReorderedHandleData& handle)> cleanup_callback_;

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);
};

}  // namespace reordering
}  // namespace wsclean

#endif
