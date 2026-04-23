#include "reorderedmsprovider.h"

#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <memory>
#include <vector>

#include <aocommon/logger.h>

#include <schaapcommon/reordering/filewriter.h>

#include <casacore/measures/Measures/MEpoch.h>

#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/ms/MeasurementSets/MSDataDescription.h>

#include "operations.h"

#include "msreaders/reorderedmsreader.h"

#include "rowproviders/directmsrowprovider.h"
#include "rowproviders/msrowprovider.h"
#include "rowproviders/noisemsrowprovider.h"

using aocommon::Logger;

using casacore::ArrayColumn;
using casacore::MSMainEnums;
using casacore::ScalarColumn;

using schaapcommon::reordering::ChannelRange;
using schaapcommon::reordering::GetMaxChannels;
using schaapcommon::reordering::GetMetaFilename;
using schaapcommon::reordering::GetPartPrefix;
using schaapcommon::reordering::HandleData;
using schaapcommon::reordering::MetaHeader;
using schaapcommon::reordering::PartHeader;

namespace wsclean {

ReorderedMsProvider::ReorderedMsProvider(const ReorderedHandle& handle, size_t part_index,
                                         aocommon::PolarizationEnum polarization)
    : handle_(handle),
      part_index_(part_index),
      current_output_row_(0),
      polarization_(polarization),
      polarization_count_in_file_(aocommon::Polarization::GetVisibilityCount(polarization_)) {
  const size_t meta_file_index = handle_.data_->metadata_indices_[part_index_];
  std::ifstream meta_file(
      GetMetaFilename(handle.data_->ms_path_, handle.data_->temporary_directory_, meta_file_index));
  if (!meta_file) {
    throw std::runtime_error("Error opening meta file for ms " + handle.data_->ms_path_);
  }

  meta_header_.Read(meta_file);
  std::vector<char> ms_path(meta_header_.filename_length + 1, char(0));
  meta_file.read(ms_path.data(), meta_header_.filename_length);
  Logger::Info << "Opening reordered part " << part_index << " for " << ms_path.data() << '\n';
  std::string part_prefix =
      GetPartPrefix(ms_path.data(), part_index, polarization, handle.data_->temporary_directory_);

  std::ifstream data_file(part_prefix + ".tmp", std::ios::in);
  if (!data_file.good())
    throw std::runtime_error("Error opening temporary data file '" + part_prefix + ".tmp'");
  part_header_.Read(data_file);
  if (!data_file.good())
    throw std::runtime_error("Error reading header from file '" + part_prefix + ".tmp'");

  if (part_header_.has_model) {
    model_file_path_ = part_prefix + "-m.tmp";
    model_file_ = MappedFile(model_file_path_);
  }
  meta_file.close();
  data_file.close();

  if (!meta_header_.data_desc_id.HasValue()) reader_.emplace(this);
}

ReorderedMsProvider::~ReorderedMsProvider() = default;

std::unique_ptr<MSReader> ReorderedMsProvider::MakeReader() {
  return std::make_unique<ReorderedMsReader>(this);
}

std::string ReorderedMsProvider::PartDescription() const {
  return handle_.data_->ms_path_ + ":" + std::to_string(part_index_);
}

void ReorderedMsProvider::NextOutputRow() {
  ++current_output_row_;
  if (reader_) {
    current_output_position_ +=
        SelectedBands()[reader_->CurrentDataDescId()].ChannelCount() * polarization_count_in_file_;
    reader_->NextInputRow();
  }
}

void ReorderedMsProvider::ResetWritePosition() {
  current_output_row_ = 0;
  current_output_position_ = 0;
  if (!meta_header_.data_desc_id.HasValue()) reader_.emplace(this);
}

void ReorderedMsProvider::WriteModel(const std::complex<float>* buffer, bool add_to_ms) {
#ifndef NDEBUG
  if (!part_header_.has_model) throw std::runtime_error("Reordered MS initialized without model");
#endif
  size_t n_channels;
  uint64_t position;
  if (reader_) {
    n_channels = SelectedBands()[reader_->CurrentDataDescId()].ChannelCount();
    position = current_output_position_;
  } else {
    n_channels = part_header_.max_channel_count;
    position = n_channels * polarization_count_in_file_ * current_output_row_;
  }

  std::complex<float>* model_write_ptr =
      reinterpret_cast<std::complex<float>*>(model_file_.Data()) + position;

  // In case the value was not sampled in this pass, it has been set to infinite
  // and should not overwrite the current value in the set.
  if (add_to_ms) {
    for (size_t i = 0; i != n_channels * polarization_count_in_file_; ++i) {
      if (std::isfinite(buffer[i].real())) model_write_ptr[i] += buffer[i];
    }
  } else {
    for (size_t i = 0; i != n_channels * polarization_count_in_file_; ++i) {
      if (std::isfinite(buffer[i].real())) model_write_ptr[i] = buffer[i];
    }
  }
}

void ReorderedMsProvider::ResetModelColumn() {
  model_file_.Close();
  std::remove(model_file_path_.c_str());
  schaapcommon::reordering::AllocateFile(model_file_path_, NMaxChannels() * NPolarizations() *
                                                               meta_header_.selected_row_count *
                                                               sizeof(std::complex<float>));
  model_file_ = MappedFile(model_file_path_);
}

double ReorderedMsProvider::Interval() { return GetMsInterval(*MS()); }

ObservationInfo ReorderedMsProvider::GetObservationInfo() {
  return ReadObservationInfo(*MS(), handle_.data_->selection_.FieldIds().front());
}

std::vector<std::string> ReorderedMsProvider::GetAntennaNames() {
  return wsclean::GetAntennaNames(MS()->antenna());
}

}  // namespace wsclean
