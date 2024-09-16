#include "reorderedmsreader.h"
#include "../reorderedmsprovider.h"

namespace wsclean {

ReorderedMsReader::ReorderedMsReader(ReorderedMsProvider* reordered_ms)
    : MSReader(reordered_ms),
      current_input_row_(0),
      read_ptr_row_offset_(0),
      meta_ptr_row_offset_(0),
      weight_ptr_row_offset_(0) {
  meta_file_.open(schaapcommon::reordering::GetMetaFilename(
                      reordered_ms->handle_.data_->ms_path_,
                      reordered_ms->handle_.data_->temporary_directory_,
                      reordered_ms->part_header_.data_desc_id),
                  std::ios::in);
  std::vector<char> ms_path(reordered_ms->meta_header_.filename_length + 1,
                            char(0));
  // meta and data header were read in ReorderedMs constructor
  meta_file_.seekg(schaapcommon::reordering::MetaHeader::BINARY_SIZE,
                   std::ios::beg);
  meta_file_.read(ms_path.data(), reordered_ms->meta_header_.filename_length);
  std::string part_prefix = schaapcommon::reordering::GetPartPrefix(
      ms_path.data(), reordered_ms->part_index_, reordered_ms->polarization_,
      reordered_ms->part_header_.data_desc_id,
      reordered_ms->handle_.data_->temporary_directory_);
  data_file_.open(part_prefix + ".tmp", std::ios::in);
  if (!data_file_.good())
    throw std::runtime_error("Error opening temporary data file in '" +
                             part_prefix + ".tmp'");
  data_file_.seekg(schaapcommon::reordering::PartHeader::BINARY_SIZE,
                   std::ios::beg);

  weight_file_.open(part_prefix + "-w.tmp", std::ios::in);
  if (!weight_file_.good())
    throw std::runtime_error("Error opening temporary data weight file '" +
                             part_prefix + "-w.tmp'");
}

bool ReorderedMsReader::CurrentRowAvailable() {
  const ReorderedMsProvider& reordered_ms =
      static_cast<const ReorderedMsProvider&>(*_msProvider);
  return current_input_row_ < reordered_ms.meta_header_.selected_row_count;
}

void ReorderedMsReader::NextInputRow() {
  const ReorderedMsProvider& reordered_ms =
      static_cast<const ReorderedMsProvider&>(*_msProvider);

  ++current_input_row_;
  if (current_input_row_ < reordered_ms.meta_header_.selected_row_count) {
    read_ptr_row_offset_ += 1;
    meta_ptr_row_offset_ += 1;
    weight_ptr_row_offset_ += 1;
  }
}

void ReorderedMsReader::ReadMeta(double& u, double& v, double& w) {
  if (meta_ptr_row_offset_ != 0)
    meta_file_.seekg(meta_ptr_row_offset_ *
                         schaapcommon::reordering::MetaRecord::BINARY_SIZE,
                     std::ios::cur);
  meta_ptr_row_offset_ = -1;

  schaapcommon::reordering::MetaRecord record;
  record.Read(meta_file_);
  u = record.u;
  v = record.v;
  w = record.w;
}

void ReorderedMsReader::ReadMeta(MSProvider::MetaData& meta_data) {
  if (meta_ptr_row_offset_ != 0)
    meta_file_.seekg(meta_ptr_row_offset_ *
                         schaapcommon::reordering::MetaRecord::BINARY_SIZE,
                     std::ios::cur);
  meta_ptr_row_offset_ = -1;

  schaapcommon::reordering::MetaRecord record;
  record.Read(meta_file_);
  meta_data.uInM = record.u;
  meta_data.vInM = record.v;
  meta_data.wInM = record.w;
  meta_data.fieldId = record.field_id;
  meta_data.antenna1 = record.antenna1;
  meta_data.antenna2 = record.antenna2;
  meta_data.time = record.time;
}

void ReorderedMsReader::ReadData(std::complex<float>* buffer) {
  const ReorderedMsProvider& reordered_ms =
      static_cast<const ReorderedMsProvider&>(*_msProvider);

  const int64_t n_visibilities = reordered_ms.part_header_.channel_count *
                                 reordered_ms.polarization_count_in_file_;
  if (read_ptr_row_offset_ != 0) {
    // Data file position was moved forward already, so seek back by one block
    data_file_.seekg(
        read_ptr_row_offset_ * (n_visibilities * sizeof(std::complex<float>)),
        std::ios::cur);
  }
  read_ptr_row_offset_ = -1;
#ifndef NDEBUG
  const size_t pos = size_t(data_file_.tellg()) -
                     schaapcommon::reordering::PartHeader::BINARY_SIZE;
  if (pos !=
      current_input_row_ * n_visibilities * sizeof(std::complex<float>)) {
    std::ostringstream s;
    s << "Not on right pos: " << pos << " instead of "
      << current_input_row_ * n_visibilities * sizeof(std::complex<float>)
      << " (row " << (pos / (n_visibilities * sizeof(std::complex<float>)))
      << " instead of " << current_input_row_ << ")";
    throw std::runtime_error(s.str());
  }
#endif
  data_file_.read(reinterpret_cast<char*>(buffer),
                  n_visibilities * sizeof(std::complex<float>));
}

void ReorderedMsReader::ReadModel(std::complex<float>* buffer) {
  const ReorderedMsProvider& reordered_ms =
      static_cast<ReorderedMsProvider&>(*_msProvider);

#ifndef NDEBUG
  if (!reordered_ms.part_header_.has_model)
    throw std::runtime_error("Reordered MS initialized without model");
#endif
  const size_t row_length = reordered_ms.part_header_.channel_count *
                            reordered_ms.polarization_count_in_file_ *
                            sizeof(std::complex<float>);
  std::copy_n(reordered_ms.model_file_.Data() + row_length * current_input_row_,
              row_length, reinterpret_cast<char*>(buffer));
}

void ReorderedMsReader::ReadWeights(float* buffer) {
  const ReorderedMsProvider& reordered_ms =
      static_cast<const ReorderedMsProvider&>(*_msProvider);

  const int64_t n_visibilities = reordered_ms.part_header_.channel_count *
                                 reordered_ms.polarization_count_in_file_;
  if (weight_ptr_row_offset_ != 0) {
    // jump to the previous block of weights
    weight_file_.seekg(
        weight_ptr_row_offset_ * (n_visibilities * sizeof(float)),
        std::ios::cur);
  }
  weight_file_.read(reinterpret_cast<char*>(buffer),
                    n_visibilities * sizeof(float));
  weight_ptr_row_offset_ = -1;
}

void ReorderedMsReader::WriteImagingWeights(const float* buffer) {
  const ReorderedMsProvider& reordered_ms =
      static_cast<const ReorderedMsProvider&>(*_msProvider);

  if (imaging_weights_file_ == nullptr) {
    std::string part_prefix = schaapcommon::reordering::GetPartPrefix(
        reordered_ms.handle_.data_->ms_path_, reordered_ms.part_index_,
        reordered_ms.polarization_, reordered_ms.part_header_.data_desc_id,
        reordered_ms.handle_.data_->temporary_directory_);
    imaging_weights_file_.reset(
        new std::fstream(part_prefix + "-imgw.tmp",
                         std::ios::in | std::ios::out | std::ios::binary));
  }
  const size_t n_vis = reordered_ms.part_header_.channel_count *
                       reordered_ms.polarization_count_in_file_;
  imaging_weight_buffer_.resize(n_vis);
  const size_t chunkSize = n_vis * sizeof(float);
  imaging_weights_file_->seekg(chunkSize * current_input_row_, std::ios::beg);
  imaging_weights_file_->read(
      reinterpret_cast<char*>(imaging_weight_buffer_.data()),
      n_vis * sizeof(float));
  for (size_t i = 0; i != n_vis; ++i) {
    if (std::isfinite(buffer[i])) imaging_weight_buffer_[i] = buffer[i];
  }
  imaging_weights_file_->seekp(chunkSize * current_input_row_, std::ios::beg);
  imaging_weights_file_->write(
      reinterpret_cast<const char*>(imaging_weight_buffer_.data()),
      n_vis * sizeof(float));
}

}  // namespace wsclean
