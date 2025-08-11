#include "reorderedmsreader.h"
#include "../reorderedmsprovider.h"

namespace wsclean {

ReorderedMsReader::ReorderedMsReader(ReorderedMsProvider* reordered_ms)
    : MSReader(reordered_ms),
      current_input_row_(0),
      read_ptr_row_offset_(0),
      meta_ptr_row_offset_(0),
      weight_ptr_row_offset_(0) {
  const size_t meta_file_index =
      reordered_ms->handle_.data_->metadata_indices_[reordered_ms->part_index_];
  const std::string meta_filename = schaapcommon::reordering::GetMetaFilename(
      reordered_ms->handle_.data_->ms_path_,
      reordered_ms->handle_.data_->temporary_directory_, meta_file_index);
  meta_file_.open(meta_filename, std::ios::in);
  std::vector<char> ms_path(reordered_ms->meta_header_.filename_length + 1,
                            char(0));
  // meta and data header were read in ReorderedMs constructor
  meta_file_.seekg(schaapcommon::reordering::MetaHeader::BINARY_SIZE,
                   std::ios::beg);
  meta_file_.read(ms_path.data(), reordered_ms->meta_header_.filename_length);
  if (!meta_file_.good())
    throw std::runtime_error("Error opening temporary metadata file: '" +
                             meta_filename);
  const std::string part_prefix = schaapcommon::reordering::GetPartPrefix(
      ms_path.data(), reordered_ms->part_index_, reordered_ms->polarization_,
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

inline uint64_t ReorderedMsReader::NVisibilitiesPerRow() const {
  ReorderedMsProvider& reordered_ms =
      static_cast<ReorderedMsProvider&>(*ms_provider_);
  const size_t n_channels =
      reordered_ms.SelectedBands()[reordered_ms.data_desc_id_].ChannelCount();
  return n_channels * reordered_ms.polarization_count_in_file_;
}

bool ReorderedMsReader::CurrentRowAvailable() {
  const ReorderedMsProvider& reordered_ms =
      static_cast<const ReorderedMsProvider&>(*ms_provider_);
  return current_input_row_ < reordered_ms.meta_header_.selected_row_count;
}

void ReorderedMsReader::NextInputRow() {
  const ReorderedMsProvider& reordered_ms =
      static_cast<const ReorderedMsProvider&>(*ms_provider_);

  ++current_input_row_;
  if (current_input_row_ < reordered_ms.meta_header_.selected_row_count) {
    read_ptr_row_offset_ += 1;
    meta_ptr_row_offset_ += 1;
    weight_ptr_row_offset_ += 1;
  }
}
void ReorderedMsReader::ReadMeta(MSProvider::MetaData& metadata) {
  if (meta_ptr_row_offset_ != 0)
    meta_file_.seekg(meta_ptr_row_offset_ *
                         schaapcommon::reordering::MetaRecord::BINARY_SIZE,
                     std::ios::cur);
  meta_ptr_row_offset_ = -1;

  schaapcommon::reordering::MetaRecord record;
  record.Read(meta_file_);
  metadata.u_in_m = record.u;
  metadata.v_in_m = record.v;
  metadata.w_in_m = record.w;
  metadata.time = record.time;
  const ReorderedMsProvider& reordered_ms =
      static_cast<const ReorderedMsProvider&>(*ms_provider_);
  metadata.data_desc_id = reordered_ms.data_desc_id_;
  metadata.field_id = record.field_id;
  metadata.antenna1 = record.antenna1;
  metadata.antenna2 = record.antenna2;
}

void ReorderedMsReader::ReadData(std::complex<float>* buffer) {
  const int64_t n_visibilities = NVisibilitiesPerRow();
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
  ReorderedMsProvider& reordered_ms =
      static_cast<ReorderedMsProvider&>(*ms_provider_);

#ifndef NDEBUG
  if (!reordered_ms.part_header_.has_model)
    throw std::runtime_error("Reordered MS initialized without model");
#endif
  const size_t row_length = NVisibilitiesPerRow() * sizeof(std::complex<float>);
  std::copy_n(reordered_ms.model_file_.Data() + row_length * current_input_row_,
              row_length, reinterpret_cast<char*>(buffer));
}

void ReorderedMsReader::ReadWeights(float* buffer) {
  const int64_t n_visibilities = NVisibilitiesPerRow();
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
  ReorderedMsProvider& reordered_ms =
      static_cast<ReorderedMsProvider&>(*ms_provider_);

  if (imaging_weights_file_ == nullptr) {
    std::string part_prefix = schaapcommon::reordering::GetPartPrefix(
        reordered_ms.handle_.data_->ms_path_, reordered_ms.part_index_,
        reordered_ms.polarization_,
        reordered_ms.handle_.data_->temporary_directory_);
    imaging_weights_file_.reset(
        new std::fstream(part_prefix + "-imgw.tmp",
                         std::ios::in | std::ios::out | std::ios::binary));
  }
  const size_t n_vis = NVisibilitiesPerRow();
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
