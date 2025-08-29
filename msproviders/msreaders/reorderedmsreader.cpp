#include "reorderedmsreader.h"
#include "../reorderedmsprovider.h"

namespace wsclean {

ReorderedMsReader::ReorderedMsReader(ReorderedMsProvider* reordered_ms)
    : MSReader(reordered_ms) {
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

  const size_t n_max_visibilities =
      reordered_ms->NMaxChannels() * reordered_ms->polarization_count_in_file_;
  data_buffer_.resize(n_max_visibilities);
  weight_buffer_.resize(n_max_visibilities);

  if (reordered_ms->meta_header_.data_desc_id) {
    metadata_.data_desc_id = *reordered_ms->meta_header_.data_desc_id;
  }
  CacheMeta();
}

inline uint64_t ReorderedMsReader::NVisibilitiesPerRow(
    size_t data_desc_id) const {
  ReorderedMsProvider& reordered_ms =
      static_cast<ReorderedMsProvider&>(*ms_provider_);
  const size_t n_channels =
      reordered_ms.SelectedBands()[data_desc_id].ChannelCount();
  return n_channels * reordered_ms.polarization_count_in_file_;
}

bool ReorderedMsReader::CurrentRowAvailable() {
  const ReorderedMsProvider& reordered_ms =
      static_cast<const ReorderedMsProvider&>(*ms_provider_);
  return current_input_row_ < reordered_ms.meta_header_.selected_row_count;
}

void ReorderedMsReader::NextInputRow() {
  current_value_position_ += NVisibilitiesPerRow(metadata_.data_desc_id);
  ++current_input_row_;
  has_data_ = false;
  has_weights_ = false;
  CacheMeta();
}

void ReorderedMsReader::ReadMeta(MSProvider::MetaData& metadata) {
  metadata = metadata_;
}

void ReorderedMsReader::CacheMeta() {
  if (Provider().IsRegular()) {
    schaapcommon::reordering::MetaRecord record;
    record.Read(meta_file_);
    metadata_.u_in_m = record.u;
    metadata_.v_in_m = record.v;
    metadata_.w_in_m = record.w;
    metadata_.time = record.time;
    // data_desc_id is already set in constructor
    metadata_.field_id = record.field_id;
    metadata_.antenna1 = record.antenna1;
    metadata_.antenna2 = record.antenna2;
  } else {
    schaapcommon::reordering::BdaMetaRecord record;
    record.Read(meta_file_);
    metadata_.u_in_m = record.u;
    metadata_.v_in_m = record.v;
    metadata_.w_in_m = record.w;
    metadata_.time = record.time;
    metadata_.data_desc_id = record.data_desc_id;
    metadata_.field_id = record.field_id;
    metadata_.antenna1 = record.antenna1;
    metadata_.antenna2 = record.antenna2;
  }
}

void ReorderedMsReader::ReadData(std::complex<float>* buffer) {
  const int64_t n_visibilities = NVisibilitiesPerRow(metadata_.data_desc_id);
  if (!has_data_) {
    data_file_.read(reinterpret_cast<char*>(data_buffer_.data()),
                    n_visibilities * sizeof(std::complex<float>));
    has_data_ = true;
  }
  std::copy_n(data_buffer_.data(), n_visibilities, buffer);
}

void ReorderedMsReader::ReadModel(std::complex<float>* buffer) {
  ReorderedMsProvider& reordered_ms =
      static_cast<ReorderedMsProvider&>(*ms_provider_);
  assert(reordered_ms.part_header_.has_model);
  const int64_t buffer_size =
      NVisibilitiesPerRow(metadata_.data_desc_id) * sizeof(std::complex<float>);
  std::copy_n(reordered_ms.model_file_.Data() +
                  current_value_position_ * sizeof(std::complex<float>),
              buffer_size, reinterpret_cast<char*>(buffer));
}

void ReorderedMsReader::ReadWeights(float* buffer) {
  const int64_t n_visibilities = NVisibilitiesPerRow(metadata_.data_desc_id);
  if (!has_weights_) {
    weight_file_.read(reinterpret_cast<char*>(weight_buffer_.data()),
                      n_visibilities * sizeof(float));
    has_weights_ = true;
  }
  std::copy_n(weight_buffer_.data(), n_visibilities, buffer);
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
  const size_t n_vis = NVisibilitiesPerRow(metadata_.data_desc_id);
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
