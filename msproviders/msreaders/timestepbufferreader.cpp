#include "timestepbufferreader.h"

namespace wsclean {

TimestepBufferReader::TimestepBufferReader(TimestepBuffer* timestep_buffer)
    : MSReader(timestep_buffer),
      ms_reader_(timestep_buffer->ms_provider_->MakeReader()),
      buffer_position_(0) {
  readTimeblock();
}

bool TimestepBufferReader::CurrentRowAvailable() {
  return !buffer_.empty() || ms_reader_->CurrentRowAvailable();
}

void TimestepBufferReader::NextInputRow() {
  ++buffer_position_;
  if (buffer_position_ == buffer_.size()) {
    readTimeblock();
  }
}

void TimestepBufferReader::ReadMeta(double& u, double& v, double& w) {
  MSProvider::MetaData& m = buffer_[buffer_position_].metadata;
  u = m.u_in_m;
  v = m.v_in_m;
  w = m.w_in_m;
}

void TimestepBufferReader::ReadMeta(MSProvider::MetaData& metadata) {
  metadata = buffer_[buffer_position_].metadata;
}

void TimestepBufferReader::ReadData(std::complex<float>* buffer) {
  std::copy(buffer_[buffer_position_].data.begin(),
            buffer_[buffer_position_].data.end(), buffer);
}

void TimestepBufferReader::ReadModel(std::complex<float>* buffer) {
  std::copy(buffer_[buffer_position_].model.begin(),
            buffer_[buffer_position_].model.end(), buffer);
}

void TimestepBufferReader::ReadWeights(float* buffer) {
  std::copy(buffer_[buffer_position_].weights.begin(),
            buffer_[buffer_position_].weights.end(), buffer);
}

void TimestepBufferReader::WriteImagingWeights(const float* buffer) {
  ms_reader_->WriteImagingWeights(buffer);
}

void TimestepBufferReader::readTimeblock() {
  // Beware that the ms_provider_ data member is a TimestepBuffer,
  // which in turn has its own ms_provider_
  TimestepBuffer& tstepbuffer = static_cast<TimestepBuffer&>(*ms_provider_);

  buffer_position_ = 0;
  buffer_.clear();
  MSProvider::MetaData metadata;
  size_t dataSize = tstepbuffer.ms_provider_->NPolarizations() *
                    tstepbuffer.ms_provider_->NChannels();

  if (ms_reader_->CurrentRowAvailable()) {
    ms_reader_->ReadMeta(metadata);
    double blockTime = metadata.time, curTime = blockTime;
    size_t writePos = 0;
    do {
      if (buffer_.size() <= writePos) {
        buffer_.emplace_back();
        TimestepBuffer::RowData& newRow = buffer_.back();
        newRow.data.resize(dataSize);
        if (tstepbuffer.read_model_) newRow.model.resize(dataSize);
        newRow.weights.resize(dataSize);
      }
      TimestepBuffer::RowData& row = buffer_[writePos];
      row.metadata = metadata;
      ms_reader_->ReadData(row.data.data());
      if (tstepbuffer.read_model_) ms_reader_->ReadModel(row.model.data());
      ms_reader_->ReadWeights(row.weights.data());
      row.row_id = ms_reader_->RowId();

      ms_reader_->NextInputRow();
      ++writePos;
      if (ms_reader_->CurrentRowAvailable()) {
        ms_reader_->ReadMeta(metadata);
        curTime = metadata.time;
      }
    } while (ms_reader_->CurrentRowAvailable() && blockTime == curTime);
    buffer_.resize(writePos);
  }
}

}  // namespace wsclean
