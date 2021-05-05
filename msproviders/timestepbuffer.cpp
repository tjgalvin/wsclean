#include "timestepbuffer.h"
#include "msreaders/timestepbufferreader.h"

std::unique_ptr<MSReader> TimestepBuffer::MakeReader() {
  std::unique_ptr<MSReader> reader(new TimestepBufferReader(this));
  return reader;
}