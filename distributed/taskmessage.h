#ifndef TASK_MESSAGE_H
#define TASK_MESSAGE_H

#include <aocommon/io/serialistream.h>
#include <aocommon/io/serialostream.h>

struct TaskMessage {
  enum Type { Finish, GriddingRequest, GriddingResult } type;
  size_t bodySize = 0;

  constexpr static size_t kSerializedSize = 12;
  void Serialize(aocommon::SerialOStream& stream) const {
    stream.UInt32(type).UInt64(bodySize);
  }

  void Unserialize(aocommon::SerialIStream& stream) {
    stream.UInt32(type).UInt64(bodySize);
  }
};

#endif
