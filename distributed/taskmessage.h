#ifndef TASK_MESSAGE_H
#define TASK_MESSAGE_H

#include <aocommon/io/serialistream.h>
#include <aocommon/io/serialostream.h>

struct TaskMessage {
  enum class Type {
    kInvalid,
    kFinish,
    kGriddingRequest,
    kGriddingResult,
    kLockRequest,
    kLockGrant,
    kLockRelease
  } type;
  union {
    size_t bodySize;  // For kGridding* types.
    size_t lockId;    // For kLock* types.
  };

  TaskMessage() : type(Type::kInvalid), bodySize(0) {}
  TaskMessage(Type type_, size_t payload) : type(type_), bodySize(payload) {}

  constexpr static size_t kSerializedSize = 12;

  void Serialize(aocommon::SerialOStream& stream) const {
    stream.UInt32(static_cast<std::uint32_t>(type)).UInt64(bodySize);
  }

  void Unserialize(aocommon::SerialIStream& stream) {
    stream.UInt32(type).UInt64(bodySize);
  }
};

#endif
