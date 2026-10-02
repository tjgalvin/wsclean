#ifndef TASK_MESSAGE_H
#define TASK_MESSAGE_H

#include <aocommon/io/serialistream.h>
#include <aocommon/io/serialostream.h>

namespace wsclean {

struct TaskMessage {
  enum class Type {
    kInvalid,
    /**
     * A 'kStart' message corresponds to a GriddingTaskManager::Start() call.
     * The main node sends it to all workers before sending Predict tasks, so
     * the workers can initialize their (local) writer locks.
     */
    kStart,
    /**
     * Using a 'kFinish' message, the main node signals the workers that all
     * work is done and the workers can exit.
     */
    kFinish,
    /** Message from the main node to a worker, containing a gridding task. */
    kGriddingRequest,
    /**
     * Message from a worker to the main node, containing the result of a
     * gridding task.
     */
    kGriddingResult,
  } type;
  union {
    // For kStart type.
    struct {
      uint32_t n_writer_groups;
      uint32_t n_parallel_tasks;
    };
    // For kGridding* types.
    uint64_t body_size;
  };
  static_assert(sizeof(body_size) == (sizeof(uint32_t) * 2));

  TaskMessage() : type(Type::kInvalid), body_size(0) {}
  TaskMessage(Type type_, size_t payload) : type(type_), body_size(payload) {}
  // Constructor for kStart
  TaskMessage(Type type_, uint32_t n_writer_groups_, uint32_t n_parallel_tasks_) : type(type_) {
    assert(type_ == Type::kStart);
    n_writer_groups = n_writer_groups_;
    n_parallel_tasks = n_parallel_tasks_;
  }

  constexpr static size_t kSerializedSize = 12;

  void Serialize(aocommon::SerialOStream& stream) const {
    stream.UInt32(static_cast<std::uint32_t>(type)).UInt64(body_size);
  }

  void Unserialize(aocommon::SerialIStream& stream) { stream.UInt32(type).UInt64(body_size); }
};

}  // namespace wsclean

#endif
