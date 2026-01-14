#include "mpischeduler.h"

#include <algorithm>
#include <cassert>
#include <memory>

#include "griddingresult.h"

#include "../main/settings.h"

#include "../distributed/mpibig.h"
#include "../distributed/taskmessage.h"

#include <aocommon/logger.h>
#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

#include <mpi.h>

using aocommon::Logger;

namespace wsclean {

namespace {
constexpr int kMainNode = 0;
constexpr int kTag = 0;
constexpr int kSlotsPerNode = 1;
}  // namespace

MPIScheduler::MPIScheduler(const Settings& settings)
    : GriddingTaskManager(settings), local_scheduler_(settings) {
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Set enough room for one pre-emptive send over and above what a node is
  // capable of processing in parallel. As 'available_room_' can become
  // negative; for compound tasks this will allow pre-emptive sending of an
  // additional compound task even if its size is much larger than 1.
  room_per_node_ = settings.parallelGridding + 1;
  available_room_.assign(world_size, room_per_node_);
  node_is_sending_.assign(world_size, false);
  if (!settings.masterDoesWork) {
    available_room_[0] = 0;
  }
  local_scheduler_.SetWriterLockManager(*this);

  // Add 1 as we want 1 thread per MPI node and won't be making use of the
  // current/main thread.
  send_thread_pool_.SetNThreads(GetSettings().nMpiNodes + 1);
  send_thread_pool_.StartParallelExecution([&](size_t) {
    std::function<void()> operation;
    while (send_task_queue_.Pop(operation)) {
      operation();
    }
  });
}

void MPIScheduler::Run(GriddingTask&& task,
                       std::function<void(GriddingResult&)> finish_callback) {
  if (!is_running_) {
    is_finishing_ = false;
    if (available_room_.size() > 1)
      receive_thread_ = std::thread([&]() { ReceiveLoop(); });
    is_running_ = true;
  }
  const size_t node = GetNode(task, std::move(finish_callback));
  // Use shared_ptr to work around std::function restriction on moveable types.
  // With C++23 this can become std::move_only_function instead.
  std::shared_ptr<GriddingTask> shared_task =
      std::make_shared<GriddingTask>(std::move(task));
  {
    std::unique_lock<std::mutex> lock(mutex_);
    while (node_is_sending_[node] == true) {
      notify_.wait(lock);
    }
    node_is_sending_[node] = true;
    send_task_queue_.Emplace([=, this]() {
      SendToNode(node, std::move(shared_task));
      node_is_sending_[node] = false;
      notify_.notify_all();
    });
  }

  std::lock_guard<std::mutex> lock(mutex_);
  ProcessReadyList_UNSYNCHRONIZED();
}

void MPIScheduler::Finish() {
  if (is_running_) {
    // Subtract 1 as the current/main thread never waits for tasks so won't be
    // registered as idling.
    send_task_queue_.WaitForIdle(send_thread_pool_.NThreads() - 1);

    Logger::Debug << "Finishing scheduler.\n";
    local_scheduler_.Finish();
    {
      std::unique_lock<std::mutex> lock(mutex_);
      is_finishing_ = true;
      notify_.notify_all();
    }

    // As long as receive tasks are running, wait and keep processing
    // the ready list
    Logger::Debug << "Waiting for all tasks to finish processing.\n";
    {
      std::unique_lock<std::mutex> lock(mutex_);
      ProcessReadyList_UNSYNCHRONIZED();
      while (AWorkerIsRunning_UNSYNCHRONIZED()) {
        notify_.wait(lock);
        ProcessReadyList_UNSYNCHRONIZED();
      }
    }

    if (available_room_.size() > 1) receive_thread_.join();

    is_running_ = false;

    // The while loop above ignores the work thread, which might
    // be gridding on the master node. Therefore, the master thread
    // might have added an item to the ready list. Therefore,
    // the ready list should once more be processed.
    // A lock is no longer required, because all threads have stopped.
    ProcessReadyList_UNSYNCHRONIZED();
  }
}

void MPIScheduler::Start(size_t n_writer_groups) {
  GriddingTaskManager::Start(n_writer_groups);

  const TaskMessage message(TaskMessage::Type::kStart, n_writer_groups);
  aocommon::SerialOStream message_stream;
  message.Serialize(message_stream);
  assert(message_stream.size() == TaskMessage::kSerializedSize);

  for (size_t rank = 1; rank < available_room_.size(); ++rank) {
    assert(available_room_[rank] > 0 &&
           available_room_[rank] == room_per_node_);
    MPI_Send(message_stream.data(), TaskMessage::kSerializedSize, MPI_BYTE,
             rank, kTag, MPI_COMM_WORLD);
  }

  if (GetSettings().masterDoesWork) {
    local_scheduler_.Start(n_writer_groups);
  }
}

void MPIScheduler::SendToNode(size_t node,
                              std::shared_ptr<GriddingTask> wrapped_task) {
  GriddingTask& task = *wrapped_task;

  if (node == 0) {
    Logger::Info << "Running gridding task " << task.unique_id
                 << " at main node.\n";

    local_scheduler_.Run(std::move(task), [this](GriddingResult& result) {
      Logger::Info << "Main node has finished gridding task "
                   << result.unique_id << ".\n";
      StoreResult(std::move(result), 0);
    });
  } else {
    aocommon::SerialOStream payload_stream;
    // To use MPI_Send_Big, a uint64_t need to be reserved
    payload_stream.UInt64(0);
    task.Serialize(payload_stream);

    const TaskMessage message(TaskMessage::Type::kGriddingRequest,
                              payload_stream.size());
    aocommon::SerialOStream task_message_stream;
    message.Serialize(task_message_stream);
    assert(task_message_stream.size() == TaskMessage::kSerializedSize);

    Logger::Info << "Sending gridding task " << task.unique_id << " to node "
                 << node << " (size: " << payload_stream.size() << ").\n";
    MPI_Send(task_message_stream.data(), task_message_stream.size(), MPI_BYTE,
             node, 0, MPI_COMM_WORLD);
    MPI_Send_Big(payload_stream.data(), payload_stream.size(), node, 0,
                 MPI_COMM_WORLD, GetSettings().maxMpiMessageSize);
  }
}

size_t MPIScheduler::GetNode(const GriddingTask& task,
                             std::function<void(GriddingResult&)>&& callback) {
  // Determine the target node using the channel to node mapping.
  const size_t node = GetSettings().channelToNode[task.outputChannelIndex];

  // Wait until available_room_[node] becomes larger than 0.
  std::unique_lock<std::mutex> lock(mutex_);
  while (available_room_[node] <= 0) {
    notify_.wait(lock);
  }
  const size_t task_size = task.num_parallel_gridders_;
  available_room_[node] -= task_size;
  task_size_[task.unique_id] = task_size;
  notify_.notify_all();  // Notify receiveLoop(). It should stop waiting.

  // Store the callback function.
  assert(callbacks_.count(task.unique_id) == 0);
  callbacks_.emplace(task.unique_id, std::move(callback));

  return node;
}

void MPIScheduler::ReceiveLoop() {
  std::unique_lock<std::mutex> lock(mutex_);
  while (!is_finishing_ || AWorkerIsRunning_UNSYNCHRONIZED()) {
    if (!AWorkerIsRunning_UNSYNCHRONIZED()) {
      notify_.wait(lock);
    } else {
      lock.unlock();

      TaskMessage message;
      MPI_Status status;
      aocommon::UVector<unsigned char> buffer(TaskMessage::kSerializedSize);
      MPI_Recv(buffer.data(), TaskMessage::kSerializedSize, MPI_BYTE,
               MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
      aocommon::SerialIStream stream(std::move(buffer));
      message.Unserialize(stream);

      const int node = status.MPI_SOURCE;
      switch (message.type) {
        case TaskMessage::Type::kGriddingResult:
          ProcessGriddingResult(node, message.body_size);
          break;
        default:
          throw std::runtime_error("Invalid message sent by node " +
                                   std::to_string(node));
      }

      lock.lock();
    }
  }
  Logger::Info << "All worker nodes have finished their gridding tasks.\n";
}

void MPIScheduler::ProcessReadyList_UNSYNCHRONIZED() {
  while (!ready_list_.empty()) {
    // Call the callback for this finished task
    GriddingResult& result = ready_list_.back();
    // Copy the task id, since callbacks may adjust the result.
    const size_t task_id = result.unique_id;
    callbacks_[task_id](result);
    ready_list_.pop_back();
    callbacks_.erase(task_id);
  }
}

bool MPIScheduler::AWorkerIsRunning_UNSYNCHRONIZED() {
  for (size_t i = 1; i != available_room_.size(); ++i) {
    if (available_room_[i] < room_per_node_) {
      return true;
    }
  }
  return false;
}

void MPIScheduler::ProcessGriddingResult(size_t node, size_t body_size) {
  aocommon::UVector<unsigned char> buffer(body_size);
  MPI_Status status;
  MPI_Recv_Big(buffer.data(), body_size, node, 0, MPI_COMM_WORLD, &status,
               GetSettings().maxMpiMessageSize);
  GriddingResult result;
  aocommon::SerialIStream stream(std::move(buffer));
  stream.UInt64();  // storage for MPI_Recv_Big
  result.Unserialize(stream);
  StoreResult(std::move(result), node);
}

void MPIScheduler::StoreResult(GriddingResult&& result, int node) {
  std::lock_guard<std::mutex> lock(mutex_);
  const size_t task_id = result.unique_id;
  available_room_[node] += task_size_[task_id];
  task_size_.erase(task_id);
  ready_list_.emplace_back(std::move(result));
  notify_.notify_all();
}

}  // namespace wsclean
