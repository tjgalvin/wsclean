#include "mpischeduler.h"

#include "griddingresult.h"

#include "../io/logger.h"

#include "../main/settings.h"

#include "../distributed/mpibig.h"
#include "../distributed/taskmessage.h"

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

#include <mpi.h>

#include <cassert>

MPIScheduler::MPIScheduler(const Settings &settings)
    : GriddingTaskManager(settings),
      _masterDoesWork(settings.masterDoesWork),
      _isRunning(false),
      _isFinishing(false),
      _mutex(),
      _receiveThread(),
      _workThread(),
      _readyList(),
      _nodes(),
      _writerLockQueues() {
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  _nodes.assign(world_size,
                std::make_pair(NodeState::kAvailable,
                               std::function<void(GriddingResult &)>()));
  if (!settings.masterDoesWork && world_size <= 1)
    throw std::runtime_error(
        "Master was told not to work, but no other workers available");
}

MPIScheduler::~MPIScheduler() { Finish(); }

void MPIScheduler::Run(GriddingTask &&task,
                       std::function<void(GriddingResult &)> finishCallback) {
  if (!_isRunning) {
    _isFinishing = false;
    if (_nodes.size() > 1)
      _receiveThread = std::thread([&]() { receiveLoop(); });
    _isRunning = true;
  }
  send(std::move(task), finishCallback);

  std::lock_guard<std::mutex> lock(_mutex);
  processReadyList_UNSYNCHRONIZED();
}

void MPIScheduler::Finish() {
  Logger::Info << "Finishing scheduler.\n";
  if (_isRunning) {
    std::unique_lock<std::mutex> lock(_mutex);
    _isFinishing = true;
    _notify.notify_all();

    // As long as receive tasks are running, wait and keep processing
    // the ready list
    processReadyList_UNSYNCHRONIZED();
    while (receiveTasksAreRunning_UNSYNCHRONIZED()) {
      _notify.wait(lock);
      processReadyList_UNSYNCHRONIZED();
    }

    lock.unlock();

    if (_nodes.size() > 1) _receiveThread.join();

    if (_workThread.joinable()) _workThread.join();
    _isRunning = false;

    // The while loop above ignores the work thread, which might
    // be gridding on the master node. Therefore, the master thread
    // might have added an item to the ready list. Therefore,
    // the ready list should once more be processed.
    // A lock is no longer required, because all threads have stopped.
    processReadyList_UNSYNCHRONIZED();
  }
}

void MPIScheduler::Start(size_t nWriterGroups) {
  GriddingTaskManager::Start(nWriterGroups);
  _writerLockQueues.resize(nWriterGroups);
}

WriterLockManager::LockGuard MPIScheduler::GetLock(size_t writerGroupIndex) {
  return LockGuard(_writerLock);
}

void MPIScheduler::runTaskOnNode0(GriddingTask &&task) {
  GriddingResult result = RunDirect(std::move(task));
  Logger::Info << "Master node has finished a gridding task.\n";
  std::unique_lock<std::mutex> lock(_mutex);
  _readyList.emplace_back(std::move(result), _nodes[0].second);
  _nodes[0].first = NodeState::kAvailable;
  lock.unlock();
  _notify.notify_all();
}

void MPIScheduler::send(GriddingTask &&task,
                        const std::function<void(GriddingResult &)> &callback) {
  int node = findAndSetNodeState(NodeState::kAvailable,
                                 std::make_pair(NodeState::kBusy, callback));
  Logger::Info << "Sending gridding task to node : " << node << '\n';

  if (node == 0) {
    if (_workThread.joinable()) _workThread.join();
    _workThread =
        std::thread(&MPIScheduler::runTaskOnNode0, this, std::move(task));
  } else {
    aocommon::SerialOStream payloadStream;
    // To use MPI_Send_Big, a uint64_t need to be reserved
    payloadStream.UInt64(0);
    task.Serialize(payloadStream);

    TaskMessage message;
    message.type = TaskMessage::Type::kGriddingRequest;
    message.bodySize = payloadStream.size();

    aocommon::SerialOStream taskMessageStream;
    message.Serialize(taskMessageStream);
    assert(taskMessageStream.size() == TaskMessage::kSerializedSize);

    MPI_Send(taskMessageStream.data(), taskMessageStream.size(), MPI_BYTE, node,
             0, MPI_COMM_WORLD);
    MPI_Send_Big(payloadStream.data(), payloadStream.size(), node, 0,
                 MPI_COMM_WORLD);
  }
}

int MPIScheduler::findAndSetNodeState(
    MPIScheduler::NodeState currentState,
    std::pair<MPIScheduler::NodeState, std::function<void(GriddingResult &)>>
        newState) {
  std::unique_lock<std::mutex> lock(_mutex);
  do {
    size_t iterEnd = _masterDoesWork ? _nodes.size() : _nodes.size() - 1;
    for (size_t i = 0; i != iterEnd; ++i) {
      const int node = _nodes.size() - i - 1;
      if (_nodes[node].first == currentState) {
        _nodes[node] = newState;
        _notify.notify_all();
        return node;
      }
    }
    _notify.wait(lock);
  } while (true);
}

void MPIScheduler::receiveLoop() {
  std::unique_lock<std::mutex> lock(_mutex);
  while (!_isFinishing || receiveTasksAreRunning_UNSYNCHRONIZED()) {
    if (!receiveTasksAreRunning_UNSYNCHRONIZED()) {
      _notify.wait(lock);
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
          processGriddingResult(node, message.bodySize);
          break;
        case TaskMessage::Type::kLockRequest:
          processLockRequest(node, message.lockId);
          break;
        case TaskMessage::Type::kLockRelease:
          processLockRelease(node, message.lockId);
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

void MPIScheduler::processReadyList_UNSYNCHRONIZED() {
  while (!_readyList.empty()) {
    // Call the callback for this finished task
    _readyList.back().second(_readyList.back().first);
    _readyList.pop_back();
  }
}

bool MPIScheduler::receiveTasksAreRunning_UNSYNCHRONIZED() {
  for (size_t i = 1; i != _nodes.size(); ++i)
    if (_nodes[i].first == NodeState::kBusy) return true;
  return false;
}

void MPIScheduler::processGriddingResult(int node, size_t bodySize) {
  aocommon::UVector<unsigned char> buffer(bodySize);
  MPI_Status status;
  MPI_Recv_Big(buffer.data(), bodySize, node, 0, MPI_COMM_WORLD, &status);
  GriddingResult result;
  aocommon::SerialIStream stream(std::move(buffer));
  stream.UInt64();  // storage for MPI_Recv_Big
  result.Unserialize(stream);

  std::lock_guard<std::mutex> lock(_mutex);
  _readyList.emplace_back(std::move(result), _nodes[node].second);
  _nodes[node].first = NodeState::kAvailable;
  _notify.notify_all();
}

void MPIScheduler::processLockRequest(int node, size_t lockId) {
  if (lockId >= _writerLockQueues.size()) {
    throw std::runtime_error("Node " + std::to_string(node) +
                             " requests invalid lock " +
                             std::to_string(lockId));
  }

  std::unique_lock<std::mutex> lock(_mutex);
  if (_nodes[node].first != NodeState::kBusy) {
    throw std::runtime_error("Non-busy node " + std::to_string(node) +
                             " requests lock " + std::to_string(lockId));
  }
  _writerLockQueues[lockId].PushBack(node);
  if (_writerLockQueues[lockId].Size() == 1) {
    lock.unlock();
    grantLock(node, lockId);
  }
}

void MPIScheduler::processLockRelease(int node, size_t lockId) {
  if (lockId >= _writerLockQueues.size()) {
    throw std::runtime_error("Node " + std::to_string(node) +
                             " releases invalid lock id " +
                             std::to_string(lockId));
  }

  std::unique_lock<std::mutex> lock(_mutex);
  if (_writerLockQueues[lockId].Empty() ||
      _writerLockQueues[lockId][0] != node) {
    throw std::runtime_error("Node " + std::to_string(node) +
                             " releases not-granted lock id " +
                             std::to_string(lockId));
  }

  _writerLockQueues[lockId].PopFront();
  if (!_writerLockQueues[lockId].Empty()) {
    const int waiting_node = _writerLockQueues[lockId][0];
    lock.unlock();
    grantLock(waiting_node, lockId);
  }
}

void MPIScheduler::grantLock(int node, size_t lockId) {
  const TaskMessage message(TaskMessage::Type::kLockGrant, lockId);
  aocommon::SerialOStream taskMessageStream;
  message.Serialize(taskMessageStream);
  assert(taskMessageStream.size() == TaskMessage::kSerializedSize);

  // Using asynchronous MPI_ISend is possible here, however, a synchronous
  // MPI_Send is much simpler. Since the message is small and the receiver
  // is already waiting for the message, the overhead of synchronous
  // communication should be limited.
  MPI_Send(taskMessageStream.data(), taskMessageStream.size(), MPI_BYTE, node,
           0, MPI_COMM_WORLD);
}
