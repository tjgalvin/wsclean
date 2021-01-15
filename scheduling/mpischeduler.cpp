#include "mpischeduler.h"

#include "griddingresult.h"

#include "../io/logger.h"

#include "../distributed/mpibig.h"
#include "../distributed/taskmessage.h"

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>

#include <mpi.h>

MPIScheduler::MPIScheduler(const class Settings &settings)
    : GriddingTaskManager(settings), _isRunning(false), _isSendFinished(false) {
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  _nodes.assign(
      world_size,
      std::make_pair(AvailableNode, std::function<void(GriddingResult &)>()));
  _taskList.resize(world_size);
}

MPIScheduler::~MPIScheduler() { Finish(); }

void MPIScheduler::Run(GriddingTask &&task,
                       std::function<void(GriddingResult &)> finishCallback) {
  if (!_isRunning) {
    _taskList.clear();
    _isSendFinished = false;
    _isFinishing = false;
    _sendThread = std::thread([&]() { sendLoop(); });
    _receiveThread = std::thread([&]() { receiveLoop(); });
    _isRunning = true;
  }
  _taskList.write(
      std::pair<GriddingTask, std::function<void(GriddingResult &)>>(
          std::move(task), finishCallback));

  std::unique_lock<std::mutex> lock(_mutex);
  while (!_readyList.empty()) {
    // Call callbacks for any finished tasks
    _readyList.back().second(_readyList.back().first);
    _readyList.pop_back();
  }
  lock.unlock();
}

void MPIScheduler::Finish() {
  Logger::Info << "Finishing scheduler.\n";
  if (_isRunning) {
    std::unique_lock<std::mutex> lock(_mutex);
    _isFinishing = true;
    _notify.notify_all();
    lock.unlock();

    _taskList.write_end();
    _sendThread.join();
    _receiveThread.join();
    if (_workThread.joinable()) _workThread.join();
    _taskList.clear();

    while (!_readyList.empty()) {
      // Call callbacks for any finished tasks
      _readyList.back().second(_readyList.back().first);
      _readyList.pop_back();
    }

    _isRunning = false;
  }
}

void MPIScheduler::node0gridder(GriddingTask task) {
  GriddingResult result = RunDirect(std::move(task));
  Logger::Info << "Master node is done gridding.\n";
  std::unique_lock<std::mutex> lock(_mutex);
  _readyList.emplace_back(std::move(result), _nodes[0].second);
  _nodes[0].first = AvailableNode;
  lock.unlock();
  _notify.notify_all();
}

void MPIScheduler::sendLoop() {
  std::pair<GriddingTask, std::function<void(GriddingResult &)>> taskPair;
  while (_taskList.read(taskPair)) {
    const GriddingTask &task = taskPair.first;
    aocommon::SerialOStream stream;
    // To use MPI_Send_Big, a uint64_t need to be reserved
    stream.UInt64(0);
    task.Serialize(stream);

    int node = findAndSetNodeState(AvailableNode,
                                   std::make_pair(BusyNode, taskPair.second));
    Logger::Info << "Sending gridding task to : " << node << '\n';

    if (node == 0) {
      if (_workThread.joinable()) _workThread.join();
      _workThread = std::thread(&MPIScheduler::node0gridder, this,
                                std::move(taskPair.first));
    } else {
      TaskMessage message;
      message.type = TaskMessage::GriddingRequest;
      message.bodySize = stream.size();
      MPI_Send(&message, sizeof(TaskMessage), MPI_BYTE, node, 0,
               MPI_COMM_WORLD);
      MPI_Send_Big(stream.data(), stream.size(), node, 0, MPI_COMM_WORLD);
    }
  }

  std::unique_lock<std::mutex> lock(_mutex);
  _isSendFinished = true;
}

int MPIScheduler::findAndSetNodeState(
    MPIScheduler::NodeState currentState,
    std::pair<MPIScheduler::NodeState, std::function<void(GriddingResult &)>>
        newState) {
  std::unique_lock<std::mutex> lock(_mutex);
  do {
    for (size_t i = 0; i != _nodes.size(); ++i) {
      int node = _nodes.size() - i - 1;
      if (_nodes[node].first == currentState) {
        _notify.notify_all();
        _nodes[node] = newState;
        return node;
      }
    }
    _notify.wait(lock);
  } while (true);
}

void MPIScheduler::receiveLoop() {
  std::unique_lock<std::mutex> lock(_mutex);
  while (!_isFinishing || anyReceiveTasks_NeedLock()) {
    if (!anyReceiveTasks_NeedLock()) {
      _notify.wait(lock);
    } else {
      lock.unlock();

      TaskMessage message;
      MPI_Status status;
      MPI_Recv(&message, sizeof(TaskMessage), MPI_BYTE, MPI_ANY_SOURCE, 0,
               MPI_COMM_WORLD, &status);
      int node = status.MPI_SOURCE;
      if (message.type != TaskMessage::GriddingResult)
        throw std::runtime_error("Invalid message sent by node " +
                                 std::to_string(node));

      aocommon::UVector<unsigned char> buffer(message.bodySize);
      MPI_Recv_Big(buffer.data(), message.bodySize, node, 0, MPI_COMM_WORLD,
                   &status);
      GriddingResult result;
      aocommon::SerialIStream stream(std::move(buffer));
      stream.UInt64();  // storage for MPI_Recv_Big
      result.Unserialize(stream);

      lock.lock();
      _readyList.emplace_back(std::move(result), _nodes[node].second);
      _nodes[node].first = AvailableNode;
      lock.unlock();

      _notify.notify_all();

      lock.lock();
    }
  }
  Logger::Info << "Receive loop finished.\n";
}

bool MPIScheduler::anyReceiveTasks_NeedLock() {
  if (!_isSendFinished) return true;
  for (size_t i = 1; i != _nodes.size(); ++i)
    if (_nodes[i].first == BusyNode) return true;
  return false;
}
