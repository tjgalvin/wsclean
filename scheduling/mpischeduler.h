#ifndef MPI_SCHEDULER_H
#define MPI_SCHEDULER_H

#include <mutex>
#include <thread>
#include <condition_variable>

#ifdef HAVE_MPI

#include "griddingtaskmanager.h"
#include "threadedscheduler.h"

class MPIScheduler final : public GriddingTaskManager {
 public:
  MPIScheduler(const class WSCleanSettings &settings);
  ~MPIScheduler();

  void Run(GriddingTask &task,
           std::function<void(GriddingResult &)> finishCallback) override;
  void Finish() override;

 private:
  enum NodeState { AvailableNode, BusyNode };

  void sendLoop();
  void receiveLoop();

  void node0gridder(GriddingTask task);

  int findAndSetNodeState(
      MPIScheduler::NodeState currentState,
      std::pair<MPIScheduler::NodeState, std::function<void(GriddingResult &)>>
          newState);
  bool anyReceiveTasks_NeedLock();

  bool _isRunning, _isFinishing, _isSendFinished;
  std::condition_variable _notify;
  std::mutex _mutex;
  std::thread _sendThread, _receiveThread, _workThread;
  aocommon::Lane<std::pair<GriddingTask, std::function<void(GriddingResult &)>>>
      _taskList;
  std::vector<std::pair<GriddingResult, std::function<void(GriddingResult &)>>>
      _readyList;
  std::vector<std::pair<NodeState, std::function<void(GriddingResult &)>>>
      _nodes;
};

#endif  // HAVE_MPI

#endif  // MPI_SCHEDULER_H
