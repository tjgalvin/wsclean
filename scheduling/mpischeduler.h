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
  MPIScheduler(const class Settings &settings);
  ~MPIScheduler();

  void Run(GriddingTask &&task,
           std::function<void(GriddingResult &)> finishCallback) override;

  void Finish() override;

 private:
  enum NodeState { AvailableNode, BusyNode };

  /**
   * Send a task to a worker node or run it on the master
   * If all nodes are busy, the call will block until a node is
   * available.
   */
  void send(GriddingTask &&task,
            const std::function<void(GriddingResult &)> &callback);

  /**
   * Wait until results are available and push these to the 'ready list'.
   * The loop ends when Finish() is called and all tasks are finished.
   * This function runs in a separate thread.
   */
  void receiveLoop();

  /**
   * Directly run the given task. This is a blocking call.
   */
  void runTaskOnNode0(GriddingTask &&task);

  /**
   * This "atomically" finds a node with a certain state and assigns a new value
   * to it. The return value is the index of the node that matched the state.
   * Searching is done from the last node to the first, so that the master node
   * is selected last.
   */
  int findAndSetNodeState(
      MPIScheduler::NodeState currentState,
      std::pair<MPIScheduler::NodeState, std::function<void(GriddingResult &)>>
          newState);

  /**
   * If any results are available, call the callback functions and remove these
   * results from the ready list. This function
   * should be called by the main thread only, so that the user of the MPI
   * scheduler does not need to synchronize.
   *
   * This function is UNSYNCHRONIZED: the caller should
   * hold the mutex locked while calling it.
   */
  void processReadyList_UNSYNCHRONIZED();

  /**
   * Return true if any tasks are still running.
   * Remember that the return value is independent of the state of the
   * master node: when the master node is gridding, it will nevertheless
   * return false if the other nodes are not running tasks.
   *
   * This function is UNSYNCHRONIZED: the caller should
   * hold the mutex locked while calling it.
   */
  bool receiveTasksAreRunning_UNSYNCHRONIZED();

  const bool _masterDoesWork;
  bool _isRunning, _isFinishing;
  std::condition_variable _notify;
  std::mutex _mutex;
  std::thread _receiveThread, _workThread;
  std::vector<std::pair<GriddingResult, std::function<void(GriddingResult &)>>>
      _readyList;
  std::vector<std::pair<NodeState, std::function<void(GriddingResult &)>>>
      _nodes;
};

#endif  // HAVE_MPI

#endif  // MPI_SCHEDULER_H
