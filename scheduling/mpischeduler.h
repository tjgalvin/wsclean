#ifndef SCHEDULING_MPI_SCHEDULER_H_
#define SCHEDULING_MPI_SCHEDULER_H_

#ifdef HAVE_MPI

#include "griddingtaskmanager.h"
#include "threadedscheduler.h"

#include <aocommon/queue.h>
#include <aocommon/threadpool.h>
#include <aocommon/uvector.h>

#include <mutex>
#include <thread>
#include <condition_variable>

namespace wsclean {

class MPIScheduler final : public GriddingTaskManager {
 public:
  MPIScheduler(const Settings& settings);
  ~MPIScheduler() final {
    send_task_queue_.Finish();
    Finish();
    send_thread_pool_.FinishParallelExecution();
  }

  /**
   * Main Run function.
   * Either sends the task to another MPI node or runs it locally.
   */
  void Run(GriddingTask&& task,
           std::function<void(GriddingResult&)> finish_callback) final;

  void Finish() final;

  void Start(size_t n_writer_groups) final;

  std::unique_ptr<WriterLock> GetLock(size_t writer_group_index) final {
    // Since WSClean uses a static outputchannel-to-node mapping,
    // synchronisation of writes only needs to happen within a node.
    return local_scheduler_.GetLock(writer_group_index);
  }

 private:
  /**
   * Send a task to a worker node or run it on the master
   */
  void SendToNode(size_t node, std::shared_ptr<GriddingTask> task);

  /**
   * Wait until results are available and push these to the 'ready list'.
   * The loop ends when Finish() is called and all tasks are finished.
   * This function runs in a separate thread.
   */
  void ReceiveLoop();

  /**
   * Gets a node index for executing a (compound) task according
   * to the channel to index mapping.
   * Updates the number of tasks assigned to the node and stores
   * the callback function.
   * @return The index of the node executing the task.
   */
  size_t GetNode(const GriddingTask& task,
                 std::function<void(GriddingResult&)>&& callback);

  /**
   * If any results are available, call the callback functions and remove these
   * results from the ready list. This function
   * should be called by the main thread only, so that the user of the MPI
   * scheduler does not need to synchronize.
   *
   * This function is UNSYNCHRONIZED: the caller should
   * hold the mutex locked while calling it.
   */
  void ProcessReadyList_UNSYNCHRONIZED();

  /**
   * Return true if any tasks are still running on worker nodes.
   * Remember that the return value is independent of the state of the
   * main node: when the main node is gridding, it will nevertheless
   * return false if the other nodes are not running tasks.
   *
   * This function is UNSYNCHRONIZED: the caller should
   * hold the mutex locked while calling it.
   */
  bool AWorkerIsRunning_UNSYNCHRONIZED();

  void ProcessGriddingResult(size_t node, size_t body_size);
  /**
   * Stores 'result' in _readyList and updates the available slots of 'node'.
   */
  void StoreResult(GriddingResult&& result, int node);

  bool is_running_ = false;
  bool is_finishing_ = false;
  std::condition_variable notify_;
  std::mutex mutex_;
  std::thread receive_thread_;
  /** Stores results of ready tasks. */
  std::vector<GriddingResult> ready_list_;
  /** Stores callbacks, indexed by task id. */
  std::map<size_t, std::function<void(GriddingResult&)>> callbacks_;

  /**
   * Available execution room for tasks for each node.
   * This is intially set to `room_per_node_` and is reduced/increased as tasks
   * are sent/completed.
   * A compound task; with multiple facets, counts as multiple tasks,
   * with the quantity determined by `task.num_parallel_gridders_`.
   * When room is <= 0 no further tasks are sent until room becomes > 0.
   * The room can become negative if the scheduler has sent more tasks to the
   * node than `room_per_node_` suggests it is capable of executing in parallel.
   * This situation occurs if a compound task is being sent pre-emptively. When
   * the node finishes a task, it can then immediately start with the
   * pre-emptively sent task instead of waiting for a new task. Or even start to
   * process the second in an overlapped manner as the first task starts to
   * complete.
   */
  std::vector<int> available_room_;
  int room_per_node_;
  /** Stores task size, indexed by task id. */
  std::map<size_t, size_t> task_size_;

  /**
   * Prevent more than one concurrent send per node at a time.
   */
  aocommon::UVector<bool> node_is_sending_;

  /**
   * The lower-level local scheduler on an MPI node.
   * Using the threaded scheduler ensures that gridding uses a separate thread.
   */
  ThreadedScheduler local_scheduler_;

  /**
   * Internal task queue to enable sending multiple MPI tasks in parallel.
   */
  aocommon::ThreadPool send_thread_pool_;
  aocommon::TaskQueue<std::function<void()>> send_task_queue_;
};

}  // namespace wsclean

#endif  // HAVE_MPI

#endif  // MPI_SCHEDULER_H
