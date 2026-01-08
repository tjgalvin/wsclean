#ifndef SCHEDULING_THREADED_SCHEDULER_H_
#define SCHEDULING_THREADED_SCHEDULER_H_

#include <atomic>
#include <memory>
#include <mutex>
#include <thread>

#include <aocommon/taskqueue.h>

#include "griddingtaskmanager.h"

#include "../structures/resources.h"

namespace wsclean {

class ThreadedScheduler final : public GriddingTaskManager {
 public:
  ThreadedScheduler(const class Settings& settings);
  ~ThreadedScheduler();

  void Run(GriddingTask&& task,
           std::function<void(GriddingResult&)> finishCallback) override;
  void Finish() override;

  void Start(size_t nWriterGroups) override;

  std::unique_ptr<WriterLock> GetLock(size_t writer_group_index) override;

 private:
  class ThreadedWriterLock final : public WriterLock {
   public:
    explicit ThreadedWriterLock(ThreadedScheduler& scheduler,
                                size_t writer_group_index)
        : scheduler_{scheduler}, writer_group_index_{writer_group_index} {
      scheduler.writer_group_locks_[writer_group_index].lock();
    }

    ~ThreadedWriterLock() override {
      scheduler_.writer_group_locks_[writer_group_index_].unlock();
    }

   private:
    ThreadedScheduler& scheduler_;
    size_t writer_group_index_;
  };

  friend class ThreadedWriterLock;

  /// Contains all data for a single task.
  struct TaskData {
    GriddingTask task;
    GriddingResult result;
    std::mutex result_mutex;
    std::atomic<std::size_t> finished_facet_count;
    std::function<void(GriddingResult&)> callback;
  };

  /// Process all tasks from the queue as they become available.
  /// Calls @ref ProcessTask() for each task.
  void ProcessQueue();
  /// Perform processing for an individual task.
  void ProcessTask(
      size_t task_id, const std::vector<size_t>& facet_indexes,
      TaskData& task_data, std::shared_ptr<CompletionSignal>& lock_excess_tasks,
      std::unique_ptr<MSGridderManagerScheduler> scheduler = nullptr);
  /// Attempt to immediately pop a task from the queue and call @ref
  /// ProcessTask() on it using the passed in resources of a currently running
  /// task. This allows for the resources of tasks that are near completion to
  /// be reused in an overlapped way by the next task to improve throughput.
  bool TryStealTask(
      std::shared_ptr<CompletionSignal>& parent_task_completion_signal,
      std::shared_ptr<CompletionSignal>& signal_stolen_task_completion,
      std::unique_ptr<MSGridderManagerScheduler>& parent_scheduler,
      size_t n_work_size);
  /// Execute post-processing callback on any tasks that have completed and
  /// queued as ready.
  void ProcessReadyList();

  /// The first value of the pair is the unique id for the tasks, for looking up
  /// the task and its callback into task_data_map_.
  /// The second value contains the indices of the facet(s) for the sub-task(s).
  using TaskQueueType =
      aocommon::TaskQueue<std::pair<size_t, std::vector<size_t>>>;

  /// Protects task_data_map_, latest_exception_ and ready_list_.
  std::mutex mutex_;
  std::vector<std::thread> thread_list_;
  /// Stores the data for the pending tasks.
  /// The map is indexed by the unique_id of the task.
  std::map<std::size_t, TaskData> task_data_map_;
  /// FIFO queue for the tasks.
  /// Each sub-task, which processes a single facet, has an entry.
  TaskQueueType task_queue_;
  /// Stores the unique_id for tasks that are complete, so the main thread
  /// can execute the callback.
  std::vector<size_t> ready_list_;
  /// Stores the latest exception that occurred while running a task in a
  /// thread, so the scheduler can rethrow it in the main thread.
  std::exception_ptr latest_exception_;
  std::vector<std::mutex> writer_group_locks_;

  const Resources resources_per_task_;
};

}  // namespace wsclean

#endif
