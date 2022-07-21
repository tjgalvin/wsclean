#ifndef THREADED_SCHEDULER_H
#define THREADED_SCHEDULER_H

#include "griddingtaskmanager.h"

#include "../structures/resources.h"

#include <mutex>
#include <thread>

class ThreadedScheduler final : public GriddingTaskManager {
 public:
  ThreadedScheduler(const class Settings& settings);
  ~ThreadedScheduler();

  void Run(GriddingTask&& task,
           std::function<void(GriddingResult&)> finishCallback) override;
  void Finish() override;

  void Start(size_t nWriterGroups) override;

  LockGuard GetLock(size_t writerGroupIndex) override;

 private:
  class ThreadedWriterLock final : public WriterLock {
   public:
    void lock() override { _mutex.lock(); }
    void unlock() override { _mutex.unlock(); }

   private:
    std::mutex _mutex;
  };

  void processQueue();

  std::mutex mutex_;
  std::vector<std::thread> thread_list_;
  aocommon::Lane<std::pair<GriddingTask, std::function<void(GriddingResult&)>>>
      task_list_;
  std::vector<std::pair<GriddingResult, std::function<void(GriddingResult&)>>>
      ready_list_;
  std::vector<ThreadedWriterLock> writer_group_locks_;

  const Resources resources_per_task_;
};

#endif
