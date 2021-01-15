#ifndef THREADED_SCHEDULER_H
#define THREADED_SCHEDULER_H

#include "griddingtaskmanager.h"

#include <mutex>
#include <thread>

class ThreadedScheduler final : public GriddingTaskManager {
 public:
  ThreadedScheduler(const class Settings& settings);
  ~ThreadedScheduler();

  void Run(GriddingTask&& task,
           std::function<void(GriddingResult&)> finishCallback) override;
  void Finish() override;

 private:
  void processQueue();

  std::mutex _mutex;
  std::vector<std::thread> _threadList;
  aocommon::Lane<std::pair<GriddingTask, std::function<void(GriddingResult&)>>>
      _taskList;
  std::vector<std::pair<GriddingResult, std::function<void(GriddingResult&)>>>
      _readyList;
};

#endif
