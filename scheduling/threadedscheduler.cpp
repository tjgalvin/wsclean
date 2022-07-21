#include "threadedscheduler.h"
#include "../gridding/msgridderbase.h"

#include "../main/settings.h"

ThreadedScheduler::ThreadedScheduler(const Settings& settings)
    : GriddingTaskManager(settings),
      task_list_(settings.parallelGridding),
      resources_per_task_(GetResources().GetPart(settings.parallelGridding)) {}

ThreadedScheduler::~ThreadedScheduler() {
  if (!thread_list_.empty()) Finish();
}

void ThreadedScheduler::Run(
    GriddingTask&& task, std::function<void(GriddingResult&)> finishCallback) {
  // Start an extra thread if not maxed out already
  if (thread_list_.size() < _settings.parallelGridding)
    thread_list_.emplace_back(&ThreadedScheduler::processQueue, this);
  else
    task_list_
        .wait_for_empty();  // if all threads are busy, block until one
                            // available (in order not to stack too many tasks)

  std::lock_guard<std::mutex> lock(mutex_);
  while (!ready_list_.empty()) {
    // Call callbacks for any finished tasks
    ready_list_.back().second(ready_list_.back().first);
    ready_list_.pop_back();
  }

  task_list_.emplace(std::move(task), std::move(finishCallback));
}

void ThreadedScheduler::processQueue() {
  std::pair<GriddingTask, std::function<void(GriddingResult&)>> taskPair;
  while (task_list_.read(taskPair)) {
    std::unique_ptr<MSGridderBase> gridder(makeGridder(resources_per_task_));
    GriddingResult result = runDirect(std::move(taskPair.first), *gridder);

    std::lock_guard<std::mutex> lock(mutex_);
    ready_list_.emplace_back(std::move(result), taskPair.second);
  }
}

void ThreadedScheduler::Start(size_t nWriterGroups) {
  GriddingTaskManager::Start(nWriterGroups);

  if (writer_group_locks_.size() < nWriterGroups)
    writer_group_locks_ = std::vector<ThreadedWriterLock>(nWriterGroups);
}

WriterLockManager::LockGuard ThreadedScheduler::GetLock(
    size_t writerGroupIndex) {
  return LockGuard(writer_group_locks_[writerGroupIndex]);
}

void ThreadedScheduler::Finish() {
  task_list_.write_end();
  for (std::thread& t : thread_list_) t.join();
  thread_list_.clear();
  task_list_.clear();
  while (!ready_list_.empty()) {
    // Call callbacks for any finished tasks
    ready_list_.back().second(ready_list_.back().first);
    ready_list_.pop_back();
  }
}
