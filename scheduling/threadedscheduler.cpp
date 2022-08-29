#include "threadedscheduler.h"
#include "../gridding/msgridderbase.h"

#include "../main/settings.h"

#include <aocommon/logger.h>

#include <string>

ThreadedScheduler::ThreadedScheduler(const Settings& settings)
    : GriddingTaskManager(settings),
      task_list_(settings.parallelGridding),
      resources_per_task_(GetResources().GetPart(settings.parallelGridding)) {}

ThreadedScheduler::~ThreadedScheduler() {
  if (!thread_list_.empty()) {
    try {
      Finish();
    } catch (std::exception& e) {
      // Normally, the user of the ThreadedScheduler calls Finish(), which
      // empties the thread_list_. If thread_list_ is non-empty in this
      // destructor, the ThreadedScheduler is destroyed because some another
      // exception occurred. We are in a destructor, so all that can be done is
      // report the error.
      using namespace std::string_literals;
      aocommon::Logger::Error
          << "Exception caught during destruction of ThreadedScheduler:\n"s +
                 e.what() + '\n';
    }
  }
}

void ThreadedScheduler::Run(
    GriddingTask&& task, std::function<void(GriddingResult&)> finishCallback) {
  // Start an extra thread if not maxed out already
  if (thread_list_.size() < _settings.parallelGridding)
    thread_list_.emplace_back(&ThreadedScheduler::ProcessQueue, this);
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
  CheckExceptions();
}

void ThreadedScheduler::ProcessQueue() {
  std::pair<GriddingTask, std::function<void(GriddingResult&)>> taskPair;
  while (task_list_.read(taskPair)) {
    try {
      std::unique_ptr<MSGridderBase> gridder(makeGridder(resources_per_task_));
      GriddingResult result = runDirect(std::move(taskPair.first), *gridder);

      std::lock_guard<std::mutex> lock(mutex_);
      ready_list_.emplace_back(std::move(result), taskPair.second);
    } catch (std::exception&) {
      std::lock_guard<std::mutex> lock(mutex_);
      latest_exception_ = std::current_exception();
    }
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
  // All threads have joined, so no lock required
  while (!ready_list_.empty()) {
    // Call callbacks for any finished tasks
    ready_list_.back().second(ready_list_.back().first);
    ready_list_.pop_back();
  }
  CheckExceptions();
}

void ThreadedScheduler::CheckExceptions() {
  if (latest_exception_) {
    std::exception_ptr to_throw = std::move(latest_exception_);
    latest_exception_ = std::exception_ptr();
    std::rethrow_exception(to_throw);
  }
}
