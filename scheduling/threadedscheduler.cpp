#include "threadedscheduler.h"
#include "../gridding/msgridderbase.h"

#include "../main/settings.h"

ThreadedScheduler::ThreadedScheduler(const class Settings& settings)
    : GriddingTaskManager(settings), _taskList(settings.parallelGridding) {}

ThreadedScheduler::~ThreadedScheduler() {
  if (!_threadList.empty()) Finish();
}

void ThreadedScheduler::Run(
    GriddingTask&& task, std::function<void(GriddingResult&)> finishCallback) {
  // Start an extra thread if not maxed out already
  if (_threadList.size() < _settings.parallelGridding)
    _threadList.emplace_back(&ThreadedScheduler::processQueue, this);
  else
    _taskList
        .wait_for_empty();  // if all threads are busy, block until one
                            // available (in order not to stack too many tasks)

  std::lock_guard<std::mutex> lock(_mutex);
  while (!_readyList.empty()) {
    // Call callbacks for any finished tasks
    _readyList.back().second(_readyList.back().first);
    _readyList.pop_back();
  }

  _taskList.emplace(std::move(task), std::move(finishCallback));
}

void ThreadedScheduler::processQueue() {
  std::pair<GriddingTask, std::function<void(GriddingResult&)>> taskPair;
  while (_taskList.read(taskPair)) {
    std::unique_ptr<MSGridderBase> gridder(makeGridder());
    GriddingResult result = runDirect(std::move(taskPair.first), *gridder);

    std::lock_guard<std::mutex> lock(_mutex);
    _readyList.emplace_back(std::move(result), taskPair.second);
  }
}

void ThreadedScheduler::Start(size_t nWriterGroups) {
  GriddingTaskManager::Start(nWriterGroups);

  if (_writerGroupLocks.size() < nWriterGroups)
    _writerGroupLocks = std::vector<ThreadedWriterLock>(nWriterGroups);
}

WriterLockManager::LockGuard ThreadedScheduler::GetLock(
    size_t writerGroupIndex) {
  return LockGuard(_writerGroupLocks[writerGroupIndex]);
}

void ThreadedScheduler::Finish() {
  _taskList.write_end();
  for (std::thread& t : _threadList) t.join();
  _threadList.clear();
  _taskList.clear();
  while (!_readyList.empty()) {
    // Call callbacks for any finished tasks
    _readyList.back().second(_readyList.back().first);
    _readyList.pop_back();
  }
}
