#ifndef WSCLEAN_THREAD_POOL_H
#define WSCLEAN_THREAD_POOL_H

#include "system.h"

#include <functional>
#include <thread>
#include <mutex>
#include <condition_variable>

#include <queue>

class ThreadPool {
 public:
  ThreadPool() : _activeThreads(0), _finish(false) {
    init_threads(System::ProcessorCount());
  }

  ~ThreadPool() {
    std::unique_lock<std::mutex> lock(_mutex);
    _finish = true;
    _change.notify_all();
    lock.unlock();
    for (size_t i = 0; i != _threads.size(); ++i) {
      _threads[i].join();
    }
  }

  template <typename func>
  void queue(func f) {
    std::lock_guard<std::mutex> lock(_mutex);
    _queuedTasks.push(f);
    _change.notify_all();
  }

  void wait_for_all_tasks() {
    std::unique_lock<std::mutex> lock(_mutex);
    while (_activeThreads != 0 || !_queuedTasks.empty()) _change.wait(lock);
  }

  size_t size() const { return _threads.size(); }

 private:
  void init_threads(size_t n) {
    _threads.resize(n);
    for (size_t i = 0; i != n; ++i) {
      _threads[i] =
          std::thread(std::mem_fn(&ThreadPool::thread_function), this);
    }
  }

  void thread_function() {
    std::unique_lock<std::mutex> lock(_mutex);
    while (!_finish) {
      while (_queuedTasks.empty() && !_finish) _change.wait(lock);

      if (!_finish) {
        ++_activeThreads;
        auto f = _queuedTasks.front();
        _queuedTasks.pop();
        lock.unlock();
        f();
        lock.lock();
        --_activeThreads;

        _change.notify_all();
      }
    }
  }

  std::vector<std::thread> _threads;
  std::size_t _activeThreads;
  std::queue<std::function<void()>> _queuedTasks;
  bool _finish;

  std::mutex _mutex;
  std::condition_variable _change;
};

#endif
