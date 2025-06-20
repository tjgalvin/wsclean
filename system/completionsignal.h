#ifndef WSCLEAN_SYSTEM_COMPLETION_SIGNAL_H_
#define WSCLEAN_SYSTEM_COMPLETION_SIGNAL_H_

#include <mutex>

namespace wsclean {

/*
 * Single use exception safe utility class to allow multiple consumer threads to
 * wait until the producing/constructing thread signals completion. Must be used
 * only once, i.e. after the signal is sent the instance should be destroyed and
 * not used again.
 *
 * Example usage:
 * Reduce the number of active worker threads in an @ref TaskQueue and later
 * increasing them back to the original amount, without having to delete/create
 * threads or interrupt the @ref TaskQueue in any way.
 *
 * void Producer(TaskQueue& queue)
 * {
 *     ExecuteTasksWithAllThreadsAvailable(queue);
 *     // Reduce number of worker threads by 4
 *     CompletionSignal signal;
 *     for(int i=0;i<4;++i)
 *     {
 *         queue.Emplace([&]() {
 *           signal.LockUntilCompletion(block_excess_tasks);
 *         });
 *     }
 *     ExecuteTasksWithLessThreadsAvailable(queue);
 *     // Increase number of worker threads by 4
 *     signal.SignalCompletion();
 *     ExecuteMoreTasksWithAllThreadsAvailable(queue);
 *     // signal object should not be used again after this
 * }
 */
class CompletionSignal {
 public:
  CompletionSignal() : lock_(mutex_) {
    // Hold a lock on the mutex so that when tasks call WaitForCompletion()
    // they will not be able to lock the mutex and will enter a
    // paused/frozen/blocked state, where they will remain until
    // SignalCompletion() is called to free the lock on the mutex.
  }
  CompletionSignal(const CompletionSignal&) = delete;
  /*
   * Signal waiting threads to stop waiting.
   * Must be called by the same thread that created the instance of this class.
   * Must be called only once. Object should not be reused for additional @ref
   * WaitForCompletion() or @ref CompletionSignal() calls.
   */
  void SignalCompletion() {
    // Releasing the lock allows all threads that are paused on the mutex inside
    // WaitForCompletion to resume.
    lock_.unlock();
  }
  /*
   * Wait until completion is signaled by the owner thread calling @ref
   * SignalCompletion()
   */
  void WaitForCompletion() {
    // As the mutex is already locked during construction all worker threads
    // that call WaitForCompletion will be unable to lock and therefore pause
    // here until SignalCompletion() is called to free the lock
    std::unique_lock<std::mutex> lock(mutex_);
  }

 private:
  std::mutex mutex_;
  std::unique_lock<std::mutex> lock_;
};

}  // namespace wsclean

#endif  // WSCLEAN_SYSTEM_COMPLETION_SIGNAL_H_
