#ifndef WSCLEAN_SYSTEM_COMPLETION_SIGNAL_H_
#define WSCLEAN_SYSTEM_COMPLETION_SIGNAL_H_

#include <mutex>

namespace wsclean {

/*
 * Single use exception safe utility class to allow multiple consumer threads to
 * wait until completion is signalled. Can be signalled from a different thread
 * than the one in which it was created.
 * Must be used only once, i.e. after the signal is sent the instance should be
 * destroyed and not used again.
 *
 * Example usage:
 * Reduce the number of active worker threads in an @ref TaskQueue and later
 * increasing them back to the original amount, without having to delete/create
 * threads or interrupt the @ref TaskQueue in any way.
 *
 * void Producer(TaskQueue& queue)
 * {
 *     ExecuteTasks(queue); // N tasks run in parallel
 *     // Reduce number of worker threads by 4
 *     CompletionSignal signal;
 *     for(int i=0;i<4;++i)
 *     {
 *         queue.Emplace([&]() {
 *           signal.WaitForCompletion();
 *         });
 *     }
 *     ExecuteTasks(queue); // N-4 tasks run in parallel
 *     // Increase number of worker threads by 4
 *     signal.SignalCompletion();
 *     // signal object should not be used again after this
 *     ExecuteTasks(queue); // N tasks run in parallel
 * }
 */
class CompletionSignal {
 public:
  /*
   * Initialise internal signal state such that tasks calling
   * WaitForCompletion() will be paused/frozen/blocked.
   * Call @ref SignalCompletion() to signal completion and release the tasks
   * from their waiting state.
   */
  CompletionSignal() = default;
  CompletionSignal(const CompletionSignal&) = delete;
  CompletionSignal& operator=(const CompletionSignal&) = delete;

  /*
   * Signal paused/frozen/blocked tasks to stop waiting.
   * Can be called from a different thread than the one that created the
   * instance of this class. Must be called only once. Object should not be
   * reused for additional @ref WaitForCompletion() or @ref CompletionSignal()
   * calls.
   */
  void SignalCompletion() {
    signal_ = true;
    signal_.notify_all();
  }
  /*
   * Wait until completion is signaled by calling @ref SignalCompletion()
   * Note that the call does not need to come from the same thread that
   */
  void WaitForCompletion() {
    // `signal_` is set to false during construction.
    // All worker threads that call WaitForCompletion() will wait here until
    // SignalCompletion() is called to set signal to true.
    signal_.wait(false);
  }

 private:
  std::atomic<bool> signal_ = false;
};

}  // namespace wsclean

#endif  // WSCLEAN_SYSTEM_COMPLETION_SIGNAL_H_
