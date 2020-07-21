#ifndef WGRIDDER_THREADING_H
#define WGRIDDER_THREADING_H

#include <cstdlib>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <queue>
#include <atomic>
#include <functional>
#include <vector>

#ifdef POCKETFFT_PTHREADS
#include <pthread.h>
#endif

namespace threading {

namespace detail {

using namespace std;

thread_local size_t thread_id = 0;
thread_local size_t num_threads = 1;
static const size_t max_threads = max(1u, thread::hardware_concurrency());

class latch {
  atomic<size_t> num_left_;
  mutex mut_;
  condition_variable completed_;
  using lock_t = unique_lock<mutex>;

 public:
  latch(size_t n) : num_left_(n) {}

  void count_down() {
    lock_t lock(mut_);
    if (--num_left_) return;
    completed_.notify_all();
  }

  void wait() {
    lock_t lock(mut_);
    completed_.wait(lock, [this] { return is_ready(); });
  }
  bool is_ready() { return num_left_ == 0; }
};

template <typename T>
class concurrent_queue {
  queue<T> q_;
  mutex mut_;
  condition_variable item_added_;
  bool shutdown_;
  using lock_t = unique_lock<mutex>;

 public:
  concurrent_queue() : shutdown_(false) {}

  void push(T val) {
    {
      lock_t lock(mut_);
      if (shutdown_) throw runtime_error("Item added to queue after shutdown");
      q_.push(move(val));
    }
    item_added_.notify_one();
  }

  bool pop(T &val) {
    lock_t lock(mut_);
    item_added_.wait(lock, [this] { return (!q_.empty() || shutdown_); });
    if (q_.empty()) return false;  // We are shutting down

    val = std::move(q_.front());
    q_.pop();
    return true;
  }

  void shutdown() {
    {
      lock_t lock(mut_);
      shutdown_ = true;
    }
    item_added_.notify_all();
  }

  void restart() { shutdown_ = false; }
};

class thread_pool {
  concurrent_queue<function<void()>> work_queue_;
  vector<thread> threads_;

  void worker_main() {
    function<void()> work;
    while (work_queue_.pop(work)) work();
  }

  void create_threads() {
    size_t nthreads = threads_.size();
    for (size_t i = 0; i < nthreads; ++i) {
      try {
        threads_[i] = thread([this] { worker_main(); });
      } catch (...) {
        shutdown();
        throw;
      }
    }
  }

 public:
  explicit thread_pool(size_t nthreads) : threads_(nthreads) {
    create_threads();
  }

  thread_pool() : thread_pool(max_threads) {}

  ~thread_pool() { shutdown(); }

  void submit(function<void()> work) { work_queue_.push(move(work)); }

  void shutdown() {
    work_queue_.shutdown();
    for (auto &thread : threads_)
      if (thread.joinable()) thread.join();
  }

  void restart() {
    work_queue_.restart();
    create_threads();
  }
};

thread_pool &get_pool() {
  static thread_pool pool;
#ifdef POCKETFFT_PTHREADS
  static once_flag f;
  call_once(f, [] {
    pthread_atfork(
        +[] { get_pool().shutdown(); },  // prepare
        +[] { get_pool().restart(); },   // parent
        +[] { get_pool().restart(); }    // child
    );
  });
#endif

  return pool;
}

/** Map a function f over nthreads */
template <typename Func>
void thread_map(size_t nthreads, Func f) {
  if (nthreads == 0) nthreads = max_threads;

  if (nthreads == 1) {
    f();
    return;
  }

  auto &pool = get_pool();
  latch counter(nthreads);
  exception_ptr ex;
  mutex ex_mut;
  for (size_t i = 0; i < nthreads; ++i) {
    pool.submit([&f, &counter, &ex, &ex_mut, i, nthreads] {
      thread_id = i;
      num_threads = nthreads;
      try {
        f();
      } catch (...) {
        lock_guard<mutex> lock(ex_mut);
        ex = current_exception();
      }
      counter.count_down();
    });
  }
  counter.wait();
  if (ex) rethrow_exception(ex);
}

template <typename It>
class Scheduler {
 private:
  using difftype = decltype(declval<It>() - declval<It>());
  size_t nthreads_;
  mutex mut_;
  It beg_, end_;
  It cur_;
  difftype chunksize_;
  double fact_max_;
  vector<difftype> nextstart;
  typedef enum { SINGLE, STATIC, DYNAMIC } SchedMode;
  SchedMode mode;
  bool single_done;

 public:
  size_t nthreads() const { return nthreads_; }
  mutex &mut() { return mut_; }

  template <typename Func>
  void execSingle(It beg, It end, Func f) {
    mode = SINGLE;
    single_done = false;
    beg_ = beg;
    end_ = end;
    f();
  }
  template <typename Func>
  void execStatic(It beg, It end, size_t nthreads, difftype chunksize, Func f) {
    mode = STATIC;
    nthreads_ = (nthreads == 0) ? max_threads : nthreads;
    beg_ = beg;
    end_ = end;
    chunksize_ =
        (chunksize < 1) ? (end_ + nthreads_ - 1 - beg_) / nthreads_ : chunksize;
    if (chunksize_ >= end_ - beg_) return execSingle(beg, end, move(f));
    nextstart.resize(nthreads_);
    for (size_t i = 0; i < nextstart.size(); ++i)
      nextstart[i] = beg_ + i * chunksize_;
    thread_map(nthreads_, move(f));
  }
  template <typename Func>
  void execDynamic(It beg, It end, size_t nthreads, difftype chunksize_min,
                   double fact_max, Func f) {
    mode = DYNAMIC;
    nthreads_ = (nthreads == 0) ? max_threads : nthreads;
    beg_ = beg;
    end_ = end;
    chunksize_ = (chunksize_min < 1) ? 1 : chunksize_min;
    if (chunksize_ * nthreads_ >= end_ - beg_)
      return execStatic(beg, end, nthreads, 0, move(f));
    fact_max_ = fact_max;
    cur_ = beg_;
    thread_map(nthreads_, move(f));
  }
  bool getNext(It &lo, It &hi) {
    switch (mode) {
      case SINGLE: {
        if (single_done) return false;
        single_done = true;
        lo = beg_;
        hi = end_;
        return true;
      }
      case STATIC: {
        if (nextstart[thread_id] >= end_) return false;
        lo = nextstart[thread_id];
        hi = min(lo + chunksize_, end_);
        nextstart[thread_id] += nthreads_ * chunksize_;
        return true;
      }
      case DYNAMIC: {
        unique_lock<mutex> lck(mut_);
        if (cur_ >= end_) return false;
        auto rem = end_ - cur_;
        difftype tmp = difftype((fact_max_ * rem) / nthreads_);
        auto sz = min(rem, max(chunksize_, tmp));
        lo = cur_;
        cur_ += sz;
        hi = cur_;
        return true;
      }
    }
    return false;
  }
};

}  // end of namespace detail

using detail::Scheduler;

}  // end of namespace threading

#endif
