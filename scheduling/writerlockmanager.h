#ifndef WRITER_LOCK_MANAGER_H
#define WRITER_LOCK_MANAGER_H

#include <mutex>

/**
 * @brief Abstract interface class providing access to the (writer)
 * locking mechanism of the GriddingTaskManager or childs thereof.
 */
class WriterLockManager {
 protected:
  class WriterLock {
   public:
    virtual ~WriterLock() = default;
    virtual void lock() = 0;
    virtual void unlock() = 0;
  };

 public:
  virtual ~WriterLockManager(){};

#if __cplusplus >= 201703L
  using LockGuard = std::lock_guard<WriterLock>;
#else
  // In GetLock(), 'return LockGuard(...);' is only possible if LockGuard is
  // movable.
  using LockGuard = std::unique_lock<WriterLock>;
#endif

  /**
   * @brief Return a writer lock guard for the given @p writerGroupIndex
   */
  virtual LockGuard GetLock(size_t writerGroupIndex) = 0;
};
#endif