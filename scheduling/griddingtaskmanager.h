#ifndef GRIDDING_OPERATION_H
#define GRIDDING_OPERATION_H

#include "griddingtask.h"
#include "griddingresult.h"
#include "writerlockmanager.h"

#include "../structures/imageweights.h"
#include "../structures/observationinfo.h"
#include "../structures/msselection.h"

#include <aocommon/polarization.h>

#include "../msproviders/msdatadescription.h"

#include <aocommon/lane.h>

#include <cstring>
#include <functional>
#include <vector>

class MSGridderBase;
class Resources;
class Settings;

class GriddingTaskManager : protected WriterLockManager {
 public:
  virtual ~GriddingTaskManager();

  /**
   * Initialize writer groups. Call this function before scheduling Predict
   * tasks in order to initialize the writer locks.
   *
   * @param nWriterGroups The number of writer groups.
   */
  virtual void Start([[maybe_unused]] size_t nWriterGroups) {}

  LockGuard GetLock([[maybe_unused]] size_t writerGroupIndex) override {
    static DummyWriterLock dummy;
    return LockGuard(dummy);
  }

  /**
   * Add the given task to the queue of tasks to be run. After finishing
   * the task, the callback is called with the results. The callback will
   * always run in the thread of the caller.
   * Depending on the type of gridding task manager, this call might block.
   *
   * This implementation runs the task directly and blocks until done.
   */
  virtual void Run(GriddingTask&& task,
                   std::function<void(GriddingResult&)> finishCallback);

  /**
   * Block until all tasks have finished.
   */
  virtual void Finish(){};

  /**
   * Run the given task. This variant of Run() does not call a
   * callback function when finished. A gridder is created for
   * the duration of the call.
   */
  GriddingResult RunDirect(GriddingTask&& task);

  /**
   * Make the gridding task manager according to the settings.
   */
  static std::unique_ptr<GriddingTaskManager> Make(
      const class Settings& settings);

 protected:
  GriddingTaskManager(const Settings& settings);
  Resources GetResources() const;

  const Settings& _settings;

  std::unique_ptr<MSGridderBase> makeGridder(const Resources& resources) const;

  /**
   * Run the provided task with the specified gridder.
   */
  GriddingResult runDirect(GriddingTask&& task, MSGridderBase& gridder);

 private:
  class DummyWriterLock final : public WriterLock {
   public:
    void lock() override {}
    void unlock() override {}
  };

  std::unique_ptr<MSGridderBase> constructGridder(
      const Resources& resources) const;
};

#endif
