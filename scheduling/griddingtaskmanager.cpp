#include "griddingtaskmanager.h"

#include <numeric>
#include <mutex>
#include <vector>

#include "griddingtask.h"
#include "griddingresult.h"
#include "mpischeduler.h"
#include "threadedscheduler.h"

#include "../gridding/h5solutiondata.h"
#include "../gridding/msgriddermanager.h"
#include "../main/settings.h"
#include "../structures/resources.h"

namespace wsclean {

GriddingTaskManager::GriddingTaskManager(const Settings& settings)
    : settings_(settings),
      solution_data_(settings),
      writer_lock_manager_(this) {}

GriddingTaskManager::~GriddingTaskManager() {}

std::unique_ptr<GriddingTaskManager> GriddingTaskManager::Make(
    const Settings& settings) {
  if (settings.UseMpi()) {
#ifdef HAVE_MPI
    return std::make_unique<MPIScheduler>(settings);
#else
    throw std::runtime_error("MPI not available");
#endif
  } else if (settings.parallelGridding > 1) {
    return std::make_unique<ThreadedScheduler>(settings);
  } else {
    return std::make_unique<GriddingTaskManager>(settings);
  }
}

Resources GriddingTaskManager::GetResources() const {
  return Resources(
      settings_.threadCount,
      GetAvailableMemory(settings_.memFraction, settings_.absMemLimit));
}

void GriddingTaskManager::Run(
    GriddingTask&& task, std::function<void(GriddingResult&)> finishCallback) {
  std::vector<size_t> facet_indices(task.facets.size());
  std::iota(facet_indices.begin(), facet_indices.end(), 0);

  GriddingResult result;
  result.facets.resize(task.facets.size());
  std::mutex result_mutex;

  RunDirect(task, facet_indices, GetResources(), result, result_mutex);

  finishCallback(result);
}

void GriddingTaskManager::RunDirect(GriddingTask& task,
                                    const std::vector<size_t>& facet_indices,
                                    const Resources& resources,
                                    GriddingResult& result,
                                    std::mutex& result_mutex) {
  assert(!facet_indices.empty());
  assert(result.facets.size() == task.facets.size());

  // Wait tasks occupy a thread from the pool by waiting on
  // lock_excess_scheduler_tasks_ which will freeze the tasks thread until
  // lock_excess_scheduler_tasks_.SignalCompletion() is called.
  if (task.operation == GriddingTask::Wait) {
    task.lock_excess_scheduler_tasks_->WaitForCompletion();
    return;
  }

  MSGridderManager manager(settings_, solution_data_);
  manager.InitializeMS(task);
  manager.InitializeGridders(task, facet_indices, resources, result.facets,
                             writer_lock_manager_);
  if (task.operation == GriddingTask::Invert) {
    if (settings_.shared_facet_reads) {
      manager.SortFacetTasks();
      manager.BatchInvert(task.num_parallel_gridders_);
    } else {
      manager.Invert();
    }
  } else {
    manager.Predict();
  }
  const bool store_common_info = (facet_indices.front() == 0);
  if (store_common_info) {
    result.unique_id = task.unique_id;
  }
  manager.ProcessResults(result_mutex, result, store_common_info);

  if (GetSettings().shared_facet_reads) {
    // Allow wait tasks to resume
    if (task.lock_excess_scheduler_tasks_) {
      task.lock_excess_scheduler_tasks_->SignalCompletion();
      task.lock_excess_scheduler_tasks_.reset();
    }
  }
}

}  // namespace wsclean
