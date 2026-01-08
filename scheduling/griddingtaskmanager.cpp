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

GriddingTaskManager::~GriddingTaskManager() = default;

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

  RunDirect(task, facet_indices, GetResources(), result, result_mutex, {});

  finishCallback(result);
}

void GriddingTaskManager::RunDirect(
    GriddingTask& task, const std::vector<size_t>& facet_indices,
    const Resources& resources, GriddingResult& result,
    std::mutex& result_mutex,
    std::function<void(std::unique_ptr<MSGridderManagerScheduler>&)>
        signal_last_work_has_started,
    std::unique_ptr<MSGridderManagerScheduler> scheduler) {
  assert(!facet_indices.empty());
  assert(result.facets.size() == task.facets.size());
  assert(!task.msList.empty());

  bool batch_task =
      task.operation == GriddingTask::Invert && settings_.shared_facet_reads;
  batch_task |=
      task.operation == GriddingTask::Predict && settings_.shared_facet_writes;
  const size_t n_threads = resources.NCpus();

  // Select which scheduler to use.
  // If we have been explicitely passed one then use that.
  // Alternatively try to re-use one from the cache.
  // As a last resort allocate a new one.
  std::unique_ptr<MSGridderManagerScheduler> selected_scheduler = nullptr;
  if (scheduler) {
    selected_scheduler = std::move(scheduler);
  } else if (batch_task) {
    {
      std::lock_guard<std::mutex> lock(scheduler_creation_mutex_);
      if (!scheduler_cache_[n_threads].empty()) {
        selected_scheduler = std::move(scheduler_cache_[n_threads].front());
        scheduler_cache_[n_threads].pop_front();
      }
    }
    if (!selected_scheduler) {
      selected_scheduler =
          std::make_unique<MSGridderManagerScheduler>(n_threads);
    }
  }
  MSGridderManager manager(settings_, solution_data_, selected_scheduler.get());
  manager.InitializeMS(task);
  manager.InitializeGridders(task, facet_indices, resources, result.facets,
                             writer_lock_manager_);
  if (task.operation == GriddingTask::Invert) {
    if (settings_.shared_facet_reads) {
      manager.BatchInvert([&]() {
        // NB! The signal can take ownership of the scheduler.
        std::lock_guard<std::mutex> lock(scheduler_creation_mutex_);
        signal_last_work_has_started(selected_scheduler);
      });
    } else {
      manager.Invert();
    }
  } else {
    if (settings_.shared_facet_writes) {
      manager.BatchPredict([&]() {
        // NB! The signal can take ownership of the scheduler.
        std::lock_guard<std::mutex> lock(scheduler_creation_mutex_);
        signal_last_work_has_started(selected_scheduler);
      });
    } else {
      manager.Predict();
    }
  }
  const bool store_common_info = (facet_indices.front() == 0);
  if (store_common_info) {
    result.unique_id = task.unique_id;
  }
  manager.ProcessResults(result_mutex, result, store_common_info);

  // We are done with this scheduler.
  // Place it in the cache for later re-use.
  if (selected_scheduler) {
    std::lock_guard<std::mutex> lock(scheduler_creation_mutex_);
    scheduler_cache_[n_threads].push_back(std::move(selected_scheduler));
  }
}

}  // namespace wsclean
