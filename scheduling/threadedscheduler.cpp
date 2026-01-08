#include "threadedscheduler.h"
#include "../gridding/msgridder.h"
#include "../gridding/msgriddermanager.h"

#include "../main/settings.h"

#include <aocommon/logger.h>

#include <string>

using aocommon::Logger;

namespace wsclean {

ThreadedScheduler::ThreadedScheduler(const Settings& settings)
    : GriddingTaskManager{settings},
      // When using the ThreadedScheduler as the main scheduler, limit the
      // number of tasks in the queue to one per thread. When stacking too many
      // tasks, memory usage could become an issue.
      // When using the ThreadedScheduler as a local scheduler with the
      // MPIScheduler, the MPIScheduler manages the task distribution.
      // The ThreadedScheduler should always queue new tasks in that case.
      task_queue_(settings.UseMpi() ? TaskQueueType()
                                    : TaskQueueType(settings.parallelGridding)),
      resources_per_task_(GetResources().GetPart(settings.parallelGridding)) {
  Logger::Debug << "[ThreadedScheduler] Using " +
                       std::to_string(settings.parallelGridding) +
                       " workers with " +
                       std::to_string(resources_per_task_.NCpus()) +
                       " threads per worker.\n";
  for (size_t i = 0; i < settings.parallelGridding; ++i) {
    thread_list_.emplace_back(&ThreadedScheduler::ProcessQueue, this);
  }
}

ThreadedScheduler::~ThreadedScheduler() {
  try {
    Finish();
  } catch (std::exception& e) {
    // Normally, the user of the ThreadedScheduler calls Finish(), which
    // rethrows any exception caught in a thread.
    // We are in a destructor, so all that can be done is report the error.
    using namespace std::string_literals;
    Logger::Error
        << "Exception caught during destruction of ThreadedScheduler:\n"s +
               e.what() + '\n';
  }
  task_queue_.Finish();  // Make all threads exit.
  for (std::thread& thread : thread_list_) thread.join();
}

void ThreadedScheduler::Run(
    GriddingTask&& task, std::function<void(GriddingResult&)> finish_callback) {
  const std::size_t task_id = task.unique_id;
  const std::size_t facet_count = task.facets.size();
  TaskData* task_data;

  // Add an entry in task_data_map_ for the task.
  {
    std::lock_guard<std::mutex> lock(mutex_);
    assert(task_data_map_.count(task_id) == 0);
    task_data = &task_data_map_[task_id];
  }
  task_data->task = std::move(task);
  task_data->result.facets.resize(facet_count);
  task_data->callback = std::move(finish_callback);

  const bool batch_invert = GetSettings().shared_facet_reads &&
                            task.operation == GriddingTask::Invert;
  const bool batch_predict =
      (GetSettings().shared_facet_reads || GetSettings().shared_facet_writes) &&
      task.operation == GriddingTask::Predict;
  if (!batch_invert && !batch_predict) {
    // Add sub-tasks for each facet to the task queue.
    for (std::size_t facet_index = 0; facet_index < facet_count;
         ++facet_index) {
      task_queue_.Emplace(task_id, std::vector<size_t>{facet_index});
    }
  } else {
    // Set the amount of parallel gridders for this task
    task_data->task.num_parallel_gridders_ = std::min(
        GetSettings().parallelGridding, std::max(facet_count, size_t{1}));

    // Each facet group is only one task, the task will handle splitting
    // resources between facets internally If we are not using compound tasks
    // then we end up here for each individual facet and thereby create one task
    // per facet
    if (facet_count > 1) {
      // Wait tasks call WaitForCompletion() on
      // lock_excess_scheduler_tasks_ which will pause/freeze/block them until
      // SignalCompletion() is called
      task_data->task.lock_excess_scheduler_tasks_ =
          std::make_shared<CompletionSignal>();

      // There must be N-1 blocker tasks before we can run, giving us N threads
      // in total as requested Set up locks so that we only start operating once
      // we have the resources The scheduler will allocate us resources when the
      // dummy tasks are placed in the queue
      for (size_t wait_task_index = 0;
           wait_task_index < task_data->task.num_parallel_gridders_ - 1;
           ++wait_task_index) {
        // Exact id doesn't really matter we just need to ensure that id is
        // unique
        const size_t wait_task_id = std::numeric_limits<size_t>::max() -
                                    (task_data->task.unique_id * 1000) -
                                    wait_task_index;
        TaskData* wait_task_data;
        {
          std::lock_guard<std::mutex> lock(mutex_);
          assert(task_data_map_.count(wait_task_id) == 0);
          wait_task_data = &task_data_map_[wait_task_id];
        }

        GriddingTask wait_task;
        wait_task.operation = GriddingTask::Wait;
        wait_task.lock_excess_scheduler_tasks_ =
            task_data->task.lock_excess_scheduler_tasks_;
        wait_task_data->task = std::move(wait_task);
        task_queue_.Emplace(wait_task_id, std::vector<size_t>{0});
      }

      std::vector<size_t> facet_indexes(facet_count);
      std::iota(facet_indexes.begin(), facet_indexes.end(), 0);
      assert(!facet_indexes.empty());
      task_queue_.Emplace(task_id, facet_indexes);
    } else {
      assert(0);
      task_queue_.Emplace(task_id, std::vector<size_t>{0});
    }
  }

  ProcessReadyList();
}

bool ThreadedScheduler::TryStealTask(
    std::shared_ptr<CompletionSignal>& parent_task_completion_signal,
    std::shared_ptr<CompletionSignal>& signal_stolen_task_completion,
    std::unique_ptr<MSGridderManagerScheduler>& scheduler, size_t n_work_size) {
  assert(n_work_size > 1);

  Logger::Debug << "[TaskSteal] Attempting to steal.\n";

  std::lock_guard<std::mutex> lock(mutex_);
  std::vector<std::pair<size_t, std::vector<size_t>>> facet_tasks;
  // We aren't just taking the worker task but also all the "Wait" tasks
  // associated with it as well. Grab them all in one go.
  if (task_queue_.TryPopN(facet_tasks, n_work_size)) {
    for (size_t i = 0; i < n_work_size - 1; ++i) {
      const size_t task_id = facet_tasks[i].first;
      TaskData& task_data = task_data_map_[task_id];
      if (task_data.task.operation == GriddingTask::Wait) {
        // We take over ownership of the existing wait task and re-use that.
        // So immediately remove the new wait task as if it is completed.
        task_data_map_.erase(task_id);
      } else {
        Logger::Debug << "[TaskSteal] Fatal error, wait tasks not in expected "
                         "order in queue.\n";
        throw std::runtime_error(
            "Invalid task order encountered in "
            "ThreadedScheduler::TryStealTask");
      }
    }
    const size_t& task_id = facet_tasks[n_work_size - 1].first;
    const std::vector<size_t>& facet_indexes =
        facet_tasks[n_work_size - 1].second;
    TaskData& task_data = task_data_map_[task_id];
    if (task_data.task.operation != GriddingTask::Wait) {
      Logger::Debug << "[TaskSteal] Success.\n";
      // Take over ownership of wait tasks from parent task.
      task_data.task.lock_excess_scheduler_tasks_ =
          parent_task_completion_signal;
      parent_task_completion_signal.reset();
      signal_stolen_task_completion =
          task_data.task.lock_excess_scheduler_tasks_;
      std::thread([&, facet_indexes, task_id]() {
        ProcessTask(task_id, facet_indexes, task_data,
                    task_data.task.lock_excess_scheduler_tasks_,
                    std::move(scheduler));
      }).detach();
    } else {
      Logger::Debug << "[TaskSteal] Fatal error, worker tasks not in expected "
                       "order in queue.\n";
      throw std::runtime_error(
          "Invalid task order encountered in ThreadedScheduler::TryStealTask");
    }
    return true;
  }
  Logger::Debug << "[TaskSteal] Failed to steal, no tasks available.\n";
  return false;
}

void ThreadedScheduler::ProcessTask(
    size_t task_id, const std::vector<size_t>& facet_indexes,
    TaskData& task_data,
    std::shared_ptr<CompletionSignal>& signal_task_completion,
    std::unique_ptr<MSGridderManagerScheduler> scheduler) {
  const size_t task_size = task_data.task.num_parallel_gridders_;

  // If we steal a task the stolen task will set this signal and notify with it
  // when it reaches completion.
  std::shared_ptr<CompletionSignal> signal_stolen_task_completion = nullptr;
  try {
    assert(!facet_indexes.empty());
    // As the gridder manager will potentially be running/managing N parallel
    // gridding tasks internally instead of just a single one it is necessary
    // to allocate it the appropriate resources for all N tasks.
    Resources task_resources = resources_per_task_.GetCombined(task_size);
    if (task_data.task.operation == GriddingTask::Wait) {
      // Wait tasks occupy a thread from the pool by waiting on
      // a completion signal.
      assert(signal_task_completion);
      signal_task_completion->WaitForCompletion();
    } else {
      RunDirect(
          task_data.task, facet_indexes, task_resources, task_data.result,
          task_data.result_mutex,
          [&](std::unique_ptr<MSGridderManagerScheduler>& parent_scheduler) {
            TryStealTask(signal_task_completion, signal_stolen_task_completion,
                         parent_scheduler, task_size);
          },
          std::move(scheduler));
      if (signal_task_completion) {
        signal_task_completion->SignalCompletion();
      }
    }
  } catch (std::exception&) {
    std::lock_guard<std::mutex> lock(mutex_);
    latest_exception_ = std::current_exception();
  }

  if (task_data.task.operation == GriddingTask::Wait) {
    std::lock_guard<std::mutex> lock(mutex_);
    task_data_map_.erase(task_id);
    return;
  }

  // Extract the new value from task_data.finished_facet_count directly.
  // When extracting the new value later, another thread may also have
  // incremented the atomic value in the mean time, and multiple threads
  // will think they have processed the last facet.
  bool process = false;
  {
    std::lock_guard<std::mutex> result_lock(task_data.result_mutex);
    task_data.finished_facet_count += facet_indexes.size();
    process = (task_data.finished_facet_count == task_data.task.facets.size());
  }
  if (process) {
    if (GetSettings().UseMpi()) {
      // Execute callback immediately, from the processing thread.
      // The MPIScheduler stores the result at the main node.
      // The MPIWorkerScheduler sends the result to the main node.
      task_data.callback(task_data.result);
      std::lock_guard<std::mutex> lock(mutex_);
      task_data_map_.erase(task_id);
    } else {
      // Store the task id and execute the callback on the main thread.
      std::lock_guard<std::mutex> lock(mutex_);
      ready_list_.emplace_back(task_id);
    }
  }

  // Wait until stolen task is done, otherwise we end up with an extra
  // available thread in the pool.
  if (signal_stolen_task_completion) {
    Logger::Debug << "[ProcessTask] Wait for stolen task.\n";
    signal_stolen_task_completion->WaitForCompletion();
  }
}

void ThreadedScheduler::ProcessQueue() {
  std::pair<size_t, std::vector<size_t>> facet_task_pair;
  while (task_queue_.Pop(facet_task_pair)) {
    const size_t task_id = facet_task_pair.first;
    const std::vector<size_t>& facet_indexes = facet_task_pair.second;
    TaskData& task_data = task_data_map_[task_id];
    std::shared_ptr<CompletionSignal>& signal_task_completion =
        task_data.task.lock_excess_scheduler_tasks_;
    ProcessTask(task_id, facet_indexes, task_data, signal_task_completion);
  }
}

void ThreadedScheduler::Start(size_t nWriterGroups) {
  assert(ready_list_.empty());
  GriddingTaskManager::Start(nWriterGroups);
  if (writer_group_locks_.size() < nWriterGroups)
    writer_group_locks_ = std::vector<std::mutex>(nWriterGroups);
}

std::unique_ptr<GriddingTaskManager::WriterLock> ThreadedScheduler::GetLock(
    size_t writer_group_index) {
  assert(writer_group_index < writer_group_locks_.size());
  return std::make_unique<ThreadedWriterLock>(*this, writer_group_index);
}

void ThreadedScheduler::Finish() {
  task_queue_.WaitForIdle(GetSettings().parallelGridding);

  ProcessReadyList();
}

void ThreadedScheduler::ProcessReadyList() {
  std::vector<std::size_t> local_ready_list;
  std::exception_ptr local_exception;

  // Move the ready_list_ and latest_exception_ to local variables so we can
  // release the lock while the callbacks run / while throwing the exception.
  {
    std::lock_guard<std::mutex> lock{mutex_};

    local_ready_list = std::move(ready_list_);
    ready_list_.clear();

    local_exception = std::move(latest_exception_);
    latest_exception_ = std::exception_ptr();
  }

  // Check for exceptions before calling callbacks, since results
  // are typically invalid when an exception occurred.
  if (local_exception) std::rethrow_exception(local_exception);

  // Call callbacks for any finished tasks
  for (std::size_t task_id : local_ready_list) {
    TaskData& task_data = task_data_map_[task_id];
    task_data.callback(task_data.result);
  }

  // Remove the finished tasks from task_data_map_.
  if (!local_ready_list.empty()) {
    std::lock_guard<std::mutex> lock{mutex_};
    for (std::size_t task_id : local_ready_list) {
      task_data_map_.erase(task_id);
    }
  }
}

}  // namespace wsclean
