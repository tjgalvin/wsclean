#ifndef WSCLEAN_GRIDDING_MS_GRIDDER_MANAGER_H_
#define WSCLEAN_GRIDDING_MS_GRIDDER_MANAGER_H_

#include <mutex>
#include <vector>

#include <aocommon/lane.h>
#include <aocommon/overlappingtaskprocessor.h>
#include <aocommon/taskqueue.h>

#include "h5solutiondata.h"
#include "msgridder.h"
#include "msprovidercollection.h"

#include "../main/settings.h"
#include "../scheduling/griddingresult.h"
#include "../scheduling/griddingtask.h"
#include "../structures/resources.h"
#include "../system/completionsignal.h"

namespace wsclean {

class GriddingTaskManager;

/**
 * Manages lifetime of thread pools and other scheduling related resources that
 * are used by MSGridderManager. This allows for these resources to be reused or
 * shared between MSGridderManager instances when appropriate.
 */
class MSGridderManagerScheduler {
 public:
  MSGridderManagerScheduler(size_t n_workers);
  ~MSGridderManagerScheduler();
  MSGridderManagerScheduler(const MSGridderManagerScheduler&) = delete;
  MSGridderManagerScheduler(MSGridderManagerScheduler&&) = delete;
  MSGridderManagerScheduler& operator=(const MSGridderManagerScheduler&) =
      delete;
  MSGridderManagerScheduler& operator=(MSGridderManagerScheduler&&) = delete;

  aocommon::TaskQueue<std::function<void()>>& GetWorkerTaskQueue() {
    return worker_task_queue_;
  }
  aocommon::OverlappingTaskProcessor& GetOverlappingTaskProcessor() {
    return overlapping_task_processor_;
  }

 private:
  // Create desired number of threads and set them up to read/call tasks from a
  // task queue.
  void InitializeThreadPoolForTaskQueue(
      std::vector<std::thread>& thread_pool_,
      aocommon::TaskQueue<std::function<void()>>& task_queue_,
      size_t n_threads);
  std::vector<std::thread> worker_thread_pool_;
  std::vector<std::thread> scheduler_thread_pool_;
  aocommon::TaskQueue<std::function<void()>> scheduler_task_queue_;
  static constexpr size_t kSchedulerTaskQueueSize_ = 2;
  // All CPU intensive worker tasks should be scheduled through this queue.
  // Non-CPU intensive tasks can be queued through `overlapping_task_processor_`
  // instead so that they don't take resources away from CPU intensive tasks.
  aocommon::TaskQueue<std::function<void()>> worker_task_queue_;
  // Use this pool only to manage scheduling tasks.
  // CPU intensive worker tasks should instead use `worker_task_queue_`.
  // Doing CPU intensive tasks in this queue can lead to oversubscription as
  // `worker_thread_pool_` will already attempt to use all cores that have been
  // made available to us.
  aocommon::OverlappingTaskProcessor overlapping_task_processor_;
};

/**
 * The MSGridderManager is a middle layer between GriddingTaskManager and
 * MsGridder derived classes.
 *
 * GriddingTaskManager is solely responsible for scheduling. MsGridder
 * derived classes are responsible for gridding (inversion/predict)
 *
 * MSGridderManager is responsible for:
 *   1. Managing gridding specific data thats lifecycle/scope is larger than a
 *      single gridder.
 *      e.g. h5parm solution data. Reading a single visibility buffer that
 *      multiple gridders will all operate on.
 *   2. Setting up requirements of the gridder classes.
 *      e.g. initialising MS readers
 *   3. Initialising gridder classes
 *   4. Starting gridder classes
 *   5. Performing final processing of the results before passing them back
 *      through the task manager
 */
class MSGridderManager {
 public:
  MSGridderManager(const Settings& settings,
                   const H5SolutionData& solution_data,
                   MSGridderManagerScheduler* scheduler = nullptr)
      : settings_(settings),
        solution_data_(solution_data),
        scheduler_(scheduler),
        w_limit_(settings.wLimit / 100.0) {}
  ~MSGridderManager() = default;
  MSGridderManager(const MSGridderManager&) = delete;
  MSGridderManager(MSGridderManager&&) = delete;
  MSGridderManager& operator=(const MSGridderManager&) = delete;
  MSGridderManager& operator=(MSGridderManager&&) = delete;

  void InitializeMS(GriddingTask& task);

  void InitializeGridders(GriddingTask& task,
                          const std::vector<size_t>& facet_indices,
                          const Resources& resources,
                          std::vector<GriddingResult::FacetData>& facet_results,
                          GriddingTaskManager* writer_lock_manager);
  void Invert();
  void BatchInvert(std::function<void()> signal_last_gridding_work_has_started);
  void Predict();
  void BatchPredict(std::function<void()> signal_last_predict_work_has_started);
  void ProcessResults(std::mutex& result_mutex, GriddingResult& result,
                      bool store_common_info);

 private:
  /** Call `operation()` once per gridder, for all gridders, running as many as
   * possible in parallel.
   *
   * Make use of all cores/threads available to the manager, each
   * gridder/call consumes 1 thread.
   * This includes threads that would otherwise be assigned to the gridders once
   * they start gridding.
   *
   * NB! Do not use this with an `operation` that will use the gridders internal
   * threads (e.g. predict/invert) as doing so will cause an oversubscription of
   * threads.
   *
   * @param operation Functor taking either one argument (gridder)
   * or two (gridder, task).
   * @param wait_for_idle If true then wait for the task queue to be idle before
   * returning. If false, return immediately. Then it becomes the callers
   * responsibility to wait for the task queue when appropriate.
   */
  template <typename T>
  void ExecuteForAllGridders(T&& operation, bool wait_for_idle = true);

  /** Call `operation()` once per gridder, for all gridders, running as many as
   * possible in parallel while using `n_cores_per_gridder` threads per
   * operation.
   *
   * It makes use of all threads available to the manager, each operation
   * consumes `n_cores_per_gridder` threads, 1 to launch the operation, the
   * remainder are blocked during execution. If the manager manages 32 threads
   * and `n_cores_per_gridder` is 8, then 4 operations will run in parallel.
   *
   * The intention is to use this for operations (predict/invert) where the
   * gridders internal threads will be used. `n_cores_per_gridder` should map to
   * the number of internal gridder threads otherwise there will be an over or
   * under subscription of threads.
   *
   * NB! This call waits/blocks until all the tasks that it queues internally
   * are completed. It does not wait for the task queue itself to be
   * empty.
   * This allows for concurrent `ExecuteForAllGriddersWithNCores()` calls with
   * tasks from the subsequent calls able to start executing prior to all tasks
   * from the first call being complete. This only occurs when remaining tasks
   * from the first core are unable to make use of all threads, this allows
   * better processing throughput in such a scenario.
   *
   * Gridders are protected from concurrent calls:
   * ExecuteForAllGriddersWithNCores2::gridder[0]->operation() will not be
   * called unless ExecuteForAllGriddersWithNCores1::gridder[0]->operation() is
   * already done executing.
   * Ordering is maintained:
   * ExecuteForAllGriddersWithNCores1::gridder2->operation() will always run
   * before ExecuteForAllGriddersWithNCores2::gridder2->operation() but this can
   * only be guaranteed for up to 2 concurrent calls.
   * NB! Never make more than 2 concurrent calls if ordering is important.
   *
   * @param operation Functor taking two arguments (gridder, task).
   * @param signal_first_operation_started
   * ExecuteForAllGriddersWithNCores will call this as a signal to the calling
   * code that the first gridder has acquired a resource lock and will start its
   * processing operation.
   * @param signal_last_operation_started
   * ExecuteForAllGriddersWithNCores will call this as a signal to the calling
   * code that the last gridder has started performing its operation.
   */
  template <typename T>
  void ExecuteForAllGriddersWithNCores(
      size_t n_cores_per_gridder, T&& operation,
      std::function<void()>&& signal_first_operation_started = {},
      std::function<void()>&& signal_last_operation_started = {});

  /**
   * Pointers to data buffers that are required by @ref ReadChunkForInvert to
   * store the chunk of shared data that is computed.
   *
   * Behaviour/size is different depending on whether corrections are to be
   * applied or not. When not applying corrections there is no need to store the
   * antenna pairs or time_offsets_, and as visibilities are collapsed in place
   * less storage is required for visibilities.
   *
   * If applying corrections:
   *   antennas: `max_n_rows`
   *   uvw: `max_n_rows * 3`
   *   visibilities: `max_n_rows * n_channels * n_vis_polarizations`
   *   time_offsets: `max_n_rows`
   *
   * If not applying corrections:
   *   antennas: `0`
   *   uvw: `max_n_rows * 3`
   *   visibilities: `(max_n_rows * n_channels) + (n_channels *
   * n_vis_polarizations)`. The slight overallocation is used to allow
   * collapsing to be done in place without a copy.
   *   time_offsets: `0`
   */
  struct InversionChunkData {
    InversionChunkData(size_t n_rows, size_t n_channels,
                       size_t n_vis_polarisations, bool apply_corrections)
        : antennas(apply_corrections ? n_rows : 0), uvw(n_rows * 3) {
      // If we don't apply corrections then we collapse the visibilities when
      // storing them to save memory. In order to be able to collapse the
      // already copied data in place we have to slightly overallocate the
      // buffer by one uncollapsed row.
      const size_t visibility_size =
          apply_corrections
              ? (n_rows * (n_channels * n_vis_polarisations))
              : ((n_rows * n_channels) + (n_channels * n_vis_polarisations));
      visibilities = aocommon::UVector<std::complex<float>>(visibility_size);
    }
    InversionChunkData() = default;
    aocommon::UVector<std::pair<size_t, size_t>> antennas;
    aocommon::UVector<double> uvw;
    aocommon::UVector<std::complex<float>> visibilities;
    // per row time offset computed during @ref LoadAndApplyCorrections()<kSum>
    // and applied during @ref LoadAndApplyCorrections()<kApply>
    std::vector<size_t> time_offsets;

    size_t gridded_visibility_count;
    double visibility_weight_sum;
    double max_gridded_weight;
    double total_weight;
    size_t n_rows;
  };

  struct PredictionChunkData {
    PredictionChunkData(size_t n_rows)
        : antennas1(n_rows),
          antennas2(n_rows),
          field_ids(n_rows),
          data_desc_ids(n_rows),
          times(n_rows),
          uvws(n_rows * 3) {}
    PredictionChunkData() = default;
    aocommon::UVector<size_t> antennas1;
    aocommon::UVector<size_t> antennas2;
    aocommon::UVector<size_t> field_ids;
    aocommon::UVector<size_t> data_desc_ids;
    aocommon::UVector<double> times;
    aocommon::UVector<double> uvws;
    size_t n_rows = 0;
  };

  /**
   * When populating @ref InversionChunkData it can be optimal to not do so a
   * single row at a time, but instead to gather rows into batches, and then
   * apply necessary processing @ref ApplyWeights() @ref CalculateWeights()
   * to those batches in parallel.
   * @ref BatchRowData facilitates this with each instance representing one
   * batch in such a scenario.
   * The necessarry information is stored such that each @ref BatchRow data
   * can be processed independently and out of order in comparison to other
   * batches.
   * @ref BatchRowData allocates temporary memory for data that will be
   * discarded and not be saved as part of @ref InversionChunkData, data that
   * can't be directly written into @ref InversionChunkData until after
   * processing.
   * @ref BatchRowData uses offsets into the memory of @ref InversionChunkData
   * for data that can be written directly.
   */
  struct BatchRowData {
    BatchRowData(size_t n_rows, size_t n_row_size)
        : metadata_(n_rows),
          visibilities(n_rows,
                       aocommon::UVector<std::complex<float>>(n_row_size)),
          weights(n_rows, aocommon::UVector<float>(n_row_size)),
          model(n_rows, aocommon::UVector<std::complex<float>>(n_row_size)) {}
    BatchRowData() = default;
    std::vector<MSProvider::MetaData> metadata_;
    std::vector<aocommon::UVector<std::complex<float>>> visibilities;
    std::vector<aocommon::UVector<float>> weights;
    std::vector<aocommon::UVector<std::complex<float>>> model;
    size_t n_rows_read = 0;
    size_t time_offsets_offset = 0;
    size_t uvws_offset = 0;
    size_t visibilities_offset = 0;
  };

  /** Read and compute data from an @ref MSReader into a single @ref
   * InversionChunkData chunk which can be passed to @ref BatchInvert() for
   * gridding multiple tasks in parallel.
   * @param [in] ms_reader A @ref MSReader from which the chunk data can be
   * read. Expected to already be set up by the caller.
   * @param [in] selected_buffer Buffer of size `n_channels` containing a
   * boolean determining whether a channel is selected or filtered out.
   * @param [in, out] chunk_data A struct with pointers to buffers of size @ref
   * n_chunk_rows * `data_size` where `data_size` is different for each buffer,
   * see @ref InversionChunkData for more size information.
   * @param [in] shared_data @MsGridderData Initialised by the caller with task
   * and measurement data so that it can be used to call methods a single time
   * for the shared data, instead of these methods having to be called on each
   * individual gridder, methods called are @ref ReadVisibilities,
   * @ref CalculateWeights, @ref StoreImagingWeights and @ref ApplyWeights.
   */
  size_t ReadChunkForInvert(GainMode gain_mode, bool apply_corrections,
                            const std::vector<MsGridder*>& gridders,
                            MsProviderCollection::MsData& ms_data,
                            size_t n_chunk_rows, MSReader& ms_reader,
                            const bool* selected_buffer,
                            InversionChunkData& chunk_data,
                            MsGridderData& shared_data);
  template <GainMode Mode>
  size_t ReadChunkForInvertImplementation(
      size_t n_parms, bool apply_corrections,
      const std::vector<MsGridder*>& gridders,
      MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
      MSReader& ms_reader, const bool* selected_buffer,
      InversionChunkData& chunk_data, MsGridderData& shared_data);
  template <GainMode Mode, size_t NParms>
  size_t ReadChunkForInvertImplementation(
      bool apply_corrections, const std::vector<MsGridder*>& gridders,
      MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
      MSReader& ms_reader, const bool* selected_buffer,
      InversionChunkData& chunk_data, MsGridderData& shared_data);
  template <GainMode Mode, size_t NParms, bool ApplyCorrections>
  size_t ReadChunkForInvertImplementation(
      const std::vector<MsGridder*>& gridders,
      MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
      MSReader& ms_reader, const bool* selected_buffer,
      InversionChunkData& chunk_data, MsGridderData& shared_data);
  template <GainMode Mode, size_t NParms, bool ApplyCorrections, bool ApplyBeam>
  size_t ReadChunkForInvertImplementation(
      const std::vector<MsGridder*>& gridders,
      MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
      MSReader& ms_reader, const bool* selected_buffer,
      InversionChunkData& chunk_data, MsGridderData& shared_data);
  template <GainMode Mode, size_t NParms, bool ApplyCorrections, bool ApplyBeam,
            bool ApplyForward>
  size_t ReadChunkForInvertImplementation(
      const std::vector<MsGridder*>& gridders,
      MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
      MSReader& ms_reader, const bool* selected_buffer,
      InversionChunkData& chunk_data, MsGridderData& shared_data);
  template <GainMode Mode, size_t NParms, bool ApplyCorrections, bool ApplyBeam,
            bool ApplyForward, bool HasH5Parm>
  size_t ReadChunkForInvertImplementation(
      const std::vector<MsGridder*>& gridders,
      MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
      MSReader& ms_reader, const bool* selected_buffer,
      InversionChunkData& chunk_data, MsGridderData& shared_data);

  /**
   * Read and compute data from an @ref MSReader into a single @ref
   * InversionChunkData at a time, using @ref ReadChunkForInvert. Pass the
   * @ref InversionChunkData to the task_lane and then continue reading a new
   * @ref InversionChunkData until all data has been consumed. See @ref
   * ReadChunkForInvert for more information.
   */
  void ReadChunksForInvert(aocommon::Lane<InversionChunkData>& task_lane,
                           size_t n_max_rows_in_memory, bool apply_corrections,
                           MsProviderCollection::MsData& ms_data,
                           MsGridderData& shared_data,
                           const std::vector<MsGridder*>& gridders,
                           size_t n_vis_polarizations,
                           const bool* selected_buffer, size_t& n_total_chunks);

  /**
   * Perform gridding on a single block of data stored in @ref
   * InversionChunkData, for all gridders/facets.
   * @return The number of rows processed.
   */
  size_t GridChunk(bool apply_corrections, size_t n_vis_polarizations,
                   const aocommon::MultiBandData& bands,
                   const InversionChunkData& chunk_data,
                   const aocommon::UVector<double>& frequencies,
                   const MsProviderCollection::MsData& ms_data,
                   size_t chunk_index,
                   std::function<void()> signal_first_gridder_has_started,
                   std::function<void()> signal_last_gridder_has_started);

  /**
   * Perform gridding on a single block of data stored in @ref
   * InversionChunkData, for a single gridder/facet.
   */
  void GridChunkForFacet(MsGridder& gridder, size_t facet_index,
                         bool apply_corrections, size_t n_vis_polarizations,
                         const aocommon::BandData& band, size_t data_desc_id,
                         const InversionChunkData& chunk_data,
                         const aocommon::UVector<double>& frequencies,
                         const MsProviderCollection::MsData& ms_data,
                         size_t chunk_index);

  /**
   * Perform gridding on chunks of @ref InversionChunkData by calling @ref
   * GridChunk() sequentailly on each chunk, as they become available in the
   * task_lane, until all chunks have been processed.
   */
  void GridChunks(
      aocommon::Lane<InversionChunkData>& task_lane, bool apply_corrections,
      const aocommon::UVector<double>& frequencies,
      const aocommon::MultiBandData& bands,
      MsProviderCollection::MsData& ms_data, size_t n_vis_polarizations,
      const size_t& n_chunks,
      std::function<void()> signal_last_gridder_of_last_chunk_has_started);

  /**
   * Read data from an @ref MSReader into a single @ref PredictionChunkData at a
   * time. Pass the filled @ref PredictionChunkData to the task_lane for
   * processing and then continue reading a new @ref PredictionChunkData until
   * all data has been consumed.
   */
  void ReadChunksForPredict(
      aocommon::Lane<PredictionChunkData>& task_lane,
      size_t n_max_rows_in_memory, MsProviderCollection::MsData& ms_data,
      MsGridderData& shared_data, const std::vector<MsGridder*>& gridders,
      const aocommon::BandData band, size_t n_vis_polarizations,
      const bool* selected_buffer, size_t& n_total_chunks);

  /**
   * Perform predict on a single block of data stored in @ref
   * PredictionChunkData, for all gridders/facets
   */
  size_t PredictChunk(const PredictionChunkData& chunk_data,
                      size_t n_vis_polarizations, size_t n_antennas,
                      std::vector<std::complex<float>>& combined_visibilities,
                      const aocommon::UVector<double>& frequencies,
                      MsProviderCollection::MsData& ms_data,
                      std::function<void()> signal_first_predict_has_started,
                      std::function<void()> signal_last_predict_has_started);

  /**
   * Perform predict on a single block of data stored in @ref
   * PredictionChunkData, for a single gridder/facet.
   */
  void PredictChunkForFacet(
      MsGridder& gridder, size_t facet_index,
      const PredictionChunkData& chunk_data, size_t n_vis_polarizations,
      size_t n_antennas,
      std::vector<std::complex<float>>& combined_visibilities,
      const aocommon::UVector<double>& frequencies,
      MsProviderCollection::MsData& ms_data,
      std::vector<std::mutex>& sum_visibilities_mutexes,
      size_t mutex_chunk_size);

  /**
   * Perform predict on chunks of @ref PredictionChunkData by calling @ref
   * PredictChunk() sequentially on each chunk, as they become available in the
   * task_lane, until all chunks have been processed.
   */
  void PredictChunks(
      aocommon::Lane<PredictionChunkData>& task_lane,
      const aocommon::UVector<double>& frequencies,
      MsProviderCollection::MsData& ms_data, size_t n_vis_polarizations,
      bool add_assign_model, const size_t& n_total_chunks,
      std::function<void()> signal_last_predict_of_last_chunk_has_started);

  std::unique_ptr<MsGridder> ConstructGridder(const Resources& resources);
  struct GriddingFacetTask {
    std::unique_ptr<MsGridder> facet_gridder;
    GriddingTask::FacetData* facet_task;
    GriddingResult::FacetData* facet_result;
  };
  std::vector<GriddingFacetTask> facet_tasks_;

  inline void InitializeMSDataVectors() {
    std::vector<MsGridder*> gridders;
    gridders.reserve(facet_tasks_.size());
    for (auto& [gridder, facet_task, facet_result] : facet_tasks_) {
      gridders.push_back(gridder.get());
    }
    ms_provider_collection_.InitializeMSDataVector(gridders, w_limit_,
                                                   solution_data_.HasData());
  }

  /** Initializes 'gridder' with values that are equal for all facets. */
  void InitializeGridderForTask(MsGridder& gridder, const GriddingTask& task,
                                GriddingTaskManager* writer_lock_manager);

  /** Initializes 'gridder' with facet-specific values. */
  void InitializeGridderForFacet(bool is_first_polarization, MsGridder& gridder,
                                 GriddingTask::FacetData& facet_task);

  const Settings& settings_;
  const H5SolutionData& solution_data_;
  MSGridderManagerScheduler* scheduler_;
  MsProviderCollection ms_provider_collection_;
  int64_t available_memory_ = 0;
  size_t available_cores_ = 0;
  size_t available_cores_per_gridder_ = 0;
  /// A fractional value that, when non-zero, places a limit on the w-value of
  /// gridded visibilities. Visibilities outside the limit are skipped.
  double w_limit_ = 0.0;
};

template <typename T>
void MSGridderManager::ExecuteForAllGridders(T&& operation,
                                             bool wait_for_idle) {
  for (GriddingFacetTask& task : facet_tasks_) {
    MsGridder* gridder = task.facet_gridder.get();
    if constexpr (std::is_invocable<T, MsGridder*, GriddingFacetTask&>::value) {
      scheduler_->GetWorkerTaskQueue().Emplace(
          [=, &task]() { operation(gridder, task); });
    } else {
      scheduler_->GetWorkerTaskQueue().Emplace([=]() { operation(gridder); });
    }
  }
  if (wait_for_idle) {
    scheduler_->GetWorkerTaskQueue().WaitForIdle(available_cores_);
  }
}

template <typename T>
void MSGridderManager::ExecuteForAllGriddersWithNCores(
    size_t n_cores_per_gridder, T&& operation,
    std::function<void()>&& signal_first_operation_started,
    std::function<void()>&& signal_last_operation_started) {
  std::vector<std::shared_ptr<CompletionSignal>> signals;
  signals.reserve(facet_tasks_.size() + 1);
  // Queue tasks in a way that each task will consume `n_cores_per_gridder` task
  // slots.
  for (const GriddingFacetTask& task : facet_tasks_) {
    MsGridder* gridder = task.facet_gridder.get();

    // Avoid gridder processing concurrently.
    gridder->processing_semaphore_.acquire();

    // Signal to @ref OverlappingTaskProcessor that its okay for the next
    // overlapping task to start without any risk of a race, as we have
    // acquired the gridder semaphore. Ordering is maintained because as this
    // task holds the semaphore the other task cannot overtake it.
    if (signal_first_operation_started) {
      signal_first_operation_started();
      signal_first_operation_started = {};
    }
    // Task must consume N threads.
    // Pause N-1 threads so that they are unavailable to the pool.
    // Only once we have done so are we allowed to run the task.
    std::shared_ptr<CompletionSignal> signal =
        std::make_shared<CompletionSignal>();
    signals.push_back(signal);

    // Queue "blocker" tasks to consume extra task slots.
    for (size_t i = 0; i < n_cores_per_gridder - 1; ++i) {
      // Consume 1 thread from the pool until completion is signalled.
      scheduler_->GetWorkerTaskQueue().Emplace(
          [=]() { signal->WaitForCompletion(); });
    }
    // Queue task to perform the operation.
    const size_t index = task.facet_task->index;
    scheduler_->GetWorkerTaskQueue().Emplace([=]() {
      // We are consuming N threads from the pool while this occurs.
      operation(gridder, index);
      gridder->processing_semaphore_.release();
      // Return the additional threads to the pool.
      signal->SignalCompletion();
    });
  }
  // Signal to higher level schedulers (e.g. @ref ThreadedScheduler) that we
  // have reached a point where this ExecuteForAllGriddersWithNCores call will
  // no longer be able to fully occupy all of our threading resources. This
  // gives the scheduler an opportunity to schedule additional work on our queue
  // to make use of the otherwise idle threading resources.
  if (signal_last_operation_started) {
    std::shared_ptr<CompletionSignal> signal =
        std::make_shared<CompletionSignal>();
    signals.push_back(signal);
    scheduler_->GetWorkerTaskQueue().Emplace([=]() {
      signal_last_operation_started();
      signal->SignalCompletion();
    });
  }
  // Wait for all tasks launched by this call to complete.
  for (const auto& signal : signals) {
    signal->WaitForCompletion();
  }
}

}  // namespace wsclean

#endif  // WSCLEAN_GRIDDING_MS_GRIDDER_MANAGER_H_
