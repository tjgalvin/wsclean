#ifndef WSCLEAN_GRIDDING_MS_GRIDDER_MANAGER_H_
#define WSCLEAN_GRIDDING_MS_GRIDDER_MANAGER_H_

#include <mutex>
#include <vector>

#include <aocommon/lane.h>
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
 * The MSGridderManager is a middle layer between GriddingTaskManager and
 * MsGridder derived classes.
 *
 * GriddingTaskManager is solely responsible for scheduling MsGridder
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
                   const H5SolutionData& solution_data)
      : settings_(settings),
        solution_data_(solution_data),
        w_limit_(settings.wLimit / 100.0) {}
  ~MSGridderManager() {}

  MSGridderManager(const MSGridderManager&) = delete;
  MSGridderManager& operator=(const MSGridderManager&) = delete;

  void InitializeMS(GriddingTask& task);

  void InitializeGridders(GriddingTask& task,
                          const std::vector<size_t>& facet_indices,
                          const Resources& resources,
                          std::vector<GriddingResult::FacetData>& facet_results,
                          GriddingTaskManager* writer_lock_manager);
  void Invert();
  void BatchInvert(size_t num_parallel_gridders);
  void Predict();
  void ProcessResults(std::mutex& result_mutex, GriddingResult& result,
                      bool store_common_info);

  /**
   * Sort facet tasks by expected gridding time, longest first.
   * This enables a mninor optimization for shared reads by slightly reducing
   * the wait time of idle cores when only a few gridders are left running at
   * the end of a batch.
   */
  void SortFacetTasks();

 private:
  /** Execute `operation` for all gridders in parallel, using all cores/threads
   * available to the manager. This includes threads that would otherwise be
   * assigned to the gridders once they start gridding.
   * @param wait_for_idle If true then wait for all gridders to finish executing
   * the operation before returning, if false then return immediately and it
   * becomes the callers responsibility to wait when appropriate.
   */
  template <typename T>
  void ExecuteForAllGridders(
      aocommon::TaskQueue<std::function<void()>>& task_queue, T&& operation,
      bool wait_for_idle = true);

  /** Execute `operation` for all gridders in parallel, using fewer cores than
   * the total available to the manager. This is used for gridding, because once
   * the gridders start we only want to use `num_parallel_gridders` worth of
   * threads/cores. The remainder are used internally inside the gridders.
   */
  template <typename T>
  void ExecuteForAllGriddersWithNCores(
      aocommon::TaskQueue<std::function<void()>>& task_queue, size_t n_cores,
      T&& operation);

  /** Pointers to data buffers that are required by @ref ReadChunkForInvert for
   * computation. These should be initialised by the called to a size of
   * `n_channels * n_vis_polarizations`
   */
  struct RowData {
    std::complex<float>* model;
    float* weights;
  };
  /** Pointers to data buffers that are required by @ref ReadChunkForInvert to
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
  struct ChunkData {
    ChunkData(size_t n_rows, size_t n_channels, size_t n_vis_polarisations,
              bool apply_corrections)
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
    ChunkData() = default;
    aocommon::UVector<std::pair<size_t, size_t>> antennas;
    aocommon::UVector<double> uvw;
    aocommon::UVector<std::complex<float>> visibilities;
    // per row time offset computed during @ref ApplyCorrections()<kSum>
    // and applied during @ref ApplyCorrections()<kApply>
    std::vector<size_t> time_offsets;

    size_t gridded_visibility_count;
    size_t visibility_weight_sum;
    size_t max_gridded_weight;
    size_t total_weight;
    size_t n_rows;
  };

  /** Read and compute data from an @ref MSReader into a single @ref ChunkData
   * chunk which can be passed to @ref BatchInvert for gridding multiple
   * gridders in parallel.
   * @param [in] task_queue A task queue that is used to call @ref
   * ApplyCorrections in parallel across multiple gridders.
   * @param [in] ms_reader A @ref MSReader from which the chunk data can be
   * read. Expected to already be set up by the caller.
   * @param [in] selected_buffer Buffer of size `n_channels` containing a
   * boolean determining whether a channel is selected or filtered out.
   * @param [in] row_data Struct containing per row buffers of size `n_channels
   * n_vis_polarizations` that are used to hold/calculate temporary data while
   * populating the chunk.
   * @param [in, out] chunk_data A struct with pointers to buffers of size @ref
   * n_chunk_rows * `data_size` where `data_size` is different for each buffer,
   * see @ref ChunkData for more size information.
   * @param [in] shared_data @MsGridderData Initialised by the caller with task
   * and measurement data so that it can be used to call methods a single time
   * for the shared data, instead of these methods having to be called on each
   * individual gridder, methods called are @ref ReadVisibilities,
   * @ref CalculateWeights, @ref StoreImagingWeights and @ref ApplyWeights.
   */
  size_t ReadChunkForInvert(
      GainMode gain_mode, bool apply_corrections,
      aocommon::TaskQueue<std::function<void()>>& task_queue,
      const std::vector<MsGridder*>& gridders,
      MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
      MSReader& ms_reader, const aocommon::BandData band,
      const bool* selected_buffer, RowData& row_data, ChunkData& chunk_data,
      MsGridderData& shared_data);
  template <GainMode Mode>
  size_t ReadChunkForInvertImplementation(
      bool apply_corrections,
      aocommon::TaskQueue<std::function<void()>>& task_queue,
      const std::vector<MsGridder*>& gridders,
      MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
      MSReader& ms_reader, const aocommon::BandData band,
      const bool* selected_buffer, RowData& row_data, ChunkData& chunk_data,
      MsGridderData& shared_data);

  /**
   * Read and compute data from an @ref MSReader into a single @ref ChunkData
   * at a time, using @ref ReadChunkForInvert. Pass the ChunkData to the
   * task_lane and then continue reading a new ChunkData until all data has been
   * consumed. See @ref ReadChunkForInvert for more information.
   */
  void ReadChunksForInvert(
      aocommon::Lane<ChunkData>& task_lane,
      aocommon::TaskQueue<std::function<void()>>& task_queue,
      size_t n_max_rows_in_memory, bool apply_corrections,
      MsProviderCollection::MsData& ms_data, MsGridderData& shared_data,
      const std::vector<MsGridder*>& gridders, const aocommon::BandData band,
      size_t n_vis_polarizations, const bool* selected_buffer);

  /**
   * Perform gridding on a single block of data stored in @ref ChunkData/
   */
  void GridChunk(size_t n_rows, size_t num_parallel_gridders,
                 bool apply_corrections, ChunkData& chunk_data,
                 std::vector<MsGridder*>& gridders,
                 size_t gridded_visibility_count, size_t visibility_weight_sum,
                 size_t max_gridded_weight, size_t total_weight,
                 size_t n_vis_polarizations,
                 aocommon::TaskQueue<std::function<void()>>& task_queue,
                 const aocommon::UVector<double>& frequencies,
                 const aocommon::BandData& band,
                 MsProviderCollection::MsData& ms_data);

  /**
   * Perform gridding on chunks of @ref ChunkData by calling @ref GridChunk
   * sequentailly on each chunk, as they become available in the task_lane,
   * until all chunks have been processed.
   */
  void GridChunks(aocommon::Lane<ChunkData>& task_lane,
                  const size_t num_parallel_gridders,
                  const bool apply_corrections,
                  std::vector<MsGridder*>& gridders,
                  aocommon::TaskQueue<std::function<void()>>& task_queue,
                  const aocommon::UVector<double>& frequencies,
                  const aocommon::BandData& band,
                  MsProviderCollection::MsData& ms_data,
                  size_t n_vis_polarizations);

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
  void InitializeGridderForFacet(MsGridder& gridder,
                                 GriddingTask::FacetData& facet_task);

  const Settings& settings_;
  const H5SolutionData& solution_data_;
  MsProviderCollection ms_provider_collection_;
  int64_t available_memory_;
  size_t available_cores_;
  /// A fractional value that, when non-zero, places a limit on the w-value of
  /// gridded visibilities. Visibilities outside the limit are skipped.
  double w_limit_ = 0.0;
};

template <typename T>
void MSGridderManager::ExecuteForAllGridders(
    aocommon::TaskQueue<std::function<void()>>& task_queue, T&& operation,
    bool wait_for_idle) {
  for (const GriddingFacetTask& task : facet_tasks_) {
    MsGridder* gridder = task.facet_gridder.get();
    task_queue.Emplace([=]() { operation(gridder); });
  }
  if (wait_for_idle) {
    task_queue.WaitForIdle(available_cores_);
  }
}

template <typename T>
void MSGridderManager::ExecuteForAllGriddersWithNCores(
    aocommon::TaskQueue<std::function<void()>>& task_queue, size_t n_cores,
    T&& operation) {
  // Pause the excess threads in the pool that we don't want to use for this
  // operation.
  const size_t excess_cores = available_cores_ - n_cores;
  CompletionSignal signal;
  for (size_t i = 0; i < excess_cores; ++i) {
    task_queue.Emplace([&]() { signal.WaitForCompletion(); });
  }

  // Run the operation with the reduced quantity of available threads.
  for (const GriddingFacetTask& task : facet_tasks_) {
    MsGridder* gridder = task.facet_gridder.get();
    const size_t index = task.facet_task->index;
    task_queue.Emplace([=]() { operation(gridder, index); });
  }
  task_queue.WaitForIdle(n_cores);

  // Restore the quantity of available threads to as it was before entering this
  // method.
  signal.SignalCompletion();

  // We must let all threads process the signal before its destructor is called.
  task_queue.WaitForIdle(available_cores_);
}

}  // namespace wsclean

#endif  // WSCLEAN_GRIDDING_MS_GRIDDER_MANAGER_H_
