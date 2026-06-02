#include "msgriddermanager.h"

#include <functional>
#include <mutex>
#include <new>
#include <vector>

#include <aocommon/logger.h>
#include <aocommon/taskqueue.h>
#include <aocommon/uvector.h>

#include <schaapcommon/facets/facet.h>
#include <schaapcommon/h5parm/h5parm.h>
#include <schaapcommon/h5parm/jonesparameters.h>
#include <schaapcommon/h5parm/soltab.h>
#include <schaapcommon/threading/threadpool.h>

#include "directmsgridder.h"
#include "griddingutils.h"
#include "h5solutiondata.h"
#include "msgridder.h"
#include "msprovidercollection.h"
#include "wsmsgridder.h"

#include "../idg/averagebeam.h"
#include "../idg/idgmsgridder.h"
#include "../idg/facetidgmsgridder.h"
#include "../main/settings.h"
#include "../structures/numberlist.h"
#include "../structures/resources.h"
#include "../wgridder/wgriddingmsgridder.h"
#include "../wtowers/wtowersmsgridder.h"

using aocommon::Logger;

namespace wsclean {

namespace {

// Prevent false sharing between threads by isolating mutexes and atomics on their own cache lines.
//
// When multiple threads frequently update adjacent variables in an array, they compete for
// ownership of the same cache line (false sharing), causing severe performance degradation.
// Using `std::hardware_destructive_interference_size` ensures each variable is
// spaced according to the CPU's cache line size (typically 64 bytes).
//
// We use wrapper structs because `std::mutex` and `std::atomic` have
// internal alignment requirements that can conflict with direct `alignas` application.
// The wrapper allows the compiler to satisfy both the cache line alignment and the type's own
// alignment requirements through padding.
//
// Important: These types must be allocated using std::make_unique_for_overwrite<T[]>(size)
// rather than std::vector<T>. std::vector does not guarantee that over-aligned types
// (alignas > alignof(std::max_align_t)) will respect the requested alignment. The C++
// standard only requires vector to honor alignments up to the platform's maximum scalar
// alignment, which is typically 16 bytes. Using make_unique_for_overwrite ensures the
// allocator respects our cache line alignment (64+ bytes) and avoids default-initialization
// overhead for large arrays.
//
// These wrappers were introduced to resolve bottlenecks in Predict/Invert operations,
// resulting in a measured 2.46x speedup for Predict on four nodes running with a 189 GB
// simulation.
#if defined(__cpp_lib_hardware_interference_size) && __cpp_lib_hardware_interference_size >= 201603
inline constexpr size_t kInterferenceSize = std::hardware_destructive_interference_size;
#else
#warning std::hardware_destructive_interference_size not supported; guessing cache size.
inline constexpr size_t kInterferenceSize = 64;
#endif
struct alignas(kInterferenceSize) CacheAlignedAtomic {
  std::atomic<size_t> value{0};
};
struct alignas(kInterferenceSize) CacheAlignedMutex {
  std::mutex value;
};

/**
 * When populating @ref InversionChunkData it can be optimal to not do so a
 * single row at a time, but instead to gather rows into batches, and then
 * apply necessary processing @ref ApplyWeights() @ref CalculateWeights()
 * to those batches in parallel.
 * @ref BatchRowData facilitates this with each instance representing one
 * batch in such a scenario.
 * The necessary information is stored such that each @ref BatchRow data
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
      : metadata(n_rows),
        visibilities(n_rows * n_row_size),
        weights(n_rows * n_row_size),
        model(n_rows * n_row_size) {}

  BatchRowData() = default;
  BatchRowData(const BatchRowData&) = delete;
  BatchRowData& operator=(const BatchRowData&) = delete;
  BatchRowData(BatchRowData&& other) noexcept = default;
  BatchRowData& operator=(BatchRowData&&) = default;

  aocommon::UVector<MSProvider::MetaData> metadata;
  uint32_t data_desc_id;
  aocommon::UVector<std::complex<float>> visibilities;
  aocommon::UVector<float> weights;
  aocommon::UVector<std::complex<float>> model;

  size_t n_rows = 0;
  size_t time_offsets_out_offset = 0;
  size_t uvws_out_offset = 0;
  size_t visibilities_out_offset = 0;
  size_t index = 0;
  std::vector<std::tuple<double, size_t, size_t>> beam_cache_offsets;
};

/**
 * Simple view class to assist with iterating over the data in a @ref BatchRowData
 */
struct BatchRowView {
 public:
  BatchRowView(BatchRowData* data, size_t model_stride, size_t weights_stride,
               size_t visibilies_stride, size_t uvws_stride, size_t n_rows_to_read)
      : data_(data),
        metadata_(data->metadata.data()),
        model_(data->model.data()),
        visibilities_(data->visibilities.data()),
        weights_(data->weights.data()),
        model_stride_(model_stride),
        visibilities_stride_(visibilies_stride),
        weights_stride_(weights_stride),
        n_rows_to_read_(n_rows_to_read) {
    data_->n_rows = 0;
    current_chunk_row_ = data_->uvws_out_offset / uvws_stride;
  }
  bool End() {
    if (data_->n_rows >= n_rows_to_read_) return true;
    return false;
  }
  void NextRow() {
    ++data_->n_rows;
    ++current_chunk_row_;
    ++metadata_;
    model_ += model_stride_;
    weights_ += weights_stride_;
    visibilities_ += visibilities_stride_;
  }
  bool ReadRow(MSReader& ms_reader, MsGridderData& gridder_data) {
    if (!ms_reader.CurrentRowAvailable() || End()) {
      return false;
    }

    ms_reader.ReadMeta(*metadata_);
    gridder_data.ReadVisibilities(ms_reader, visibilities_, weights_, model_);
    // TODO: Can't use bulk inversion with BDA yet
    data_->data_desc_id = metadata_->data_desc_id;
    return true;
  }

  BatchRowData* data_;
  MSProvider::MetaData* metadata_;
  std::complex<float>* model_;
  std::complex<float>* visibilities_;
  float* weights_;

  size_t current_chunk_row_;

 private:
  const size_t model_stride_;
  const size_t visibilities_stride_;
  const size_t weights_stride_;
  const size_t n_rows_to_read_;
};

/**
 * Simple view class to assist with iterating over the data in a @ref InversionChunkData
 */
template <bool ApplyCorrections>
struct InversionChunkView {
 public:
  InversionChunkView(InversionChunkData* data, const BatchRowData& rows, size_t uvws_stride,
                     size_t visibilies_stride)
      : data_(data),
        time_offsets_(data->time_offsets.data() + rows.time_offsets_out_offset),
        uvws_(data->uvw.data() + rows.uvws_out_offset),
        visibilities_(data->visibilities.data() + rows.visibilities_out_offset),
        uvws_stride_(uvws_stride),
        visibilities_stride_(visibilies_stride) {}
  void NextRow() {
    if constexpr (ApplyCorrections) {
      ++time_offsets_;
    }
    uvws_ += uvws_stride_;
    visibilities_ += visibilities_stride_;
  }

  InversionChunkData* data_;
  size_t* time_offsets_;
  double* uvws_;
  std::complex<float>* visibilities_;

 private:
  const size_t uvws_stride_;
  const size_t visibilities_stride_;
};

/**
 * Perform sequential operations required on each row:
 * 1. Prepare metadata for parallel correction computations.
 * 2. Pre-allocate space for beam response updates.
 * 3. Cache data for @ref ComputeBeamResponseForRowInOrder to compute response updates.
 * 4. Cache h5parm and beam response solution offsets.
 */
template <bool ApplyBeam, bool HasH5Parm>
void PrepareCorrectionMetadataForRow(const std::vector<MsGridder*>& gridders,
                                     const aocommon::MultiBandData& bands,
                                     std::vector<size_t>& time_offsets,
                                     VisibilityModifier& shared_modifier,
                                     MSProvider::MetaData& metadata, BatchRowData& rows) {
  size_t time_offset = time_offsets.back();
  // Cache the h5parm solutions
  for (MsGridder* gridder : gridders) {
    gridder->LoadCorrections<false, HasH5Parm>(bands, metadata.time, metadata.field_id,
                                               time_offset);
  }
  time_offsets.push_back(time_offset);
#ifdef HAVE_EVERYBEAM
  // Cache the offsets at which we will need to compute the beam response.
  // Don't yet actually compute the beam response.
  if constexpr (ApplyBeam) {
    bool has_time_update = shared_modifier.ShouldUpdateCachedBeamResponseForTime(metadata.time);
    int offset = 0;
    for (MsGridder* gridder : gridders) {
      offset =
          gridder->GetVisibilityModifier().ReserveCacheChunkSpaceForBeamResponse(has_time_update);
    }
    if (offset != -1) {
      rows.beam_cache_offsets.emplace_back(double(metadata.time), size_t(metadata.field_id),
                                           size_t(offset));
    }
  }
#endif
}

/**
 * Compute the beam responses previously pre-allocated in @ref PrepareCorrectionMetadataForRow.
 * This can be called in parallel and allow parallelism across facets however a strict ordering is
 * maintained on a per facet basis.
 */
template <bool ApplyBeam>
void ComputeBeamResponseForRowInOrder(const std::vector<MsGridder*>& gridders,
                                      const aocommon::MultiBandData& bands, BatchRowData& rows,
                                      CacheAlignedAtomic* index_process_order) {
#ifdef HAVE_EVERYBEAM
  if constexpr (ApplyBeam) {
    for (size_t gridder_index = 0; gridder_index < gridders.size(); ++gridder_index) {
      auto& atomic_value = index_process_order[gridder_index].value;
      size_t current = atomic_value.load();
      while (current != rows.index) {
        atomic_value.wait(current);
        current = atomic_value.load();
      }
      for (const auto& [meta_time, field_id, offset] : rows.beam_cache_offsets) {
        gridders[gridder_index]->GetVisibilityModifier().ComputeBeamResponseForCache(
            offset, meta_time, field_id, bands);
      }
      atomic_value++;
      atomic_value.notify_all();
    }
  }
#endif
}

/*
 * The input for this function is the predicted facet visibilities for a single
 * facet.
 * 1. Expand the visibilities.
 * 2. Apply corrections to the visibilities.
 * 3. Cumulatively add the expanded and corrected visibilities with those of
 * other facets.
 */
void ExpandAndCombineFacetVisibilities(
    size_t n_antennas, size_t n_channels, size_t n_rows, size_t n_vis_polarizations,
    const aocommon::MultiBandData& bands, const size_t* data_desc_ids, const size_t* antennas1,
    const size_t* antennas2, const size_t* field_ids, const double* times, const double* uvw,
    const std::complex<float>* facet_visibilities, MsGridder& gridder,
    std::complex<float>* combined_visibilities, CacheAlignedMutex* sum_visibilities_mutexes,
    size_t mutex_chunk_size) {
  aocommon::UVector<std::complex<float>> visibilities_scratch(n_channels * n_vis_polarizations);
  const size_t n_chunks = n_rows / mutex_chunk_size + 1;
  for (size_t chunk = 0; chunk < n_chunks; ++chunk) {
    std::lock_guard<std::mutex> sum_lock(sum_visibilities_mutexes[chunk].value);
    const size_t start = chunk * mutex_chunk_size;
    const size_t end = std::min(start + mutex_chunk_size, n_rows);
    for (size_t row = start; row < end; ++row) {
      switch (n_vis_polarizations) {
        case 1:
          std::copy_n(facet_visibilities, n_channels, visibilities_scratch.data());
          break;
        case 2:
          internal::ExpandData<2>(n_channels, facet_visibilities, visibilities_scratch.data(),
                                  gridder.Polarization());
          break;
        case 4:
          internal::ExpandData<4>(n_channels, facet_visibilities, visibilities_scratch.data(),
                                  gridder.Polarization());
          break;
      }
      gridder.CorrectInstrumentalVisibilities(n_antennas, bands, *data_desc_ids,
                                              visibilities_scratch.data(), uvw, *field_ids,
                                              *antennas1, *antennas2, *times);

      // In case the value was not sampled in this pass, it has been set to
      // infinite and should not overwrite the current value in the set.
      // NB! This is only true for multi pass gridders (wstacking gridder) other
      // gridders could consider skipping this check if it becomes important for
      // performance.
      for (size_t channel = 0; channel < n_channels * n_vis_polarizations; ++channel) {
        if (std::isfinite(visibilities_scratch[channel].real())) {
          combined_visibilities[channel] += visibilities_scratch[channel];
        }
      }
      facet_visibilities += n_channels;
      combined_visibilities += n_channels * n_vis_polarizations;
      data_desc_ids++;
      field_ids++;
      antennas1++;
      antennas2++;
      times++;
    }
  }
}

/**
 * Perform predict on a single block of data stored in @ref
 * PredictionChunkData, for a single gridder/facet.
 */
void PredictChunkForFacet(MsGridder& gridder, size_t facet_index,
                          const PredictionChunkData& chunk_data, size_t n_vis_polarizations,
                          size_t n_antennas,
                          std::vector<std::complex<float>>& combined_visibilities,
                          const aocommon::UVector<double>& frequencies,
                          MsProviderCollection::MsData& ms_data,
                          CacheAlignedMutex* sum_visibilities_mutexes, size_t mutex_chunk_size) {
  Logger::Info << "Predicting facet " + std::to_string(facet_index) + "\n";

  // TODO
  if (!ms_data.ms_provider->IsRegular())
    throw std::runtime_error("Can't use bulk inversion with BDA yet");
  const size_t n_channels = frequencies.size();

  aocommon::UVector<std::complex<float>> facet_visibilities(chunk_data.n_rows * n_channels);
  gridder.PredictChunk(chunk_data.n_rows, n_channels, frequencies.data(), chunk_data.uvws.data(),
                       facet_visibilities.data());

  ExpandAndCombineFacetVisibilities(
      n_antennas, n_channels, chunk_data.n_rows, n_vis_polarizations,
      ms_data.ms_provider->SelectedBands(), chunk_data.data_desc_ids.data(),
      chunk_data.antennas1.data(), chunk_data.antennas2.data(), chunk_data.field_ids.data(),
      chunk_data.times.data(), chunk_data.uvws.data(), facet_visibilities.data(), gridder,
      combined_visibilities.data(), sum_visibilities_mutexes, mutex_chunk_size);

  Logger::Info << "Done predicting facet " + std::to_string(facet_index) + "\n";
}
}  // namespace

MSGridderManagerScheduler::MSGridderManagerScheduler(size_t n_workers)
    : overlapping_task_processor_(scheduler_task_queue_) {
  InitializeThreadPoolForTaskQueue(worker_thread_pool_, worker_task_queue_, n_workers);
  InitializeThreadPoolForTaskQueue(scheduler_thread_pool_, scheduler_task_queue_,
                                   kSchedulerTaskQueueSize_);
}

void MSGridderManagerScheduler::InitializeThreadPoolForTaskQueue(
    std::vector<std::thread>& thread_pool_, aocommon::TaskQueue<std::function<void()>>& task_queue_,
    size_t n_threads) {
  if (thread_pool_.size() == 0) {
    thread_pool_.reserve(n_threads);
    for (size_t i = 0; i < n_threads; ++i) {
      thread_pool_.emplace_back([&] {
        std::function<void()> operation;
        while (task_queue_.Pop(operation)) {
          operation();
        }
      });
    }
  }
}

MSGridderManagerScheduler::~MSGridderManagerScheduler() {
  worker_task_queue_.Finish();
  for (std::thread& thread : worker_thread_pool_) {
    thread.join();
  }
  scheduler_task_queue_.Finish();
  for (std::thread& thread : scheduler_thread_pool_) {
    thread.join();
  }
}

void MSGridderManager::InitializeMS(GriddingTask& task) {
  for (const MsListItem& item : task.msList) {
    ms_provider_collection_.Add(item.ms_description->GetProvider(),
                                item.ms_description->Selection(), item.ms_index);
  }
  ms_provider_collection_.InitializeMS();
}

void MSGridderManager::InitializeGridders(GriddingTask& task,
                                          const std::vector<size_t>& facet_indices,
                                          const Resources& resources,
                                          std::vector<GriddingResult::FacetData>& facet_results,
                                          GriddingTaskManager* writer_lock_manager) {
  available_memory_ = resources.Memory();
  available_cores_ = resources.NCpus();
  Resources per_gridder_resources = resources.GetPart(task.num_parallel_gridders_);
  available_cores_per_gridder_ = per_gridder_resources.NCpus();
  const bool is_first_polarization = task.polarization == *settings_.polarizations.begin();
  for (size_t facet_index : facet_indices) {
    assert(facet_index < task.facets.size());

    // Create a new gridder for each facet / sub-task, since gridders do not
    // support reusing them for multiple tasks.
    std::unique_ptr<MsGridder> gridder = ConstructGridder(per_gridder_resources);
    GriddingTask::FacetData* facet_task = &task.facets[facet_index];
    GriddingResult::FacetData* facet_result = &facet_results[facet_index];

    if (solution_data_.HasData()) {
      gridder->GetVisibilityModifier().SetH5Parm(
          solution_data_.GetH5Parms(), solution_data_.GetFirstSolutions(),
          solution_data_.GetSecondSolutions(), solution_data_.GetGainTypes());
    }
    InitializeGridderForTask(*gridder, task, writer_lock_manager);

    const bool has_input_average_beam(facet_task->averageBeam);
    if (has_input_average_beam) {
      assert(dynamic_cast<IdgMsGridder*>(gridder.get()));
      IdgMsGridder& idgGridder = static_cast<IdgMsGridder&>(*gridder);
      idgGridder.SetAverageBeam(std::move(facet_task->averageBeam));
    }

    InitializeGridderForFacet(is_first_polarization, *gridder, *facet_task);

    facet_tasks_.emplace_back(GriddingFacetTask{std::move(gridder), facet_task, facet_result});
  }
}

size_t MSGridderManager::ReadChunkForInvert(
    GainMode gain_mode, bool apply_corrections, const std::vector<MsGridder*>& gridders,
    MsProviderCollection::MsData& ms_data, size_t n_chunk_rows, MSReader& ms_reader,
    const bool* selected_buffer, InversionChunkData& chunk_data, MsGridderData& shared_data) {
  const size_t n_parms = gridders[0]->NumValuesPerSolution();
  switch (gain_mode) {
    case GainMode::kXX:
      return ReadChunkForInvertImplementation<GainMode::kXX>(
          n_parms, apply_corrections, gridders, ms_data, n_chunk_rows, ms_reader, selected_buffer,
          chunk_data, shared_data);
      break;
    case GainMode::kYY:
      return ReadChunkForInvertImplementation<GainMode::kYY>(
          n_parms, apply_corrections, gridders, ms_data, n_chunk_rows, ms_reader, selected_buffer,
          chunk_data, shared_data);
      break;
    case GainMode::k2VisDiagonal:
      return ReadChunkForInvertImplementation<GainMode::k2VisDiagonal>(
          n_parms, apply_corrections, gridders, ms_data, n_chunk_rows, ms_reader, selected_buffer,
          chunk_data, shared_data);
      break;
    case GainMode::kTrace:
      return ReadChunkForInvertImplementation<GainMode::kTrace>(
          n_parms, apply_corrections, gridders, ms_data, n_chunk_rows, ms_reader, selected_buffer,
          chunk_data, shared_data);
      break;
    case GainMode::kFull:
      return ReadChunkForInvertImplementation<GainMode::kFull>(
          n_parms, apply_corrections, gridders, ms_data, n_chunk_rows, ms_reader, selected_buffer,
          chunk_data, shared_data);
      break;
  }
  assert(false);
  return 0;
}

template <GainMode Mode>
size_t MSGridderManager::ReadChunkForInvertImplementation(
    size_t n_parms, bool apply_corrections, const std::vector<MsGridder*>& gridders,
    MsProviderCollection::MsData& ms_data, size_t n_chunk_rows, MSReader& ms_reader,
    const bool* selected_buffer, InversionChunkData& chunk_data, MsGridderData& shared_data) {
  if (n_parms == 2) {
    return ReadChunkForInvertImplementation<Mode, 2>(apply_corrections, gridders, ms_data,
                                                     n_chunk_rows, ms_reader, selected_buffer,
                                                     chunk_data, shared_data);
  } else {
    return ReadChunkForInvertImplementation<Mode, 4>(apply_corrections, gridders, ms_data,
                                                     n_chunk_rows, ms_reader, selected_buffer,
                                                     chunk_data, shared_data);
  }
}

template <GainMode Mode, size_t NParms>
size_t MSGridderManager::ReadChunkForInvertImplementation(
    bool apply_corrections, const std::vector<MsGridder*>& gridders,
    MsProviderCollection::MsData& ms_data, size_t n_chunk_rows, MSReader& ms_reader,
    const bool* selected_buffer, InversionChunkData& chunk_data, MsGridderData& shared_data) {
  if (apply_corrections) {
    return ReadChunkForInvertImplementation<Mode, NParms, true>(
        gridders, ms_data, n_chunk_rows, ms_reader, selected_buffer, chunk_data, shared_data);
  } else {
    return ReadChunkForInvertImplementation<Mode, NParms, false>(
        gridders, ms_data, n_chunk_rows, ms_reader, selected_buffer, chunk_data, shared_data);
  }
}

template <GainMode Mode, size_t NParms, bool ApplyCorrections>
size_t MSGridderManager::ReadChunkForInvertImplementation(const std::vector<MsGridder*>& gridders,
                                                          MsProviderCollection::MsData& ms_data,
                                                          size_t n_chunk_rows, MSReader& ms_reader,
                                                          const bool* selected_buffer,
                                                          InversionChunkData& chunk_data,
                                                          MsGridderData& shared_data) {
  if constexpr (ApplyCorrections) {
    const bool apply_beam = settings_.applyFacetBeam || settings_.gridWithBeam;
    if (apply_beam) {
      return ReadChunkForInvertImplementation<Mode, NParms, ApplyCorrections, true>(
          gridders, ms_data, n_chunk_rows, ms_reader, selected_buffer, chunk_data, shared_data);
    } else {
      return ReadChunkForInvertImplementation<Mode, NParms, ApplyCorrections, false>(
          gridders, ms_data, n_chunk_rows, ms_reader, selected_buffer, chunk_data, shared_data);
    }
  } else {
    return ReadChunkForInvertImplementation<Mode, NParms, ApplyCorrections, false, false, false>(
        gridders, ms_data, n_chunk_rows, ms_reader, selected_buffer, chunk_data, shared_data);
  }
}

template <GainMode Mode, size_t NParms, bool ApplyCorrections, bool ApplyBeam>
size_t MSGridderManager::ReadChunkForInvertImplementation(const std::vector<MsGridder*>& gridders,
                                                          MsProviderCollection::MsData& ms_data,
                                                          size_t n_chunk_rows, MSReader& ms_reader,
                                                          const bool* selected_buffer,
                                                          InversionChunkData& chunk_data,
                                                          MsGridderData& shared_data) {
  const bool apply_forward = gridders[0]->GetPsfMode() == PsfMode::kDirectionDependent;
  if (apply_forward) {
    return ReadChunkForInvertImplementation<Mode, NParms, ApplyCorrections, ApplyBeam, true>(
        gridders, ms_data, n_chunk_rows, ms_reader, selected_buffer, chunk_data, shared_data);
  } else {
    return ReadChunkForInvertImplementation<Mode, NParms, ApplyCorrections, ApplyBeam, false>(
        gridders, ms_data, n_chunk_rows, ms_reader, selected_buffer, chunk_data, shared_data);
  }
}

template <GainMode Mode, size_t NParms, bool ApplyCorrections, bool ApplyBeam, bool ApplyForward>
size_t MSGridderManager::ReadChunkForInvertImplementation(const std::vector<MsGridder*>& gridders,
                                                          MsProviderCollection::MsData& ms_data,
                                                          size_t n_chunk_rows, MSReader& ms_reader,
                                                          const bool* selected_buffer,
                                                          InversionChunkData& chunk_data,
                                                          MsGridderData& shared_data) {
  const bool has_h5_parm = gridders[0]->visibility_modifier_.HasH5Parm();
  if (has_h5_parm) {
    return ReadChunkForInvertImplementation<Mode, NParms, ApplyCorrections, ApplyBeam, ApplyForward,
                                            true>(gridders, ms_data, n_chunk_rows, ms_reader,
                                                  selected_buffer, chunk_data, shared_data);
  } else {
    return ReadChunkForInvertImplementation<Mode, NParms, ApplyCorrections, ApplyBeam, ApplyForward,
                                            false>(gridders, ms_data, n_chunk_rows, ms_reader,
                                                   selected_buffer, chunk_data, shared_data);
  }
}

template <GainMode Mode, size_t NParms, bool ApplyCorrections, bool ApplyBeam, bool ApplyForward,
          bool HasH5Parm>
size_t MSGridderManager::ReadChunkForInvertImplementation(const std::vector<MsGridder*>& gridders,
                                                          MsProviderCollection::MsData& ms_data,
                                                          size_t n_chunk_rows, MSReader& ms_reader,
                                                          const bool* selected_buffer,
                                                          InversionChunkData& chunk_data,
                                                          MsGridderData& shared_data) {
  size_t n_chunk_rows_read = 0;

  // TODO
  const aocommon::MultiBandData& bands = ms_reader.Provider().SelectedBands();
  assert(bands.BandCount() == 1);
  const aocommon::BandData& band = *bands.begin();

  const size_t row_size = band.ChannelCount() * ms_data.ms_provider->NPolarizations();
  const size_t model_stride = shared_data.DoSubtractModel() ? row_size : 0;
  const size_t uvws_stride = 3;
  const size_t weights_stride = row_size;
  const size_t visibilities_in_stride = row_size;
  size_t visibilities_out_stride = band.ChannelCount();
  // When not applying corrections we collapse the polarizations.
  if constexpr (ApplyCorrections) {
    visibilities_out_stride *= ms_data.ms_provider->NPolarizations();
  }

  // Allow reading to get a bit ahead of processing but not by too much.
  aocommon::Lane<BatchRowData> task_lane(available_cores_ * 2);

  // Calculate maximum timesteps for beam cache allocation.
  // A PreallocatedVector is used to prevent race conditions, and it is over allocated (@sa
  // StartProcessingChunk()) as min(n_chunk_rows, unique_times), because timesteps may be unevenly
  // distributed across chunks, and computing the exact count requires iterating the full rows.
  const size_t n_max_beam_time_steps = std::min(n_chunk_rows, ms_data.unique_times->size());

  // Visibility modifier needs to know when we are starting a new chunk in order
  // to track the beam correctly.
  for (MsGridder* gridder : gridders) {
    gridder->GetVisibilityModifier().StartProcessingChunk(n_chunk_rows, n_max_beam_time_steps);
  }

  // Pre-process rows serially and partition into chunks for parallel execution.
  // The serial portion is kept minimal to avoid bottlenecks, with chunk sizes
  // balanced to offset thread overhead against worker throughput.
  // Current sizing is an educated heuristic; it may require further tuning or
  // dynamic calculation based on row length and other variable inputs.
  constexpr size_t kProcessingChunkTargetSizeRows = 4096;
  // This thread performs the following:
  //     1. Continues reading new data into the pipeline.
  //     2. Executes strictly ordered pre-processing (kept to an absolute minimum).
  //     3. Offloads prepared data to the task lane for parallelized completion.
  std::thread read_rows_thread([&] {
    std::pair<size_t, size_t>* antennas = chunk_data.antennas.data();
    size_t uvws_offset = 0;
    size_t visibilities_offset = 0;
    size_t time_offsets_offset = 0;
    size_t read_index = 0;

    while (n_chunk_rows_read < n_chunk_rows && ms_reader.CurrentRowAvailable()) {
      BatchRowData rows(kProcessingChunkTargetSizeRows, row_size);
      rows.uvws_out_offset = uvws_offset;
      rows.visibilities_out_offset = visibilities_offset;
      rows.time_offsets_out_offset = time_offsets_offset;
      rows.index = read_index++;
      BatchRowView row_view(&rows, model_stride, weights_stride, visibilities_in_stride,
                            uvws_stride, kProcessingChunkTargetSizeRows);
      while (n_chunk_rows_read < n_chunk_rows && row_view.ReadRow(ms_reader, shared_data)) {
        MSProvider::MetaData& metadata = *row_view.metadata_;
        if constexpr (ApplyCorrections) {
          // Cache the baseline.
          *antennas = std::make_pair(metadata.antenna1, metadata.antenna2);
          // Cache h5parm response and time offset.
          // Cache beam response offsets and pre-allocate space for beam response computation.
          PrepareCorrectionMetadataForRow<ApplyBeam, HasH5Parm>(
              gridders, bands, chunk_data.time_offsets, shared_data.GetVisibilityModifier(),
              metadata, rows);
          ++antennas;
          ++time_offsets_offset;
        }
        ms_reader.NextInputRow();
        row_view.NextRow();
        ++n_chunk_rows_read;
        uvws_offset += uvws_stride;
        visibilities_offset += visibilities_out_stride;
      }
      task_lane.write(std::move(rows));
    }
    task_lane.write_end();
  });

  // Per gridder mutex to prevent data race on correction sums inside
  // ApplyCorrections()
  std::unique_ptr<CacheAlignedMutex[]> gridder_mutexes =
      std::make_unique_for_overwrite<CacheAlignedMutex[]>(gridders.size());
  std::unique_ptr<CacheAlignedAtomic[]> index_process_order =
      std::make_unique_for_overwrite<CacheAlignedAtomic[]>(gridders.size());
  // NB! This deliberately leads to overallocation of threads
  // As this still outperforms the alternative of not overlapping the IO.
  // Future changes should implement task stealing which would
  // fix this overallocation.
  std::vector<std::thread> thread_pool_process;
  thread_pool_process.reserve(available_cores_);
  for (size_t i = 0; i < available_cores_; ++i) {
    thread_pool_process.emplace_back([&] {
      aocommon::UVector<float> image_weights_temp(band.ChannelCount());

      // Per thread local counters to prevent contention and race conditions.
      // Propagate to global counters once per thread at end of processing loop.
      size_t local_gridded_visibility_count = 0;
      double local_total_weight = 0.0;
      double local_max_gridded_weight = 0.0;
      double local_visibility_weight_sum = 0.0;
      BatchRowData rows;
      while (task_lane.read(rows)) {
        // For this part of processing ordering still matters.
        // index_process_order is used per facet to enforce ordering while
        // allowing some parallelism.
        ComputeBeamResponseForRowInOrder<ApplyBeam>(gridders, bands, rows,
                                                    index_process_order.get());

        // For the remaining processing batches can be processed out of order.
        // Map/copy from row view to chunk view.
        InversionChunkView<ApplyCorrections> chunk_view_out(&chunk_data, rows, uvws_stride,
                                                            visibilities_out_stride);
        BatchRowView row_view_in(&rows, model_stride, weights_stride, visibilities_in_stride,
                                 uvws_stride, rows.n_rows);
        while (!row_view_in.End()) {
          chunk_view_out.uvws_[0] = row_view_in.metadata_->u_in_m;
          chunk_view_out.uvws_[1] = row_view_in.metadata_->v_in_m;
          chunk_view_out.uvws_[2] = row_view_in.metadata_->w_in_m;

          shared_data.ModifyVisibilities<true>(row_view_in.visibilities_, chunk_view_out.uvws_,
                                               band, row_view_in.model_);
          shared_data.CalculateWeights(row_view_in.weights_, image_weights_temp.data(),
                                       chunk_view_out.uvws_, band, selected_buffer);

          // Sum the corrections and apply the weights.
          // We store the appropriate time_offset to be used later along with
          // other required info when we apply the corrections
          if constexpr (ApplyCorrections) {
            size_t& time_offset = *chunk_view_out.time_offsets_;
            for (size_t gridder_index = 0; gridder_index < gridders.size(); ++gridder_index) {
              std::lock_guard<std::mutex> lock(gridder_mutexes[gridder_index].value);
              const BeamResponseMap empty_beam_response;
              const BeamResponseMap& beam_response = [&]() -> const BeamResponseMap& {
#ifdef HAVE_EVERYBEAM
                if constexpr (ApplyBeam) {
                  return gridders[gridder_index]->GetVisibilityModifier().GetCachedBeamResponseMap(
                      chunk_data.index, row_view_in.current_chunk_row_);
                }
#endif
                return empty_beam_response;
              }();
              gridders[gridder_index]
                  ->ApplyCorrections<Mode, NParms, ModifierBehaviour::kSum, ApplyBeam, ApplyForward,
                                     HasH5Parm>(
                      ms_data.antenna_names.size(), row_view_in.visibilities_, band,
                      rows.data_desc_id, row_view_in.weights_, row_view_in.metadata_->antenna1,
                      row_view_in.metadata_->antenna2, time_offset, image_weights_temp.data(),
                      beam_response);
            }
          }

          shared_data.ApplyWeights<Mode>(row_view_in.visibilities_, band.ChannelCount(),
                                         row_view_in.weights_, image_weights_temp.data(),
                                         local_gridded_visibility_count, local_total_weight,
                                         local_max_gridded_weight, local_visibility_weight_sum);

          // When not applying corrections we collapse the polarizations.
          // When applying correction we need to keep them.
          if constexpr (!ApplyCorrections) {
            if (ms_data.ms_provider->NPolarizations() == 2) {
              internal::CollapseData<2>(band.ChannelCount(), row_view_in.visibilities_,
                                        shared_data.Polarization());
            } else if (ms_data.ms_provider->NPolarizations() == 4) {
              internal::CollapseData<4>(band.ChannelCount(), row_view_in.visibilities_,
                                        shared_data.Polarization());
            }
          }
          std::copy_n(row_view_in.visibilities_, visibilities_out_stride,
                      chunk_view_out.visibilities_);
          chunk_view_out.NextRow();
          row_view_in.NextRow();
        }
      }
      shared_data.AddVisibilityCounts(local_gridded_visibility_count, local_total_weight,
                                      local_max_gridded_weight, local_visibility_weight_sum);
    });
  }
  read_rows_thread.join();
  for (std::thread& thread : thread_pool_process) {
    thread.join();
  }

  // Visibility modifier needs to know when we have finished a chunk in order to
  // track the beam correctly.
  for (MsGridder* gridder : gridders) {
    gridder->GetVisibilityModifier().FinishProcessingChunk();
  }

  return n_chunk_rows_read;
}

void MSGridderManager::Invert() {
  InitializeMSDataVectors();

  for (const GriddingFacetTask& task : facet_tasks_) {
    const std::unique_ptr<MsGridder>& gridder = task.facet_gridder;
    gridder->CalculateOverallMetaData();
    gridder->StartInversion();
    const size_t n_inversion_passes = gridder->GetNInversionPasses();
    for (size_t pass_index = 0; pass_index < n_inversion_passes; ++pass_index) {
      gridder->StartInversionPass(pass_index);
      for (MsProviderCollection::MsData& ms_data : ms_provider_collection_.ms_data_vector_) {
        gridder->StartMeasurementSet(ms_provider_collection_.Count(), ms_data, false);
        ms_data.total_rows_processed += gridder->GridMeasurementSet(ms_data);
      }
      gridder->FinishInversionPass(pass_index);
    }
    gridder->FinishInversion();
  }
}

size_t MSGridderManager::GridChunk(bool apply_corrections, size_t n_vis_polarizations,
                                   const aocommon::MultiBandData& bands,
                                   const InversionChunkData& chunk_data,
                                   const aocommon::UVector<double>& frequencies,
                                   const MsProviderCollection::MsData& ms_data, size_t chunk_index,
                                   std::function<void()> signal_first_gridder_has_started,
                                   std::function<void()> signal_last_gridder_has_started) {
  Logger::Info << "Gridding " + std::to_string(chunk_data.n_rows) + " rows for " +
                      std::to_string(facet_tasks_.size()) + " facets using " +
                      std::to_string(available_cores_) + " threads " +
                      std::to_string(available_cores_per_gridder_) + " threads per gridder...\n";

  // TODO
  if (bands.BandCount() != 1) throw std::runtime_error("Can't use bulk inversion with BDA yet");
  const aocommon::BandData& band = *bands.begin();
  const size_t data_desc_id = *bands.DataDescIds().begin();

  ExecuteForAllGriddersWithNCores(
      available_cores_per_gridder_,
      [&](MsGridder* gridder, size_t facet_index) {
        GridChunkForFacet(*gridder, facet_index, apply_corrections, n_vis_polarizations, band,
                          data_desc_id, chunk_data, frequencies, ms_data, chunk_index);
      },
      std::move(signal_first_gridder_has_started), std::move(signal_last_gridder_has_started));
  Logger::Info << "Finished gridding " + std::to_string(chunk_data.n_rows) + " rows for " +
                      std::to_string(facet_tasks_.size()) + " facets.\n";
  return chunk_data.n_rows * facet_tasks_.size();
}

void MSGridderManager::GridChunkForFacet(MsGridder& gridder, size_t facet_index,
                                         bool apply_corrections, size_t n_vis_polarizations,
                                         const aocommon::BandData& band, size_t data_desc_id,
                                         const InversionChunkData& chunk_data,
                                         const aocommon::UVector<double>& frequencies,
                                         const MsProviderCollection::MsData& ms_data,
                                         size_t chunk_index) {
  Logger::Info << "Gridding facet " + std::to_string(facet_index) + "\n";

  const std::vector<std::complex<float>>& parm_response =
      gridder.GetVisibilityModifier().GetCachedParmResponse(
          ms_data.original_ms_index)[data_desc_id];

  std::shared_ptr<BeamResponseCacheChunk> beam_response =
      gridder.GetVisibilityModifier().TakeCachedBeamResponseChunk(chunk_index);

  gridder.gridded_visibility_count_ = chunk_data.gridded_visibility_count;
  gridder.visibility_weight_sum_ = chunk_data.visibility_weight_sum;
  gridder.max_gridded_weight_ = chunk_data.max_gridded_weight;
  gridder.total_weight_ = chunk_data.total_weight;

  gridder.GridSharedMeasurementSetChunk(
      apply_corrections, n_vis_polarizations, chunk_data.n_rows, chunk_data.uvw.data(),
      frequencies.data(), band, data_desc_id, chunk_data.antennas.data(),
      chunk_data.visibilities.data(),
      apply_corrections ? chunk_data.time_offsets.data() + 1 : nullptr,
      ms_data.antenna_names.size(), parm_response, *beam_response.get());

  Logger::Info << "Done gridding facet " + std::to_string(facet_index) + "\n";
}

void MSGridderManager::ReadChunksForInvert(aocommon::Lane<InversionChunkData>& task_lane,
                                           size_t n_max_rows_in_memory, bool apply_corrections,
                                           MsProviderCollection::MsData& ms_data,
                                           MsGridderData& shared_data,
                                           const std::vector<MsGridder*>& gridders,
                                           size_t n_vis_polarizations, const bool* selected_buffer,
                                           size_t& n_total_chunks) {
  // TODO
  const aocommon::MultiBandData& bands = ms_data.ms_provider->SelectedBands();
  assert(bands.BandCount() == 1);
  const aocommon::BandData& band = *bands.begin();

  // Processing more data per call minimizes total FFTs; aim for maximum memory utilization.
  Logger::Info << "Max " << n_max_rows_in_memory << " rows fit in memory.\n";
  // Enable pipelining by overlapping I/O and processing across multiple chunks.
  // The first chunk is smaller to reduce initial latency; currently using an
  // educated heuristic of 1.5, to be refined based on future benchmarks.
  const std::vector<size_t> chunks_to_process =
      CalculateChunkSizes(ms_data.ms_provider->NRows(), n_max_rows_in_memory, 1.5);
  n_total_chunks = chunks_to_process.size();
  Logger::Info << "Reading " + std::to_string(n_total_chunks) + " chunk(s) with sizes: {" +
                      NumberList::ToString(chunks_to_process) + "}\n";
  std::unique_ptr<MSReader> ms_reader = ms_data.ms_provider->MakeReader();

  size_t chunk_index = 0;
  while (ms_reader->CurrentRowAvailable()) {
    const size_t target_chunk_size = chunks_to_process[chunk_index];
    Logger::Info << "Loading " + std::to_string(target_chunk_size) + " rows into memory chunk " +
                        std::to_string(chunk_index + 1) + ".\n";

    InversionChunkData chunk_data(target_chunk_size, band.ChannelCount(), n_vis_polarizations,
                                  apply_corrections, chunk_index);
    if (apply_corrections) {
      chunk_data.time_offsets.reserve(target_chunk_size + 1);
      chunk_data.time_offsets.push_back(0);
    }

    const size_t n_rows =
        ReadChunkForInvert(shared_data.GetGainMode(), apply_corrections, gridders, ms_data,
                           target_chunk_size, *ms_reader, selected_buffer, chunk_data, shared_data);

    chunk_data.gridded_visibility_count = shared_data.gridded_visibility_count_;
    chunk_data.visibility_weight_sum = shared_data.visibility_weight_sum_;
    chunk_data.max_gridded_weight = shared_data.max_gridded_weight_;
    chunk_data.total_weight = shared_data.total_weight_;
    chunk_data.n_rows = n_rows;
    Logger::Debug << "Done loading chunk " + std::to_string(chunk_index + 1) + ".\n";

    task_lane.write(std::move(chunk_data));
    ++chunk_index;
  }
  // Visibility modifier needs to know when we are done processing chunks, so it
  // can free caches.
  for (const auto& gridder : gridders) {
    gridder->GetVisibilityModifier().FinishChunkedProcessing();
  };
  Logger::Info << "All gridding rows loaded.\n";
  task_lane.write_end();
}

void MSGridderManager::GridChunks(
    aocommon::Lane<InversionChunkData>& task_lane, bool apply_corrections,
    const aocommon::UVector<double>& frequencies, const aocommon::MultiBandData& bands,
    MsProviderCollection::MsData& ms_data, size_t n_vis_polarizations, const size_t& n_chunks,
    std::function<void()> signal_last_gridder_of_last_chunk_has_started) {
  // Process up to two chunks at a time in parallel.
  // If data for a second predict chunk becomes available while the first is
  // still processing then the second predict can make use of cores that would
  // otherwise be idle when the last few facets of the first predict are
  // finishing up.
  scheduler_->GetOverlappingTaskProcessor().Process<InversionChunkData>(
      task_lane,
      [&](InversionChunkData&& chunk_data, size_t chunk_index,
          std::binary_semaphore& first_gridder_started) mutable {
        Logger::Debug << "Gridding chunk " + std::to_string(chunk_index + 1) + "/" +
                             std::to_string(n_chunks) + ".\n";
        const bool is_last_chunk = chunk_index == n_chunks - 1;
        ms_data.total_rows_processed += GridChunk(
            apply_corrections, n_vis_polarizations, bands, chunk_data, frequencies, ms_data,
            chunk_index, [&]() { first_gridder_started.release(); },
            is_last_chunk ? signal_last_gridder_of_last_chunk_has_started : nullptr);
        Logger::Debug << "Done gridding chunk " + std::to_string(chunk_index + 1) + "/" +
                             std::to_string(n_chunks) + ".\n";
      },
      "Gridding");
}

void MSGridderManager::BatchInvert(std::function<void()> signal_last_gridding_work_has_started) {
  assert(facet_tasks_.size() > 1);
  assert(scheduler_);
  InitializeMSDataVectors();

  MsProviderCollection& providers = ms_provider_collection_;

  std::vector<MsGridder*> gridders;
  gridders.reserve(facet_tasks_.size());
  for (const GriddingFacetTask& task : facet_tasks_) {
    const std::unique_ptr<MsGridder>& gridder = task.facet_gridder;
    gridders.emplace_back(gridder.get());
  }

  ExecuteForAllGridders([=](MsGridder* gridder) {
    gridder->CalculateOverallMetaData();
    gridder->StartInversion();
  });
  const size_t n_inversion_passes = gridders[0]->GetNInversionPasses();
  for (size_t pass_index = 0; pass_index < n_inversion_passes; ++pass_index) {
    ExecuteForAllGridders([=](MsGridder* gridder) { gridder->StartInversionPass(pass_index); });
    size_t ms_index = 0;
    for (MsProviderCollection::MsData& ms_data : providers.ms_data_vector_) {
      const bool is_last_ms_data = (ms_index == providers.ms_data_vector_.size() - 1);
      Logger::Debug << "Starting gridding of MS " + std::to_string(ms_index + 1) + "/" +
                           std::to_string(providers.ms_data_vector_.size()) + "\n";

      MsGridderData shared_data(settings_);
      shared_data.CopyTaskData((*gridders[0]), solution_data_, ms_data);

      ExecuteForAllGridders(
          [&](MsGridder* gridder) {
            gridder->StartMeasurementSet(providers.Count(), ms_data, false);
          },
          false);
      shared_data.StartMeasurementSet(providers.Count(), ms_data, false);
      scheduler_->GetWorkerTaskQueue().WaitForIdle(available_cores_);

      const aocommon::MultiBandData& bands(ms_data.ms_provider->SelectedBands());
      const size_t n_max_channels = bands.MaxBandChannels();
      const size_t n_vis_polarizations = ms_data.ms_provider->NPolarizations();

      // We need to sum constant memory usage up across all gridders as each
      // gridder has its own internal memory usage based on image size
      size_t constant_mem = 0;
      for (MsGridder* gridder : gridders) {
        constant_mem += gridder->CalculateConstantMemory();
      }
      // We incur these additional per row memory overheads with data that we
      // have to cache for later in order to apply the corrections
      size_t additional_per_vis_row_mem = 0;
      bool apply_corrections = gridders[0]->WillApplyCorrections();
      if (apply_corrections) {
        // For each row we have to store an antenna pair and a solution time
        // offset
        additional_per_vis_row_mem = sizeof(size_t) * 3;
      }
      // Per visibility memory is only calculated once as its shared across
      // gridders.
      const size_t per_row_uvw_memory_consumption = sizeof(double) * 3;
      const size_t n_max_rows_in_memory = gridders[0]->CalculateMaxRowsInMemory(
          available_memory_, constant_mem, additional_per_vis_row_mem,
          per_row_uvw_memory_consumption, n_max_channels,
          apply_corrections ? n_vis_polarizations : 1);

      // TODO
      if (bands.BandCount() != 1) throw std::runtime_error("Can't use bulk inversion with BDA yet");
      const aocommon::BandData& band = *bands.begin();

      aocommon::UVector<double> frequencies(n_max_channels);
      for (size_t i = 0; i != frequencies.size(); ++i) {
        frequencies[i] = band.ChannelFrequency(i);
      }
      aocommon::UVector<bool> selected_buffer(n_max_channels, true);

      // Iterate over data in chunks until all visibilities
      // have been gridded.
      aocommon::Lane<InversionChunkData> task_lane(1);
      size_t n_total_chunks = 0;
      std::thread read_chunks_thread([&]() {
        ReadChunksForInvert(task_lane, n_max_rows_in_memory, apply_corrections, ms_data,
                            shared_data, gridders, n_vis_polarizations, selected_buffer.data(),
                            n_total_chunks);
      });
      // Only signal completion on the last MS to avoid premature workflow continuation
      GridChunks(task_lane, apply_corrections, frequencies, bands, ms_data, n_vis_polarizations,
                 n_total_chunks, is_last_ms_data ? signal_last_gridding_work_has_started : nullptr);
      read_chunks_thread.join();

      Logger::Debug << "Finished gridding of MS " + std::to_string(ms_index + 1) + "/" +
                           std::to_string(providers.ms_data_vector_.size()) + "\n";
      ++ms_index;
    }
    ExecuteForAllGridders([=](MsGridder* gridder) { gridder->FinishInversionPass(pass_index); });
  }
  ExecuteForAllGridders([](MsGridder* gridder) { gridder->FinishInversion(); });
}

void MSGridderManager::Predict() {
  InitializeMSDataVectors();

  for (const GriddingFacetTask& task : facet_tasks_) {
    const std::unique_ptr<MsGridder>& gridder = task.facet_gridder;
    gridder->CalculateOverallMetaData();
    gridder->StartPredict(std::move(task.facet_task->modelImages));
    const size_t n_predict_passes = gridder->GetNPredictPasses();
    for (size_t pass_index = 0; pass_index < n_predict_passes; ++pass_index) {
      gridder->StartPredictPass(pass_index);
      for (MsProviderCollection::MsData& ms_data : ms_provider_collection_.ms_data_vector_) {
        gridder->StartMeasurementSet(ms_provider_collection_.Count(), ms_data, true);
        ms_data.total_rows_processed += gridder->PredictMeasurementSet(ms_data);
      }
      gridder->FinishPredictPass(pass_index);
    }
    gridder->FinishPredict();
  }
}

size_t MSGridderManager::PredictChunk(const PredictionChunkData& chunk_data,
                                      size_t n_vis_polarizations, size_t n_antennas,
                                      std::vector<std::complex<float>>& combined_visibilities,
                                      const aocommon::UVector<double>& frequencies,
                                      MsProviderCollection::MsData& ms_data,
                                      std::function<void()> signal_first_predict_has_started,
                                      std::function<void()> signal_last_predict_has_started) {
  Logger::Info << "Predicting " + std::to_string(chunk_data.n_rows) + " rows for " +
                      std::to_string(facet_tasks_.size()) + " facets using " +
                      std::to_string(available_cores_) + " threads " +
                      std::to_string(available_cores_per_gridder_) + " threads per gridder...\n";

  // As all facets are summing their values into combined_visibilities we have
  // to protect it with a mutex to avoid data corruption.
  // Chunk into ranges and use a mutex per range instead of a single mutex to
  // strike a balance between gridder lock contention and time spent waiting.
  constexpr size_t kMutexChunkSize = 2048;
  std::unique_ptr<CacheAlignedMutex[]> sum_visibilities_mutexes =
      std::make_unique_for_overwrite<CacheAlignedMutex[]>((chunk_data.n_rows / kMutexChunkSize) +
                                                          1);

  ExecuteForAllGriddersWithNCores(
      available_cores_per_gridder_,
      [&](MsGridder* gridder, size_t facet_index) {
        PredictChunkForFacet(*gridder, facet_index, chunk_data, n_vis_polarizations, n_antennas,
                             combined_visibilities, frequencies, ms_data,
                             sum_visibilities_mutexes.get(), kMutexChunkSize);
      },
      std::move(signal_first_predict_has_started), std::move(signal_last_predict_has_started));

  Logger::Info << "Finished Predicting " + std::to_string(chunk_data.n_rows) + " rows for " +
                      std::to_string(facet_tasks_.size()) + " facets.\n";

  return chunk_data.n_rows * facet_tasks_.size();
}

void MSGridderManager::PredictChunks(
    aocommon::Lane<PredictionChunkData>& task_lane, const aocommon::UVector<double>& frequencies,
    MsProviderCollection::MsData& ms_data, size_t n_vis_polarizations, bool add_assign_model,
    const size_t& n_chunks, std::function<void()> signal_last_predict_of_last_chunk_has_started) {
  // TODO
  if (!ms_data.ms_provider->IsRegular())
    throw std::runtime_error("Can't use bulk inversion with BDA yet");
  const aocommon::BandData& band = *ms_data.ms_provider->SelectedBands().begin();

  // Process up to two chunks at a time in parallel.
  // If data for a second predict chunk becomes available while the first is
  // still processing then the second predict can make use of cores that would
  // otherwise be idle when the last few facets of the first predict are
  // finishing up.
  scheduler_->GetOverlappingTaskProcessor().Process<PredictionChunkData>(
      task_lane,
      [&](PredictionChunkData&& chunk_data, size_t chunk_index,
          std::binary_semaphore& first_predict_started) mutable {
        std::vector<std::complex<float>> combined_visibilities(
            chunk_data.n_rows * band.ChannelCount() * n_vis_polarizations);
        // Predict the per facet visibilities; expand and apply corrections
        // then combine them together.
        Logger::Debug << "Predicting chunk " + std::to_string(chunk_index + 1) + "/" +
                             std::to_string(n_chunks) + ".\n";
        bool is_last_chunk = chunk_index == n_chunks - 1;
        ms_data.total_rows_processed += PredictChunk(
            chunk_data, n_vis_polarizations, ms_data.antenna_names.size(), combined_visibilities,
            frequencies, ms_data, [&]() { first_predict_started.release(); },
            is_last_chunk ? signal_last_predict_of_last_chunk_has_started : nullptr);
        Logger::Debug << "Done predicting chunk " + std::to_string(chunk_index + 1) + "/" +
                             std::to_string(n_chunks) + ".\n";

        // Do a single write for the combined/expanded chunk of visibilities
        // of all facets.
        Logger::Debug << "Writing predicted chunk " + std::to_string(chunk_index + 1) + "/" +
                             std::to_string(n_chunks) + ".\n";
        std::complex<float>* visibilities = combined_visibilities.data();
        const size_t stride = band.ChannelCount() * n_vis_polarizations;
        for (size_t row = 0; row != chunk_data.n_rows; ++row) {
          ms_data.ms_provider->WriteModel(visibilities, add_assign_model);
          ms_data.ms_provider->NextOutputRow();
          visibilities += stride;
        }
        Logger::Debug << "Done writing predicted chunk " + std::to_string(chunk_index + 1) + "/" +
                             std::to_string(n_chunks) + ".\n";
      },
      "Predict");
}

void MSGridderManager::ReadChunksForPredict(
    aocommon::Lane<PredictionChunkData>& task_lane, size_t n_max_rows_in_memory,
    MsProviderCollection::MsData& ms_data, MsGridderData& shared_data,
    const std::vector<MsGridder*>& gridders, const aocommon::BandData band,
    size_t n_vis_polarizations, const bool* selected_buffer, size_t& n_total_chunks) {
  // Processing more data per call minimizes total FFTs; aim for maximum memory utilization.
  Logger::Info << "Max " << n_max_rows_in_memory << " rows fit in memory.\n";
  // Enable pipelining by overlapping I/O and processing across multiple chunks.
  // The first chunk is smaller to reduce initial latency; currently using an
  // educated heuristic of 1.5, to be refined based on future benchmarks.
  const std::vector<size_t> chunks_to_process =
      CalculateChunkSizes(ms_data.ms_provider->NRows(), n_max_rows_in_memory, 1.5);
  n_total_chunks = chunks_to_process.size();
  Logger::Info << "Reading " + std::to_string(chunks_to_process.size()) +
                      " chunk(s) with sizes: {" + NumberList::ToString(chunks_to_process) + "}\n";

  // Set provider up for writing
  ms_data.ms_provider->ReopenRW();
  ms_data.ms_provider->ResetWritePosition();
  std::unique_ptr<MSReader> ms_reader = ms_data.ms_provider->MakeReader();
  size_t chunk_index = 0;
  while (ms_reader->CurrentRowAvailable()) {
    const size_t target_chunk_size = chunks_to_process[chunk_index];
    Logger::Info << "Loading " + std::to_string(target_chunk_size) + " rows into memory chunk " +
                        std::to_string(chunk_index + 1) + ".\n";
    PredictionChunkData chunk_data(target_chunk_size);
    while (ms_reader->CurrentRowAvailable() && chunk_data.n_rows < target_chunk_size) {
      MSProvider::MetaData metadata;
      shared_data.ReadPredictMetaData(metadata);
      chunk_data.uvws[chunk_data.n_rows * 3] = metadata.u_in_m;
      chunk_data.uvws[chunk_data.n_rows * 3 + 1] = metadata.v_in_m;
      chunk_data.uvws[chunk_data.n_rows * 3 + 2] = metadata.w_in_m;
      chunk_data.antennas1[chunk_data.n_rows] = metadata.antenna1;
      chunk_data.antennas2[chunk_data.n_rows] = metadata.antenna2;
      chunk_data.field_ids[chunk_data.n_rows] = metadata.field_id;
      chunk_data.data_desc_ids[chunk_data.n_rows] = metadata.data_desc_id;
      chunk_data.times[chunk_data.n_rows] = metadata.time;
      chunk_data.n_rows++;
      ms_reader->NextInputRow();
    }
    Logger::Debug << "Done loading chunk " + std::to_string(chunk_index + 1) + ".\n";

    task_lane.write(std::move(chunk_data));
    ++chunk_index;
  }
  Logger::Info << "All predict rows loaded.\n";
  task_lane.write_end();
}

void MSGridderManager::BatchPredict(std::function<void()> signal_last_predict_work_has_started) {
  assert(facet_tasks_.size() > 1);
  assert(scheduler_);
  InitializeMSDataVectors();

  MsProviderCollection& providers = ms_provider_collection_;

  std::vector<MsGridder*> gridders;
  gridders.reserve(facet_tasks_.size());
  for (const GriddingFacetTask& task : facet_tasks_) {
    const std::unique_ptr<MsGridder>& gridder = task.facet_gridder;
    gridders.emplace_back(gridder.get());
  }

  ExecuteForAllGridders([=](MsGridder* gridder, GriddingFacetTask& task) {
    gridder->CalculateOverallMetaData();
    gridder->StartPredict(std::move(task.facet_task->modelImages));
  });

  const size_t n_predict_passes = gridders[0]->GetNPredictPasses();
  for (size_t pass_index = 0; pass_index < n_predict_passes; ++pass_index) {
    ExecuteForAllGridders([=](MsGridder* gridder) { gridder->StartPredictPass(pass_index); });
    size_t ms_index = 0;
    for (MsProviderCollection::MsData& ms_data : providers.ms_data_vector_) {
      const bool is_last_ms_data = (ms_index == providers.ms_data_vector_.size() - 1);
      Logger::Debug << "Starting prediction of MS " + std::to_string(ms_index + 1) + "/" +
                           std::to_string(providers.ms_data_vector_.size()) + "\n";
      MsGridderData shared_data(settings_);
      shared_data.CopyTaskData((*gridders[0]), solution_data_, ms_data);
      ExecuteForAllGridders(
          [&](MsGridder* gridder) {
            gridder->StartMeasurementSet(providers.Count(), ms_data, true);
          },
          false);
      shared_data.StartMeasurementSet(providers.Count(), ms_data, true);
      scheduler_->GetWorkerTaskQueue().WaitForIdle(available_cores_);

      // TODO
      if (!ms_data.ms_provider->IsRegular())
        throw std::runtime_error("Can't use bulk inversion with BDA yet");
      const aocommon::BandData& band = *ms_data.ms_provider->SelectedBands().begin();

      const size_t n_channels = band.ChannelCount();
      const size_t n_vis_polarizations = ms_data.ms_provider->NPolarizations();

      // Each gridder has its own constant memory usage:
      //   * Memory to store the grid
      //   * A copy of the dirty image
      //   * Memory to store predicted visibilities.
      size_t constant_mem = 0;
      for (MsGridder* gridder : gridders) {
        constant_mem += gridder->CalculateConstantMemory();
      }
      // Additionally there is per row memory usage across all gridders:
      //   * One set of expanded/corrected visibilities
      //   * One set of metadata
      size_t shared_mem_per_row = sizeof(std::complex<float>) * n_channels * n_vis_polarizations;
      shared_mem_per_row += sizeof(double) * 3;  // uvw
      shared_mem_per_row += sizeof(double) * 3;  // time
      shared_mem_per_row += sizeof(size_t) * 3;  // antenna1, antenna2, field_id
      // Which we approximate to a per gridder overhead instead.
      // This will not be perfect but should be close enough.
      double additional_per_vis_row_mem = shared_mem_per_row / gridders.size();
      // uvw is already accounted for in shared memory so don't count it again.
      const size_t per_row_uvw_memory_consumption = 0;
      // Compute maximum rows for a single gridder based on the above
      // assumptions/approximations.
      const size_t n_max_overall_rows_in_memory = gridders[0]->CalculateMaxRowsInMemory(
          available_memory_, constant_mem, additional_per_vis_row_mem,
          per_row_uvw_memory_consumption, n_channels, 1);
      // Finally approximate that again to multiple gridders to get maximum rows
      // per gridder.
      const size_t n_max_rows_in_memory = n_max_overall_rows_in_memory / gridders.size();

      aocommon::UVector<double> frequencies(n_channels);
      for (size_t i = 0; i != frequencies.size(); ++i) {
        frequencies[i] = band.ChannelFrequency(i);
      }
      aocommon::UVector<bool> selected_buffer(n_channels, true);

      // Iterate over data in chunks until all visibilities have been predicted.
      aocommon::Lane<PredictionChunkData> task_lane(1);
      size_t n_total_chunks = 0;
      std::thread read_chunks_thread([&] {
        ReadChunksForPredict(task_lane, n_max_rows_in_memory, ms_data, shared_data, gridders, band,
                             n_vis_polarizations, selected_buffer.data(), n_total_chunks);
      });
      // Only signal completion on the last MS to avoid premature workflow continuation
      PredictChunks(task_lane, frequencies, ms_data, n_vis_polarizations,
                    gridders[0]->ShouldAddAssignModel(), n_total_chunks,
                    is_last_ms_data ? signal_last_predict_work_has_started : nullptr);
      read_chunks_thread.join();
      Logger::Debug << "Finished prediction of MS " + std::to_string(ms_index + 1) + "/" +
                           std::to_string(providers.ms_data_vector_.size()) + "\n";
      ++ms_index;
    }
    ExecuteForAllGridders([=](MsGridder* gridder) { gridder->FinishPredictPass(pass_index); });
  }
  ExecuteForAllGridders([](MsGridder* gridder) { gridder->FinishPredict(); });
}

void MSGridderManager::ProcessResults(std::mutex& result_mutex, GriddingResult& result,
                                      bool store_common_info) {
  for (auto& [gridder, facet_task, facet_result] : facet_tasks_) {
    // Add facet-specific result values to the result.
    facet_result->images = gridder->ResultImages();
    facet_result->actualWGridSize = gridder->ActualWGridSize();
    facet_result->averageCorrection = gridder->GetAverageCorrection();
    facet_result->averageBeamCorrection = gridder->GetAverageBeamCorrection();
    facet_result->cache = gridder->AcquireMetaDataCache();

    // The gridder resets visibility counters in each gridding invocation,
    // so they only contain the statistics of that invocation.
    facet_result->imageWeight = gridder->ImageWeight();
    facet_result->normalizationFactor = gridder->NormalizationFactor();
    facet_result->effectiveGriddedVisibilityCount = gridder->EffectiveGriddedVisibilityCount();
    {
      std::lock_guard<std::mutex> result_lock(result_mutex);
      result.griddedVisibilityCount = gridder->GriddedVisibilityCount();
      result.visibilityWeightSum = gridder->VisibilityWeightSum();
    }

    // If the average beam already exists on input, IDG will not recompute it,
    // so in that case there is no need to return the unchanged average beam.
    const bool has_input_average_beam(facet_task->averageBeam);
    IdgMsGridder* idgGridder = dynamic_cast<IdgMsGridder*>(gridder.get());
    if (idgGridder && !has_input_average_beam) {
      facet_result->averageBeam = idgGridder->ReleaseAverageBeam();
    }

    if (store_common_info) {
      // Store result values that are equal for all facets.
      result.startTime = ms_provider_collection_.StartTime();
      result.beamSize = gridder->BeamSize();
    }
  }
}

std::unique_ptr<MsGridder> MSGridderManager::ConstructGridder(const Resources& resources) {
  switch (settings_.gridderType) {
    case GridderType::IDG:
      return std::make_unique<IdgMsGridder>(settings_, resources, ms_provider_collection_);
    case GridderType::FacetIDG:
      return std::make_unique<FacetIdgMsGridder>(settings_, resources, ms_provider_collection_);
    case GridderType::WGridder:
      return std::make_unique<WGriddingMSGridder>(settings_, resources, ms_provider_collection_,
                                                  false);
    case GridderType::TunedWGridder:
      return std::make_unique<WGriddingMSGridder>(settings_, resources, ms_provider_collection_,
                                                  true);
    case GridderType::WTowers:
#ifdef BUILD_WTOWERS
      return std::make_unique<WTowersMsGridder>(settings_, resources, ms_provider_collection_);
#else
      throw std::runtime_error("w-towers gridder is not available");
#endif
    case GridderType::DirectFT:
      switch (settings_.directFTPrecision) {
        case DirectFTPrecision::Float:
          return std::make_unique<DirectMSGridder<float>>(settings_, resources,
                                                          ms_provider_collection_);
        case DirectFTPrecision::Double:
          return std::make_unique<DirectMSGridder<double>>(settings_, resources,
                                                           ms_provider_collection_);
        case DirectFTPrecision::LongDouble:
          return std::make_unique<DirectMSGridder<long double>>(settings_, resources,
                                                                ms_provider_collection_);
      }
      break;
    case GridderType::WStacking:
      return std::make_unique<WSMSGridder>(settings_, resources, ms_provider_collection_);
  }
  return {};
}

void MSGridderManager::InitializeGridderForTask(MsGridder& gridder, const GriddingTask& task,
                                                GriddingTaskManager* writer_lock_manager) {
  gridder.SetGridMode(settings_.gridMode);

  gridder.SetFacetGroupIndex(task.facetGroupIndex);
  gridder.SetImagePadding(settings_.imagePadding);
  gridder.SetPhaseCentreDec(task.observationInfo.phaseCentreDec);
  gridder.SetPhaseCentreRA(task.observationInfo.phaseCentreRA);

  if (settings_.hasShift) {
    double main_image_dl = 0.0;
    double main_image_dm = 0.0;
    aocommon::ImageCoordinates::RaDecToLM(
        settings_.shiftRA, settings_.shiftDec, task.observationInfo.phaseCentreRA,
        task.observationInfo.phaseCentreDec, main_image_dl, main_image_dm);
    gridder.SetMainImageDL(main_image_dl);
    gridder.SetMainImageDM(main_image_dm);
  }

  gridder.SetPolarization(task.polarization);
  gridder.SetIsComplex(task.polarization == aocommon::Polarization::XY ||
                       task.polarization == aocommon::Polarization::YX);
  gridder.SetIsFirstTask(task.isFirstTask);
  gridder.SetImageWeights(task.imageWeights.get());
  if (task.operation == GriddingTask::Invert) {
    if (task.imagePSF) {
      if (settings_.ddPsfGridWidth > 1 || settings_.ddPsfGridHeight > 1) {
        gridder.SetPsfMode(PsfMode::kDirectionDependent);
      } else {
        gridder.SetPsfMode(PsfMode::kSingle);
      }
    } else {
      gridder.SetPsfMode(PsfMode::kNone);
    }
    gridder.SetDoSubtractModel(task.subtractModel);
    gridder.SetStoreImagingWeights(task.storeImagingWeights);
  } else {
    gridder.SetWriterLockManager(writer_lock_manager);
  }
}

void MSGridderManager::InitializeGridderForFacet(bool is_first_polarization, MsGridder& gridder,
                                                 GriddingTask::FacetData& facet_task) {
  const schaapcommon::facets::Facet* facet = facet_task.facet.get();
  gridder.SetAddAssignModel(facet, is_first_polarization);
  gridder.SetIsFacet(facet);
  if (facet) {
    gridder.SetFacetIndex(facet_task.index);
    gridder.SetImageWidth(facet->GetUntrimmedBoundingBox().Width());
    gridder.SetImageHeight(facet->GetUntrimmedBoundingBox().Height());
    gridder.SetTrimSize(facet->GetTrimmedBoundingBox().Width(),
                        facet->GetTrimmedBoundingBox().Height());
    gridder.GetVisibilityModifier().SetFacetDirection(facet->RA(), facet->Dec());
  } else {
    gridder.SetImageWidth(settings_.paddedImageWidth);
    gridder.SetImageHeight(settings_.paddedImageHeight);
    gridder.SetTrimSize(settings_.trimmedImageWidth, settings_.trimmedImageHeight);
  }
  gridder.SetLShift(facet_task.l_shift);
  gridder.SetMShift(facet_task.m_shift);

  std::unique_ptr<MetaDataCache> cache = std::move(facet_task.cache);
  if (!cache) cache = std::make_unique<MetaDataCache>();
  gridder.SetMetaDataCache(std::move(cache));
}

}  // namespace wsclean
