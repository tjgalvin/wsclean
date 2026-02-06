#include "msgriddermanager.h"

#include <functional>
#include <mutex>
#include <vector>

#include <aocommon/logger.h>
#include <aocommon/taskqueue.h>
#include <aocommon/threadpool.h>
#include <aocommon/uvector.h>
#include <schaapcommon/facets/facet.h>
#include <schaapcommon/h5parm/h5parm.h>
#include <schaapcommon/h5parm/jonesparameters.h>
#include <schaapcommon/h5parm/soltab.h>

#include "directmsgridder.h"
#include "h5solutiondata.h"
#include "msgridder.h"
#include "msprovidercollection.h"
#include "wsmsgridder.h"

#include "../idg/averagebeam.h"
#include "../idg/idgmsgridder.h"
#include "../idg/facetidgmsgridder.h"
#include "../main/settings.h"
#include "../structures/resources.h"
#include "../wgridder/wgriddingmsgridder.h"
#include "../wtowers/wtowersmsgridder.h"

using aocommon::Logger;

namespace wsclean {

MSGridderManagerScheduler::MSGridderManagerScheduler(size_t n_workers)
    : overlapping_task_processor_(scheduler_task_queue_) {
  InitializeThreadPoolForTaskQueue(worker_thread_pool_, worker_task_queue_,
                                   n_workers);
  InitializeThreadPoolForTaskQueue(
      scheduler_thread_pool_, scheduler_task_queue_, kSchedulerTaskQueueSize_);
}

void MSGridderManagerScheduler::InitializeThreadPoolForTaskQueue(
    std::vector<std::thread>& thread_pool_,
    aocommon::TaskQueue<std::function<void()>>& task_queue_, size_t n_threads) {
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
                                item.ms_description->Selection(),
                                item.ms_index);
  }
  ms_provider_collection_.InitializeMS();
}

void MSGridderManager::InitializeGridders(
    GriddingTask& task, const std::vector<size_t>& facet_indices,
    const Resources& resources,
    std::vector<GriddingResult::FacetData>& facet_results,
    GriddingTaskManager* writer_lock_manager) {
  available_memory_ = resources.Memory();
  available_cores_ = resources.NCpus();
  Resources per_gridder_resources =
      resources.GetPart(task.num_parallel_gridders_);
  available_cores_per_gridder_ = per_gridder_resources.NCpus();
  const bool is_first_polarization =
      task.polarization == *settings_.polarizations.begin();
  for (size_t facet_index : facet_indices) {
    assert(facet_index < task.facets.size());

    // Create a new gridder for each facet / sub-task, since gridders do not
    // support reusing them for multiple tasks.
    std::unique_ptr<MsGridder> gridder =
        ConstructGridder(per_gridder_resources);
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

    facet_tasks_.emplace_back(
        GriddingFacetTask{std::move(gridder), facet_task, facet_result});
  }
}

size_t MSGridderManager::ReadChunkForInvert(
    GainMode gain_mode, bool apply_corrections,
    const std::vector<MsGridder*>& gridders,
    MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
    MSReader& ms_reader, const bool* selected_buffer,
    InversionChunkData& chunk_data, MsGridderData& shared_data) {
  const size_t n_parms = gridders[0]->NumValuesPerSolution();
  switch (gain_mode) {
    case GainMode::kXX:
      return ReadChunkForInvertImplementation<GainMode::kXX>(
          n_parms, apply_corrections, gridders, ms_data, n_chunk_rows,
          ms_reader, selected_buffer, chunk_data, shared_data);
      break;
    case GainMode::kYY:
      return ReadChunkForInvertImplementation<GainMode::kYY>(
          n_parms, apply_corrections, gridders, ms_data, n_chunk_rows,
          ms_reader, selected_buffer, chunk_data, shared_data);
      break;
    case GainMode::k2VisDiagonal:
      return ReadChunkForInvertImplementation<GainMode::k2VisDiagonal>(
          n_parms, apply_corrections, gridders, ms_data, n_chunk_rows,
          ms_reader, selected_buffer, chunk_data, shared_data);
      break;
    case GainMode::kTrace:
      return ReadChunkForInvertImplementation<GainMode::kTrace>(
          n_parms, apply_corrections, gridders, ms_data, n_chunk_rows,
          ms_reader, selected_buffer, chunk_data, shared_data);
      break;
    case GainMode::kFull:
      return ReadChunkForInvertImplementation<GainMode::kFull>(
          n_parms, apply_corrections, gridders, ms_data, n_chunk_rows,
          ms_reader, selected_buffer, chunk_data, shared_data);
      break;
  }
  assert(false);
  return 0;
}

template <GainMode Mode>
size_t MSGridderManager::ReadChunkForInvertImplementation(
    size_t n_parms, bool apply_corrections,
    const std::vector<MsGridder*>& gridders,
    MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
    MSReader& ms_reader, const bool* selected_buffer,
    InversionChunkData& chunk_data, MsGridderData& shared_data) {
  if (n_parms == 2) {
    return ReadChunkForInvertImplementation<Mode, 2>(
        apply_corrections, gridders, ms_data, n_chunk_rows, ms_reader,
        selected_buffer, chunk_data, shared_data);
  } else {
    return ReadChunkForInvertImplementation<Mode, 4>(
        apply_corrections, gridders, ms_data, n_chunk_rows, ms_reader,
        selected_buffer, chunk_data, shared_data);
  }
}

template <GainMode Mode, size_t NParms>
size_t MSGridderManager::ReadChunkForInvertImplementation(
    bool apply_corrections, const std::vector<MsGridder*>& gridders,
    MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
    MSReader& ms_reader, const bool* selected_buffer,
    InversionChunkData& chunk_data, MsGridderData& shared_data) {
  if (apply_corrections) {
    return ReadChunkForInvertImplementation<Mode, NParms, true>(
        gridders, ms_data, n_chunk_rows, ms_reader, selected_buffer, chunk_data,
        shared_data);
  } else {
    return ReadChunkForInvertImplementation<Mode, NParms, false>(
        gridders, ms_data, n_chunk_rows, ms_reader, selected_buffer, chunk_data,
        shared_data);
  }
}

template <GainMode Mode, size_t NParms, bool ApplyCorrections>
size_t MSGridderManager::ReadChunkForInvertImplementation(
    const std::vector<MsGridder*>& gridders,
    MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
    MSReader& ms_reader, const bool* selected_buffer,
    InversionChunkData& chunk_data, MsGridderData& shared_data) {
  if constexpr (ApplyCorrections) {
    const bool apply_beam = settings_.applyFacetBeam || settings_.gridWithBeam;
    if (apply_beam) {
      return ReadChunkForInvertImplementation<Mode, NParms, ApplyCorrections,
                                              true>(
          gridders, ms_data, n_chunk_rows, ms_reader, selected_buffer,
          chunk_data, shared_data);
    } else {
      return ReadChunkForInvertImplementation<Mode, NParms, ApplyCorrections,
                                              false>(
          gridders, ms_data, n_chunk_rows, ms_reader, selected_buffer,
          chunk_data, shared_data);
    }
  } else {
    return ReadChunkForInvertImplementation<Mode, NParms, ApplyCorrections,
                                            false, false, false>(
        gridders, ms_data, n_chunk_rows, ms_reader, selected_buffer, chunk_data,
        shared_data);
  }
}

template <GainMode Mode, size_t NParms, bool ApplyCorrections, bool ApplyBeam>
size_t MSGridderManager::ReadChunkForInvertImplementation(
    const std::vector<MsGridder*>& gridders,
    MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
    MSReader& ms_reader, const bool* selected_buffer,
    InversionChunkData& chunk_data, MsGridderData& shared_data) {
  const bool apply_forward =
      gridders[0]->GetPsfMode() == PsfMode::kDirectionDependent;
  if (apply_forward) {
    return ReadChunkForInvertImplementation<Mode, NParms, ApplyCorrections,
                                            ApplyBeam, true>(
        gridders, ms_data, n_chunk_rows, ms_reader, selected_buffer, chunk_data,
        shared_data);
  } else {
    return ReadChunkForInvertImplementation<Mode, NParms, ApplyCorrections,
                                            ApplyBeam, false>(
        gridders, ms_data, n_chunk_rows, ms_reader, selected_buffer, chunk_data,
        shared_data);
  }
}

template <GainMode Mode, size_t NParms, bool ApplyCorrections, bool ApplyBeam,
          bool ApplyForward>
size_t MSGridderManager::ReadChunkForInvertImplementation(
    const std::vector<MsGridder*>& gridders,
    MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
    MSReader& ms_reader, const bool* selected_buffer,
    InversionChunkData& chunk_data, MsGridderData& shared_data) {
  const bool has_h5_parm = gridders[0]->visibility_modifier_.HasH5Parm();
  if (has_h5_parm) {
    return ReadChunkForInvertImplementation<Mode, NParms, ApplyCorrections,
                                            ApplyBeam, ApplyForward, true>(
        gridders, ms_data, n_chunk_rows, ms_reader, selected_buffer, chunk_data,
        shared_data);
  } else {
    return ReadChunkForInvertImplementation<Mode, NParms, ApplyCorrections,
                                            ApplyBeam, ApplyForward, false>(
        gridders, ms_data, n_chunk_rows, ms_reader, selected_buffer, chunk_data,
        shared_data);
  }
}

template <GainMode Mode, size_t NParms, bool ApplyCorrections, bool ApplyBeam,
          bool ApplyForward, bool HasH5Parm>
size_t MSGridderManager::ReadChunkForInvertImplementation(
    const std::vector<MsGridder*>& gridders,
    MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
    MSReader& ms_reader, const bool* selected_buffer,
    InversionChunkData& chunk_data, MsGridderData& shared_data) {
  size_t n_chunk_rows_read = 0;

  // TODO
  const aocommon::MultiBandData& bands = ms_reader.Provider().SelectedBands();
  assert(bands.BandCount() == 1);
  const aocommon::BandData& band = *bands.begin();

  const size_t uvws_stride = 3;
  // When not applying corrections we collapse the polarizations.
  size_t visibilities_stride = band.ChannelCount();
  if constexpr (ApplyCorrections) {
    visibilities_stride *= ms_data.ms_provider->NPolarizations();
  }

  // Initialize the row buffer with the desired batch size.
  size_t n_buffer_rows = 1000;
  const size_t n_row_size =
      band.ChannelCount() * ms_data.ms_provider->NPolarizations();

  // Allow reading to get a bit ahead of processing but not by too much.
  aocommon::Lane<BatchRowData> task_lane(available_cores_ * 2);

  // Visibility modifier needs to know when we are starting a new chunk in order
  // to track the beam correctly.
  for (MsGridder* gridder : gridders) {
    gridder->GetVisibilityModifier().StartProcessingChunk();
  }

  std::thread read_rows_thread([&] {
    std::pair<size_t, size_t>* antennas = chunk_data.antennas.data();
    size_t uvws_offset = 0;
    size_t visibilities_offset = 0;
    size_t time_offsets_offset = 0;

    while (ms_reader.CurrentRowAvailable() &&
           n_chunk_rows_read < n_chunk_rows) {
      BatchRowData rows(n_buffer_rows, n_row_size);
      rows.uvws_offset = uvws_offset;
      rows.visibilities_offset = visibilities_offset;
      rows.time_offsets_offset = time_offsets_offset;
      while (ms_reader.CurrentRowAvailable() &&
             rows.n_rows_read < n_buffer_rows &&
             n_chunk_rows_read < n_chunk_rows) {
        auto& metadata = rows.metadata_[rows.n_rows_read];
        ms_reader.ReadMeta(metadata);

        // Read and store all visibilities and weights, we need them all in
        // memory when calling 'InlineApplyWeightsAndCorrections' so that we
        // can calculate the final visibilities to return correctly
        shared_data.ReadVisibilities(ms_reader,
                                     rows.visibilities[rows.n_rows_read].data(),
                                     rows.weights[rows.n_rows_read].data(),
                                     rows.model[rows.n_rows_read].data());

        constexpr bool kCacheEntireBeam = true;
        if constexpr (ApplyCorrections) {
          *antennas = std::make_pair(metadata.antenna1, metadata.antenna2);
          ++antennas;
          size_t time_offset = chunk_data.time_offsets.back();
          for (const auto& gridder : gridders) {
            gridder->LoadCorrections<ApplyBeam, HasH5Parm>(
                bands, metadata.time, metadata.field_id, time_offset,
                kCacheEntireBeam);
          };
          chunk_data.time_offsets.emplace_back(time_offset);
          ++time_offsets_offset;
        }

        ++rows.n_rows_read;
        ++n_chunk_rows_read;
        if (n_chunk_rows_read % 100000 == 0) {
          Logger::Debug << "n_chunk_rows_read: " +
                               std::to_string(n_chunk_rows_read) + "\n";
        }
        ms_reader.NextInputRow();
      }
      uvws_offset += uvws_stride * rows.n_rows_read;
      visibilities_offset += visibilities_stride * rows.n_rows_read;
      task_lane.write(std::move(rows));
    }
    task_lane.write_end();
  });

  // Per gridder mutex to prevent data race on currection sums inside
  // ApplyCorrections()
  std::vector<std::mutex> gridder_mutexes(gridders.size());
  // NB! This delibritely leads to overallocation of threads
  // As this still outperforms the alternative of not overlapping the IO.
  // Future changes should implement task stealing which would
  // fix this overallocation.
  std::vector<std::thread> thread_pool_process;
  thread_pool_process.reserve(available_cores_);
  for (size_t i = 0; i < available_cores_; ++i) {
    thread_pool_process.emplace_back([&] {
      BatchRowData rows;
      aocommon::UVector<float> image_weights(band.ChannelCount());
      // Per thread local counters to prevent contention and race conditions.
      // Propagate to global counters once per thread at end of processing loop.
      size_t local_gridded_visibility_count = 0;
      double local_total_weight = 0.0;
      double local_max_gridded_weight = 0.0;
      double local_visibility_weight_sum = 0.0;
      while (task_lane.read(rows)) {
        double* uvws = chunk_data.uvw.data() + rows.uvws_offset;
        std::complex<float>* visibilities =
            chunk_data.visibilities.data() + rows.visibilities_offset;
        size_t* time_offsets =
            chunk_data.time_offsets.data() + rows.time_offsets_offset;
        for (size_t n_buffer_index = 0; n_buffer_index < rows.n_rows_read;
             ++n_buffer_index) {
          MSProvider::MetaData& metadata = rows.metadata_[n_buffer_index];
          std::complex<float>* row_visibilities =
              rows.visibilities[n_buffer_index].data();
          float* row_weights = rows.weights[n_buffer_index].data();
          std::complex<float>* row_model = rows.model[n_buffer_index].data();
          uvws[0] = metadata.u_in_m;
          uvws[1] = metadata.v_in_m;
          uvws[2] = metadata.w_in_m;

          shared_data.ModifyVisibilities<true>(row_visibilities, uvws, band,
                                               row_model);
          shared_data.CalculateWeights(row_weights, image_weights.data(), uvws,
                                       band, selected_buffer);

          // Sum the corrections and apply the weights.
          // We store the appropriate time_offset to be used later along with
          // other required info when we apply the corrections
          if constexpr (ApplyCorrections) {
            size_t& time_offset = *time_offsets;
            for (size_t gridder_index = 0; gridder_index < gridders.size();
                 ++gridder_index) {
              std::lock_guard<std::mutex> lock(gridder_mutexes[gridder_index]);
              gridders[gridder_index]
                  ->ApplyCorrections<Mode, NParms, ModifierBehaviour::kSum,
                                     ApplyBeam, ApplyForward, HasH5Parm>(
                      ms_data.antenna_names.size(), row_visibilities, band,
                      metadata.data_desc_id, row_weights, metadata.antenna1,
                      metadata.antenna2, time_offset, image_weights.data());
            };
            ++time_offsets;
          }

          shared_data.ApplyWeights<Mode>(
              row_visibilities, band.ChannelCount(), row_weights,
              image_weights.data(), local_gridded_visibility_count,
              local_total_weight, local_max_gridded_weight,
              local_visibility_weight_sum);

          // When not applying corrections we collapse the polarizations.
          // When applying correction we need to keep them.
          if constexpr (!ApplyCorrections) {
            if (ms_data.ms_provider->NPolarizations() == 2) {
              internal::CollapseData<2>(band.ChannelCount(), row_visibilities,
                                        shared_data.Polarization());
            } else if (ms_data.ms_provider->NPolarizations() == 4) {
              internal::CollapseData<4>(band.ChannelCount(), row_visibilities,
                                        shared_data.Polarization());
            }
          }
          std::copy_n(row_visibilities, visibilities_stride, visibilities);
          visibilities += visibilities_stride;
          uvws += uvws_stride;
        }
      }
      shared_data.AddVisibilityCounts(
          local_gridded_visibility_count, local_total_weight,
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
      for (MsProviderCollection::MsData& ms_data :
           ms_provider_collection_.ms_data_vector_) {
        gridder->StartMeasurementSet(ms_provider_collection_.Count(), ms_data,
                                     false);
        ms_data.total_rows_processed += gridder->GridMeasurementSet(ms_data);
      }
      gridder->FinishInversionPass(pass_index);
    }
    gridder->FinishInversion();
  }
}

size_t MSGridderManager::GridChunk(
    bool apply_corrections, size_t n_vis_polarizations,
    const aocommon::MultiBandData& bands, const InversionChunkData& chunk_data,
    const aocommon::UVector<double>& frequencies,
    const MsProviderCollection::MsData& ms_data, size_t chunk_index,
    std::function<void()> signal_first_gridder_has_started,
    std::function<void()> signal_last_gridder_has_started) {
  Logger::Info << "Gridding " + std::to_string(chunk_data.n_rows) +
                      " rows for " + std::to_string(facet_tasks_.size()) +
                      " facets using " + std::to_string(available_cores_) +
                      " threads " +
                      std::to_string(available_cores_per_gridder_) +
                      " threads per gridder...\n";

  // TODO
  if (bands.BandCount() != 1)
    throw std::runtime_error("Can't use bulk inversion with BDA yet");
  const aocommon::BandData& band = *bands.begin();
  const size_t data_desc_id = *bands.DataDescIds().begin();

  ExecuteForAllGriddersWithNCores(
      available_cores_per_gridder_,
      [&](MsGridder* gridder, size_t facet_index) {
        GridChunkForFacet(*gridder, facet_index, apply_corrections,
                          n_vis_polarizations, band, data_desc_id, chunk_data,
                          frequencies, ms_data, chunk_index);
      },
      std::move(signal_first_gridder_has_started),
      std::move(signal_last_gridder_has_started));
  Logger::Info << "Finished gridding " + std::to_string(chunk_data.n_rows) +
                      " rows for " + std::to_string(facet_tasks_.size()) +
                      " facets.\n";
  return chunk_data.n_rows * facet_tasks_.size();
}

void MSGridderManager::GridChunkForFacet(
    MsGridder& gridder, size_t facet_index, bool apply_corrections,
    size_t n_vis_polarizations, const aocommon::BandData& band,
    size_t data_desc_id, const InversionChunkData& chunk_data,
    const aocommon::UVector<double>& frequencies,
    const MsProviderCollection::MsData& ms_data, size_t chunk_index) {
  Logger::Info << "Gridding facet " + std::to_string(facet_index) + "\n";

  const std::vector<std::complex<float>>& parm_response =
      gridder.GetVisibilityModifier().GetCachedParmResponse(
          ms_data.original_ms_index)[data_desc_id];

  std::shared_ptr<BeamResponseCacheChunk> beam_response =
      gridder.GetVisibilityModifier().TakeCachedBeamResponse(chunk_index);

  gridder.gridded_visibility_count_ = chunk_data.gridded_visibility_count;
  gridder.visibility_weight_sum_ = chunk_data.visibility_weight_sum;
  gridder.max_gridded_weight_ = chunk_data.max_gridded_weight;
  gridder.total_weight_ = chunk_data.total_weight;

  gridder.GridSharedMeasurementSetChunk(
      apply_corrections, n_vis_polarizations, chunk_data.n_rows,
      chunk_data.uvw.data(), frequencies.data(), band, data_desc_id,
      chunk_data.antennas.data(), chunk_data.visibilities.data(),
      apply_corrections ? chunk_data.time_offsets.data() + 1 : nullptr,
      ms_data.antenna_names.size(), parm_response, *beam_response.get());

  Logger::Info << "Done gridding facet " + std::to_string(facet_index) + "\n";
}

void MSGridderManager::ReadChunksForInvert(
    aocommon::Lane<InversionChunkData>& task_lane, size_t n_max_rows_in_memory,
    bool apply_corrections, MsProviderCollection::MsData& ms_data,
    MsGridderData& shared_data, const std::vector<MsGridder*>& gridders,
    size_t n_vis_polarizations, const bool* selected_buffer,
    size_t& n_total_chunks) {
  // TODO
  const aocommon::MultiBandData& bands = ms_data.ms_provider->SelectedBands();
  assert(bands.BandCount() == 1);
  const aocommon::BandData& band = *bands.begin();

  // We read chunks based on the maximum amount of rows we think we can fit
  // in memory at a time.
  Logger::Info << "Max " << n_max_rows_in_memory << " rows fit in memory.\n";
  const size_t n_total_rows_in_ms = ms_data.ms_provider->NRows();
  n_max_rows_in_memory = std::min(n_max_rows_in_memory, n_total_rows_in_ms);
  // We want two chunks of memory, one reading while one processes.
  // However estimate that we can read 50% of the next chunk before
  // the first is done gridding, so divide by 1.5 instead of 2.
  // NB! This should be revised in future if/when loading is faster
  // than gridding which would be ideal but is not the case currently.
  const size_t n_rows_per_chunk = n_max_rows_in_memory / 1.5;
  n_total_chunks = n_total_rows_in_ms / n_rows_per_chunk;
  size_t target_chunk_size = n_rows_per_chunk;
  // Compute the partial chunk remainder that might be left over.
  size_t n_rows_in_smaller_chunk = n_max_rows_in_memory % n_rows_per_chunk;
  if (n_rows_in_smaller_chunk > 0) {
    n_total_chunks += 1;
    target_chunk_size = n_rows_in_smaller_chunk;
  } else {
    n_rows_in_smaller_chunk = n_rows_per_chunk;
  }
  size_t chunk_index = 0;
  Logger::Info << "Reading " + std::to_string(n_total_chunks) +
                      " chunks with " +
                      std::to_string(n_rows_in_smaller_chunk) +
                      " rows in the first chunk and " +
                      std::to_string(n_rows_per_chunk) +
                      " rows per remaining chunk.\n";
  std::unique_ptr<MSReader> ms_reader = ms_data.ms_provider->MakeReader();
  while (ms_reader->CurrentRowAvailable()) {
    Logger::Info << "Loading " << target_chunk_size
                 << " rows into memory chunk " << chunk_index << ".\n";
    InversionChunkData chunk_data(target_chunk_size, band.ChannelCount(),
                                  n_vis_polarizations, apply_corrections);
    if (apply_corrections) {
      chunk_data.time_offsets.reserve(target_chunk_size + 1);
      chunk_data.time_offsets.push_back(0);
    }

    const size_t n_rows =
        ReadChunkForInvert(shared_data.GetGainMode(), apply_corrections,
                           gridders, ms_data, target_chunk_size, *ms_reader,
                           selected_buffer, chunk_data, shared_data);

    chunk_data.gridded_visibility_count = shared_data.gridded_visibility_count_;
    chunk_data.visibility_weight_sum = shared_data.visibility_weight_sum_;
    chunk_data.max_gridded_weight = shared_data.max_gridded_weight_;
    chunk_data.total_weight = shared_data.total_weight_;
    chunk_data.n_rows = n_rows;
    Logger::Debug << "Done loading chunk " << chunk_index << ".\n";

    task_lane.write(std::move(chunk_data));
    ++chunk_index;
    target_chunk_size = n_rows_per_chunk;
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
    const aocommon::UVector<double>& frequencies,
    const aocommon::MultiBandData& bands, MsProviderCollection::MsData& ms_data,
    size_t n_vis_polarizations, const size_t& n_chunks,
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
        Logger::Debug << "Gridding chunk " + std::to_string(chunk_index + 1) +
                             "/" + std::to_string(n_chunks) + ".\n";
        const bool is_last_chunk = chunk_index == n_chunks - 1;
        ms_data.total_rows_processed += GridChunk(
            apply_corrections, n_vis_polarizations, bands, chunk_data,
            frequencies, ms_data, chunk_index,
            [&]() { first_gridder_started.release(); },
            is_last_chunk ? signal_last_gridder_of_last_chunk_has_started
                          : nullptr);
        Logger::Debug << "Done gridding chunk " +
                             std::to_string(chunk_index + 1) + "/" +
                             std::to_string(n_chunks) + ".\n";
      },
      "Gridding");
}

void MSGridderManager::BatchInvert(
    std::function<void()> signal_last_gridding_work_has_started) {
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
    ExecuteForAllGridders(
        [=](MsGridder* gridder) { gridder->StartInversionPass(pass_index); });
    for (MsProviderCollection::MsData& ms_data : providers.ms_data_vector_) {
      MsGridderData shared_data(settings_);
      shared_data.CopyTaskData((*gridders[0]), solution_data_, ms_data);

      ExecuteForAllGridders(
          [&](MsGridder* gridder) {
            gridder->StartMeasurementSet(providers.Count(), ms_data, false);
          },
          false);
      shared_data.StartMeasurementSet(providers.Count(), ms_data, false);
      scheduler_->GetWorkerTaskQueue().WaitForIdle(available_cores_);

      const aocommon::MultiBandData& bands(
          ms_data.ms_provider->SelectedBands());
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
      if (bands.BandCount() != 1)
        throw std::runtime_error("Can't use bulk inversion with BDA yet");
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
        ReadChunksForInvert(task_lane, n_max_rows_in_memory, apply_corrections,
                            ms_data, shared_data, gridders, n_vis_polarizations,
                            selected_buffer.data(), n_total_chunks);
      });
      GridChunks(task_lane, apply_corrections, frequencies, bands, ms_data,
                 n_vis_polarizations, n_total_chunks,
                 signal_last_gridding_work_has_started);
      read_chunks_thread.join();
    }
    ExecuteForAllGridders(
        [=](MsGridder* gridder) { gridder->FinishInversionPass(pass_index); });
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
      for (MsProviderCollection::MsData& ms_data :
           ms_provider_collection_.ms_data_vector_) {
        gridder->StartMeasurementSet(ms_provider_collection_.Count(), ms_data,
                                     true);
        ms_data.total_rows_processed += gridder->PredictMeasurementSet(ms_data);
      }
      gridder->FinishPredictPass(pass_index);
    }
    gridder->FinishPredict();
  }
}

namespace internal {
/*
 * The input for this function is the predicted facet visibilities for a single
 * facet.
 * 1. Expand the visibilities.
 * 2. Apply corrections to the visibilities.
 * 3. Cumulatively add the expanded and corrected visibilities with those of
 * other facets.
 */
void ExpandAndCombineFacetVisibilities(
    size_t n_antennas, size_t n_channels, size_t n_rows,
    size_t n_vis_polarizations, const aocommon::MultiBandData& bands,
    const size_t* data_desc_ids, const size_t* antennas1,
    const size_t* antennas2, const size_t* field_ids, const double* times,
    const double* uvw, const std::complex<float>* facet_visibilities,
    MsGridder& gridder, std::complex<float>* combined_visibilities,
    std::vector<std::mutex>& sum_visibilities_mutexes,
    size_t mutex_chunk_size) {
  aocommon::UVector<std::complex<float>> visibilities_scratch(
      n_channels * n_vis_polarizations);
  const size_t n_chunks = n_rows / mutex_chunk_size + 1;
  for (size_t chunk = 0; chunk < n_chunks; ++chunk) {
    std::lock_guard<std::mutex> sum_lock(sum_visibilities_mutexes[chunk]);
    const size_t start = chunk * mutex_chunk_size;
    const size_t end = std::min(start + mutex_chunk_size, n_rows);
    for (size_t i = start; i < end; ++i) {
      switch (n_vis_polarizations) {
        case 1:
          std::copy_n(facet_visibilities, n_channels,
                      visibilities_scratch.data());
          break;
        case 2:
          internal::ExpandData<2>(n_channels, facet_visibilities,
                                  visibilities_scratch.data(),
                                  gridder.Polarization());
          break;
        case 4:
          internal::ExpandData<4>(n_channels, facet_visibilities,
                                  visibilities_scratch.data(),
                                  gridder.Polarization());
          break;
      }
      gridder.CorrectInstrumentalVisibilities(
          n_antennas, bands, *data_desc_ids, visibilities_scratch.data(), uvw,
          *field_ids, *antennas1, *antennas2, *times);

      // In case the value was not sampled in this pass, it has been set to
      // infinite and should not overwrite the current value in the set.
      // NB! This is only true for multi pass gridders (wstacking gridder) other
      // gridders could consider skipping this check if it becomes important for
      // performance.
      for (size_t i = 0; i < n_channels * n_vis_polarizations; ++i) {
        if (std::isfinite(visibilities_scratch[i].real())) {
          combined_visibilities[i] += visibilities_scratch[i];
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
}  // namespace internal

size_t MSGridderManager::PredictChunk(
    const PredictionChunkData& chunk_data, size_t n_vis_polarizations,
    size_t n_antennas, std::vector<std::complex<float>>& combined_visibilities,
    const aocommon::UVector<double>& frequencies,
    MsProviderCollection::MsData& ms_data,
    std::function<void()> signal_first_predict_has_started,
    std::function<void()> signal_last_predict_has_started) {
  Logger::Info << "Predicting " + std::to_string(chunk_data.n_rows) +
                      " rows for " + std::to_string(facet_tasks_.size()) +
                      " facets using " + std::to_string(available_cores_) +
                      " threads " +
                      std::to_string(available_cores_per_gridder_) +
                      " threads per gridder...\n";

  // As all facets are summing their values into combined_visibilities we have
  // to protect it with a mutex to avoid data corruption.
  // Chunk into ranges and use a mutex per range instead of a single mutex to
  // strike a balance between gridder lock contention and time spent waiting.
  constexpr size_t kMutexChunkSize = 2048;
  std::vector<std::mutex> sum_visibilities_mutexes(
      (chunk_data.n_rows / kMutexChunkSize) + 1);

  ExecuteForAllGriddersWithNCores(
      available_cores_per_gridder_,
      [&](MsGridder* gridder, size_t facet_index) {
        PredictChunkForFacet(*gridder, facet_index, chunk_data,
                             n_vis_polarizations, n_antennas,
                             combined_visibilities, frequencies, ms_data,
                             sum_visibilities_mutexes, kMutexChunkSize);
      },
      std::move(signal_first_predict_has_started),
      std::move(signal_last_predict_has_started));

  Logger::Info << "Finished Predicting " + std::to_string(chunk_data.n_rows) +
                      " rows for " + std::to_string(facet_tasks_.size()) +
                      " facets.\n";

  return chunk_data.n_rows * facet_tasks_.size();
}

void MSGridderManager::PredictChunkForFacet(
    MsGridder& gridder, size_t facet_index,
    const PredictionChunkData& chunk_data, size_t n_vis_polarizations,
    size_t n_antennas, std::vector<std::complex<float>>& combined_visibilities,
    const aocommon::UVector<double>& frequencies,
    MsProviderCollection::MsData& ms_data,
    std::vector<std::mutex>& sum_visibilities_mutexes,
    size_t mutex_chunk_size) {
  using internal::ExpandAndCombineFacetVisibilities;

  Logger::Info << "Predicting facet " + std::to_string(facet_index) + "\n";

  // TODO
  if (!ms_data.ms_provider->IsRegular())
    throw std::runtime_error("Can't use bulk inversion with BDA yet");
  const size_t n_channels = frequencies.size();

  aocommon::UVector<std::complex<float>> facet_visibilities(chunk_data.n_rows *
                                                            n_channels);
  gridder.PredictChunk(chunk_data.n_rows, n_channels, frequencies.data(),
                       chunk_data.uvws.data(), facet_visibilities.data());

  ExpandAndCombineFacetVisibilities(
      n_antennas, n_channels, chunk_data.n_rows, n_vis_polarizations,
      ms_data.ms_provider->SelectedBands(), chunk_data.data_desc_ids.data(),
      chunk_data.antennas1.data(), chunk_data.antennas2.data(),
      chunk_data.field_ids.data(), chunk_data.times.data(),
      chunk_data.uvws.data(), facet_visibilities.data(), gridder,
      combined_visibilities.data(), sum_visibilities_mutexes, mutex_chunk_size);

  Logger::Info << "Done predicting facet " + std::to_string(facet_index) + "\n";
}

void MSGridderManager::PredictChunks(
    aocommon::Lane<PredictionChunkData>& task_lane,
    const aocommon::UVector<double>& frequencies,
    MsProviderCollection::MsData& ms_data, size_t n_vis_polarizations,
    bool add_assign_model, const size_t& n_chunks,
    std::function<void()> signal_last_predict_of_last_chunk_has_started) {
  // TODO
  if (!ms_data.ms_provider->IsRegular())
    throw std::runtime_error("Can't use bulk inversion with BDA yet");
  const aocommon::BandData& band =
      *ms_data.ms_provider->SelectedBands().begin();

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
        Logger::Debug << "Predicting chunk " + std::to_string(chunk_index + 1) +
                             "/" + std::to_string(n_chunks) + ".\n";
        bool is_last_chunk = chunk_index == n_chunks - 1;
        ms_data.total_rows_processed += PredictChunk(
            chunk_data, n_vis_polarizations, ms_data.antenna_names.size(),
            combined_visibilities, frequencies, ms_data,
            [&]() { first_predict_started.release(); },
            is_last_chunk ? signal_last_predict_of_last_chunk_has_started
                          : nullptr);
        Logger::Debug << "Done predicting chunk " +
                             std::to_string(chunk_index + 1) + "/" +
                             std::to_string(n_chunks) + ".\n";

        // Do a single write for the combined/expanded chunk of visibilities
        // of all facets.
        Logger::Debug << "Writing predicted chunk " +
                             std::to_string(chunk_index + 1) + "/" +
                             std::to_string(n_chunks) + ".\n";
        std::complex<float>* visibilities = combined_visibilities.data();
        const size_t stride = band.ChannelCount() * n_vis_polarizations;
        for (size_t row = 0; row != chunk_data.n_rows; ++row) {
          ms_data.ms_provider->WriteModel(visibilities, add_assign_model);
          ms_data.ms_provider->NextOutputRow();
          visibilities += stride;
        }
        Logger::Debug << "Done writing predicted chunk " +
                             std::to_string(chunk_index + 1) + "/" +
                             std::to_string(n_chunks) + ".\n";
      },
      "Predict");
}

void MSGridderManager::ReadChunksForPredict(
    aocommon::Lane<PredictionChunkData>& task_lane, size_t n_max_rows_in_memory,
    MsProviderCollection::MsData& ms_data, MsGridderData& shared_data,
    const std::vector<MsGridder*>& gridders, const aocommon::BandData band,
    size_t n_vis_polarizations, const bool* selected_buffer,
    size_t& n_total_chunks) {
  // We read chunks based on the maximum amount of rows we think we can fit
  // in memory at a time.
  Logger::Info << "Max " << n_max_rows_in_memory << " rows fit in memory.\n";
  const size_t n_total_rows_in_ms = ms_data.ms_provider->NRows();
  n_max_rows_in_memory = std::min(n_max_rows_in_memory, n_total_rows_in_ms);
  // We want two chunks of memory, one reading while one processes.
  // However estimate that we can read 50% of the next chunk before
  // the first is done gridding, so divide by 1.5 instead of 2.
  // NB! This should be revised in future if/when loading is faster
  // than gridding which would be ideal but is not the case currently.
  const size_t n_rows_per_chunk = n_max_rows_in_memory / 1.5;
  n_total_chunks = (n_total_rows_in_ms / n_rows_per_chunk);
  size_t target_chunk_size = n_rows_per_chunk;
  // Compute the partial chunk remainder that might be left over.
  size_t n_rows_in_smaller_chunk = n_max_rows_in_memory % n_rows_per_chunk;
  if (n_rows_in_smaller_chunk > 0) {
    n_total_chunks += 1;
    target_chunk_size = n_rows_in_smaller_chunk;
  } else {
    n_rows_in_smaller_chunk = n_rows_per_chunk;
  }
  size_t chunk_index = 0;
  Logger::Info << "Reading " + std::to_string(n_total_chunks) +
                      " chunks with " +
                      std::to_string(n_rows_in_smaller_chunk) +
                      " rows in first chunk and " +
                      std::to_string(n_rows_per_chunk) +
                      " rows per remaining chunk.\n";
  // Set provider up for writing
  ms_data.ms_provider->ReopenRW();
  ms_data.ms_provider->ResetWritePosition();
  std::unique_ptr<MSReader> ms_reader = ms_data.ms_provider->MakeReader();
  while (ms_reader->CurrentRowAvailable()) {
    Logger::Info << "Loading " << target_chunk_size
                 << " rows into memory chunk " << chunk_index << ".\n";
    PredictionChunkData chunk_data(target_chunk_size);
    while (ms_reader->CurrentRowAvailable() &&
           chunk_data.n_rows < target_chunk_size) {
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
    Logger::Debug << "Done loading chunk " << chunk_index << ".\n";

    task_lane.write(std::move(chunk_data));
    ++chunk_index;
    target_chunk_size = n_rows_per_chunk;
  }
  Logger::Info << "All predict rows loaded.\n";
  task_lane.write_end();
}

void MSGridderManager::BatchPredict(
    std::function<void()> signal_last_predict_work_has_started) {
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
    ExecuteForAllGridders(
        [=](MsGridder* gridder) { gridder->StartPredictPass(pass_index); });
    for (MsProviderCollection::MsData& ms_data : providers.ms_data_vector_) {
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
      const aocommon::BandData& band =
          *ms_data.ms_provider->SelectedBands().begin();

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
      size_t shared_mem_per_row =
          sizeof(std::complex<float>) * n_channels * n_vis_polarizations;
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
      const size_t n_max_overall_rows_in_memory =
          gridders[0]->CalculateMaxRowsInMemory(
              available_memory_, constant_mem, additional_per_vis_row_mem,
              per_row_uvw_memory_consumption, n_channels, 1);
      // Finally approximate that again to multiple gridders to get maximum rows
      // per gridder.
      const size_t n_max_rows_in_memory =
          n_max_overall_rows_in_memory / gridders.size();

      aocommon::UVector<double> frequencies(n_channels);
      for (size_t i = 0; i != frequencies.size(); ++i) {
        frequencies[i] = band.ChannelFrequency(i);
      }
      aocommon::UVector<bool> selected_buffer(n_channels, true);

      // Iterate over data in chunks until all visibilities have been predicted.
      aocommon::Lane<PredictionChunkData> task_lane(1);
      size_t n_total_chunks = 0;
      std::thread read_chunks_thread([&] {
        ReadChunksForPredict(task_lane, n_max_rows_in_memory, ms_data,
                             shared_data, gridders, band, n_vis_polarizations,
                             selected_buffer.data(), n_total_chunks);
      });
      PredictChunks(task_lane, frequencies, ms_data, n_vis_polarizations,
                    gridders[0]->ShouldAddAssignModel(), n_total_chunks,
                    signal_last_predict_work_has_started);
      read_chunks_thread.join();
    }
    ExecuteForAllGridders(
        [=](MsGridder* gridder) { gridder->FinishPredictPass(pass_index); });
  }
  ExecuteForAllGridders([](MsGridder* gridder) { gridder->FinishPredict(); });
}

void MSGridderManager::ProcessResults(std::mutex& result_mutex,
                                      GriddingResult& result,
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
    facet_result->effectiveGriddedVisibilityCount =
        gridder->EffectiveGriddedVisibilityCount();
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

std::unique_ptr<MsGridder> MSGridderManager::ConstructGridder(
    const Resources& resources) {
  switch (settings_.gridderType) {
    case GridderType::IDG:
      return std::make_unique<IdgMsGridder>(settings_, resources,
                                            ms_provider_collection_);
    case GridderType::FacetIDG:
      return std::make_unique<FacetIdgMsGridder>(settings_, resources,
                                                 ms_provider_collection_);
    case GridderType::WGridder:
      return std::make_unique<WGriddingMSGridder>(
          settings_, resources, ms_provider_collection_, false);
    case GridderType::TunedWGridder:
      return std::make_unique<WGriddingMSGridder>(
          settings_, resources, ms_provider_collection_, true);
    case GridderType::WTowers:
#ifdef BUILD_WTOWERS
      return std::make_unique<WTowersMsGridder>(settings_, resources,
                                                ms_provider_collection_);
#else
      throw std::runtime_error("w-towers gridder is not available");
#endif
    case GridderType::DirectFT:
      switch (settings_.directFTPrecision) {
        case DirectFTPrecision::Float:
          return std::make_unique<DirectMSGridder<float>>(
              settings_, resources, ms_provider_collection_);
        case DirectFTPrecision::Double:
          return std::make_unique<DirectMSGridder<double>>(
              settings_, resources, ms_provider_collection_);
        case DirectFTPrecision::LongDouble:
          return std::make_unique<DirectMSGridder<long double>>(
              settings_, resources, ms_provider_collection_);
      }
      break;
    case GridderType::WStacking:
      return std::make_unique<WSMSGridder>(settings_, resources,
                                           ms_provider_collection_);
  }
  return {};
}

void MSGridderManager::InitializeGridderForTask(
    MsGridder& gridder, const GriddingTask& task,
    GriddingTaskManager* writer_lock_manager) {
  gridder.SetGridMode(settings_.gridMode);

  gridder.SetFacetGroupIndex(task.facetGroupIndex);
  gridder.SetImagePadding(settings_.imagePadding);
  gridder.SetPhaseCentreDec(task.observationInfo.phaseCentreDec);
  gridder.SetPhaseCentreRA(task.observationInfo.phaseCentreRA);

  if (settings_.hasShift) {
    double main_image_dl = 0.0;
    double main_image_dm = 0.0;
    aocommon::ImageCoordinates::RaDecToLM(settings_.shiftRA, settings_.shiftDec,
                                          task.observationInfo.phaseCentreRA,
                                          task.observationInfo.phaseCentreDec,
                                          main_image_dl, main_image_dm);
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

void MSGridderManager::InitializeGridderForFacet(
    bool is_first_polarization, MsGridder& gridder,
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
    gridder.GetVisibilityModifier().SetFacetDirection(facet->RA(),
                                                      facet->Dec());
  } else {
    gridder.SetImageWidth(settings_.paddedImageWidth);
    gridder.SetImageHeight(settings_.paddedImageHeight);
    gridder.SetTrimSize(settings_.trimmedImageWidth,
                        settings_.trimmedImageHeight);
  }
  gridder.SetLShift(facet_task.l_shift);
  gridder.SetMShift(facet_task.m_shift);

  std::unique_ptr<MetaDataCache> cache = std::move(facet_task.cache);
  if (!cache) cache = std::make_unique<MetaDataCache>();
  gridder.SetMetaDataCache(std::move(cache));
}

}  // namespace wsclean
