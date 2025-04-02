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
#include "../main/settings.h"
#include "../structures/resources.h"
#include "../wgridder/wgriddingmsgridder.h"
#include "../wtowers/wtowersmsgridder.h"

using aocommon::Logger;

namespace wsclean {

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
  for (size_t facet_index : facet_indices) {
    assert(facet_index < task.facets.size());

    // Create a new gridder for each facet / sub-task, since gridders do not
    // support reusing them for multiple tasks.
    std::unique_ptr<MsGridder> gridder =
        ConstructGridder(resources.GetPart(task.num_parallel_gridders_));
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

    InitializeGridderForFacet(*gridder, *facet_task);

    facet_tasks_.emplace_back(
        GriddingFacetTask{std::move(gridder), facet_task, facet_result});
  }
}

size_t MSGridderManager::ReadChunkForInvert(
    GainMode gain_mode, bool apply_corrections,
    const std::vector<MsGridder*>& gridders,
    MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
    MSReader& ms_reader, const aocommon::BandData band,
    const bool* selected_buffer, InversionChunkData& chunk_data,
    MsGridderData& shared_data) {
  const size_t n_parms = gridders[0]->NumValuesPerSolution();
  switch (gain_mode) {
    case GainMode::kXX:
      return ReadChunkForInvertImplementation<GainMode::kXX>(
          n_parms, apply_corrections, gridders, ms_data, n_chunk_rows,
          ms_reader, band, selected_buffer, chunk_data, shared_data);
      break;
    case GainMode::kYY:
      return ReadChunkForInvertImplementation<GainMode::kYY>(
          n_parms, apply_corrections, gridders, ms_data, n_chunk_rows,
          ms_reader, band, selected_buffer, chunk_data, shared_data);
      break;
    case GainMode::k2VisDiagonal:
      return ReadChunkForInvertImplementation<GainMode::k2VisDiagonal>(
          n_parms, apply_corrections, gridders, ms_data, n_chunk_rows,
          ms_reader, band, selected_buffer, chunk_data, shared_data);
      break;
    case GainMode::kTrace:
      return ReadChunkForInvertImplementation<GainMode::kTrace>(
          n_parms, apply_corrections, gridders, ms_data, n_chunk_rows,
          ms_reader, band, selected_buffer, chunk_data, shared_data);
      break;
    case GainMode::kFull:
      return ReadChunkForInvertImplementation<GainMode::kFull>(
          n_parms, apply_corrections, gridders, ms_data, n_chunk_rows,
          ms_reader, band, selected_buffer, chunk_data, shared_data);
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
    MSReader& ms_reader, const aocommon::BandData band,
    const bool* selected_buffer, InversionChunkData& chunk_data,
    MsGridderData& shared_data) {
  if (n_parms == 2) {
    return ReadChunkForInvertImplementation<Mode, 2>(
        apply_corrections, gridders, ms_data, n_chunk_rows, ms_reader, band,
        selected_buffer, chunk_data, shared_data);
  } else {
    return ReadChunkForInvertImplementation<Mode, 4>(
        apply_corrections, gridders, ms_data, n_chunk_rows, ms_reader, band,
        selected_buffer, chunk_data, shared_data);
  }
}

template <GainMode Mode, size_t NParms>
size_t MSGridderManager::ReadChunkForInvertImplementation(
    bool apply_corrections, const std::vector<MsGridder*>& gridders,
    MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
    MSReader& ms_reader, const aocommon::BandData band,
    const bool* selected_buffer, InversionChunkData& chunk_data,
    MsGridderData& shared_data) {
  if (apply_corrections) {
    return ReadChunkForInvertImplementation<Mode, NParms, true>(
        gridders, ms_data, n_chunk_rows, ms_reader, band, selected_buffer,
        chunk_data, shared_data);
  } else {
    return ReadChunkForInvertImplementation<Mode, NParms, false>(
        gridders, ms_data, n_chunk_rows, ms_reader, band, selected_buffer,
        chunk_data, shared_data);
  }
}

template <GainMode Mode, size_t NParms, bool ApplyCorrections>
size_t MSGridderManager::ReadChunkForInvertImplementation(
    const std::vector<MsGridder*>& gridders,
    MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
    MSReader& ms_reader, const aocommon::BandData band,
    const bool* selected_buffer, InversionChunkData& chunk_data,
    MsGridderData& shared_data) {
  if constexpr (ApplyCorrections) {
    const bool apply_beam = settings_.applyFacetBeam || settings_.gridWithBeam;
    if (apply_beam) {
      return ReadChunkForInvertImplementation<Mode, NParms, ApplyCorrections,
                                              true>(
          gridders, ms_data, n_chunk_rows, ms_reader, band, selected_buffer,
          chunk_data, shared_data);
    } else {
      return ReadChunkForInvertImplementation<Mode, NParms, ApplyCorrections,
                                              false>(
          gridders, ms_data, n_chunk_rows, ms_reader, band, selected_buffer,
          chunk_data, shared_data);
    }
  } else {
    return ReadChunkForInvertImplementation<Mode, NParms, ApplyCorrections,
                                            false, false, false>(
        gridders, ms_data, n_chunk_rows, ms_reader, band, selected_buffer,
        chunk_data, shared_data);
  }
}

template <GainMode Mode, size_t NParms, bool ApplyCorrections, bool ApplyBeam>
size_t MSGridderManager::ReadChunkForInvertImplementation(
    const std::vector<MsGridder*>& gridders,
    MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
    MSReader& ms_reader, const aocommon::BandData band,
    const bool* selected_buffer, InversionChunkData& chunk_data,
    MsGridderData& shared_data) {
  const bool apply_forward =
      gridders[0]->GetPsfMode() == PsfMode::kDirectionDependent;
  if (apply_forward) {
    return ReadChunkForInvertImplementation<Mode, NParms, ApplyCorrections,
                                            ApplyBeam, true>(
        gridders, ms_data, n_chunk_rows, ms_reader, band, selected_buffer,
        chunk_data, shared_data);
  } else {
    return ReadChunkForInvertImplementation<Mode, NParms, ApplyCorrections,
                                            ApplyBeam, false>(
        gridders, ms_data, n_chunk_rows, ms_reader, band, selected_buffer,
        chunk_data, shared_data);
  }
}

template <GainMode Mode, size_t NParms, bool ApplyCorrections, bool ApplyBeam,
          bool ApplyForward>
size_t MSGridderManager::ReadChunkForInvertImplementation(
    const std::vector<MsGridder*>& gridders,
    MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
    MSReader& ms_reader, const aocommon::BandData band,
    const bool* selected_buffer, InversionChunkData& chunk_data,
    MsGridderData& shared_data) {
  const bool has_h5_parm = gridders[0]->visibility_modifier_.HasH5Parm();
  if (has_h5_parm) {
    return ReadChunkForInvertImplementation<Mode, NParms, ApplyCorrections,
                                            ApplyBeam, ApplyForward, true>(
        gridders, ms_data, n_chunk_rows, ms_reader, band, selected_buffer,
        chunk_data, shared_data);
  } else {
    return ReadChunkForInvertImplementation<Mode, NParms, ApplyCorrections,
                                            ApplyBeam, ApplyForward, false>(
        gridders, ms_data, n_chunk_rows, ms_reader, band, selected_buffer,
        chunk_data, shared_data);
  }
}

template <GainMode Mode, size_t NParms, bool ApplyCorrections, bool ApplyBeam,
          bool ApplyForward, bool HasH5Parm>
size_t MSGridderManager::ReadChunkForInvertImplementation(
    const std::vector<MsGridder*>& gridders,
    MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
    MSReader& ms_reader, const aocommon::BandData band,
    const bool* selected_buffer, InversionChunkData& chunk_data,
    MsGridderData& shared_data) {
  size_t n_chunk_rows_read = 0;

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

  // Allow reading to get a but ahead of processing but not by too much.
  aocommon::Lane<BatchRowData> task_lane(available_cores_ * 2);

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

        if constexpr (ApplyCorrections) {
          *antennas = std::make_pair(metadata.antenna1, metadata.antenna2);
          ++antennas;
          size_t time_offset = chunk_data.time_offsets.back();
          for (const auto& gridder : gridders) {
            gridder->LoadCorrections<ApplyBeam, HasH5Parm>(
                band, metadata.time, metadata.field_id, time_offset);
          };
          chunk_data.time_offsets.emplace_back(time_offset);
          ++time_offsets_offset;
        }

        ++rows.n_rows_read;
        ++n_chunk_rows_read;
        if (n_chunk_rows_read % 100000 == 0) {
          Logger::Debug << "n_chunk_rows_read: " << n_chunk_rows_read << "\n";
        }
        ms_reader.NextInputRow();
      }
      uvws_offset += uvws_stride * rows.n_rows_read;
      visibilities_offset += visibilities_stride * rows.n_rows_read;
      task_lane.write(std::move(rows));
    }
    task_lane.write_end();
  });

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

          shared_data.CalculateWeights(uvws, row_visibilities, band,
                                       row_weights, row_model, selected_buffer,
                                       image_weights.data());

          // Sum the corrections and apply the weights.
          // We store the appropriate time_offset to be used later along with
          // other required info when we apply the corrections
          if constexpr (ApplyCorrections) {
            size_t& time_offset = *time_offsets;
            for (const auto& gridder : gridders) {
              gridder->ApplyCorrections<Mode, NParms, ModifierBehaviour::kSum,
                                        ApplyBeam, ApplyForward, HasH5Parm>(
                  ms_data.antenna_names.size(), row_visibilities, band,
                  row_weights, metadata.antenna1, metadata.antenna2,
                  time_offset, image_weights.data());
            };
            ++time_offsets;
          }

          shared_data.ApplyWeights<Mode>(row_visibilities, band.ChannelCount(),
                                         row_weights, image_weights.data());

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
    });
  }
  read_rows_thread.join();
  for (std::thread& thread : thread_pool_process) {
    thread.join();
  }
  return n_chunk_rows_read;
}

void MSGridderManager::Invert() {
  assert(facet_tasks_.size() == 1);
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
    aocommon::TaskQueue<std::function<void()>>& task_queue,
    size_t num_parallel_gridders, bool apply_corrections,
    size_t n_vis_polarizations, const aocommon::BandData& band,
    const InversionChunkData& chunk_data,
    const aocommon::UVector<double>& frequencies,
    const MsProviderCollection::MsData& ms_data) {
  Logger::Info << "Gridding " + std::to_string(chunk_data.n_rows) +
                      " rows for " + std::to_string(facet_tasks_.size()) +
                      " facets using " + std::to_string(num_parallel_gridders) +
                      " threads...\n";

  ExecuteForAllGriddersWithNCores(
      task_queue, num_parallel_gridders,
      [&](MsGridder* gridder, size_t facet_index) {
        Logger::Info << "Gridding facet " + std::to_string(facet_index) + "\n";

        const std::vector<std::complex<float>>& parm_response =
            gridder->GetVisibilityModifier().GetCachedParmResponse(
                ms_data.original_ms_index);

        gridder->gridded_visibility_count_ =
            chunk_data.gridded_visibility_count;
        gridder->visibility_weight_sum_ = chunk_data.visibility_weight_sum;
        gridder->max_gridded_weight_ = chunk_data.max_gridded_weight;
        gridder->total_weight_ = chunk_data.total_weight;

        gridder->GridSharedMeasurementSetChunk(
            apply_corrections, n_vis_polarizations, chunk_data.n_rows,
            chunk_data.uvw.data(), frequencies.data(), band,
            chunk_data.antennas.data(), chunk_data.visibilities.data(),
            apply_corrections ? chunk_data.time_offsets.data() + 1 : nullptr,
            ms_data.antenna_names.size(), parm_response);
        Logger::Info << "Done gridding facet " + std::to_string(facet_index) +
                            "\n";
      });
  Logger::Info << "Finished gridding " + std::to_string(chunk_data.n_rows) +
                      " rows for " + std::to_string(facet_tasks_.size()) +
                      " facets using " + std::to_string(num_parallel_gridders) +
                      " threads...\n";
  return chunk_data.n_rows * facet_tasks_.size();
}

void MSGridderManager::ReadChunksForInvert(
    aocommon::Lane<InversionChunkData>& task_lane, size_t n_max_rows_in_memory,
    bool apply_corrections, MsProviderCollection::MsData& ms_data,
    MsGridderData& shared_data, const std::vector<MsGridder*>& gridders,
    const aocommon::BandData band, size_t n_vis_polarizations,
    const bool* selected_buffer) {
  // We read chunks based on the maximum amount of rows we think we can fit
  // in memory at a time.
  Logger::Info << "Max " << n_max_rows_in_memory << " rows fit in memory.\n";
  SynchronizedMS ms = ms_data.ms_provider->MS();
  const size_t n_total_rows_in_ms = ms->nrow();
  n_max_rows_in_memory = std::min(n_max_rows_in_memory, n_total_rows_in_ms);
  // We want two chunks of memory, one reading while one processes.
  // However estimate that we can read 50% of the next chunk before
  // the first is done gridding, so divide by 1.5 instead of 2.
  // NB! This should be revised in future if/when loading is faster
  // than gridding which would be ideal but is not the case currently.
  const size_t n_rows_per_chunk = n_max_rows_in_memory / 1.5;
  size_t total_chunks = (n_total_rows_in_ms / n_rows_per_chunk);
  size_t target_chunk_size = n_rows_per_chunk;
  // Compute the partial chunk remainder that might be left over.
  size_t n_rows_in_smaller_chunk = n_max_rows_in_memory % n_rows_per_chunk;
  if (n_rows_in_smaller_chunk > 0) {
    total_chunks += 1;
    target_chunk_size = n_rows_in_smaller_chunk;
  } else {
    n_rows_in_smaller_chunk = n_rows_per_chunk;
  }
  size_t chunk_index = 1;
  Logger::Info << "Reading " << total_chunks << " chunks with "
               << n_rows_in_smaller_chunk << " rows in first chunk and "
               << n_rows_per_chunk << " rows per remaining chunk.\n";
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
                           band, selected_buffer, chunk_data, shared_data);

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
  Logger::Info << "All gridding rows loaded.\n";
  task_lane.write_end();
}

void MSGridderManager::GridChunks(
    aocommon::Lane<InversionChunkData>& task_lane,
    const size_t num_parallel_gridders, const bool apply_corrections,
    aocommon::TaskQueue<std::function<void()>>& task_queue,
    const aocommon::UVector<double>& frequencies,
    const aocommon::BandData& band, MsProviderCollection::MsData& ms_data,
    size_t n_vis_polarizations) {
  InversionChunkData chunk_data;
  size_t chunk_index = 1;
  while (task_lane.read(chunk_data)) {
    Logger::Info << "Gridding chunk" << chunk_index << ".\n";
    ms_data.total_rows_processed +=
        GridChunk(task_queue, num_parallel_gridders, apply_corrections,
                  n_vis_polarizations, band, chunk_data, frequencies, ms_data);
    Logger::Info << "Done gridding chunk" << chunk_index << ".\n";
    ++chunk_index;
  }
  Logger::Info << "All gridding rows processed.\n";
}

void MSGridderManager::BatchInvert(size_t num_parallel_gridders) {
  assert(facet_tasks_.size() > 1);
  InitializeMSDataVectors();

  MsProviderCollection& providers = ms_provider_collection_;

  aocommon::TaskQueue<std::function<void()>> task_queue;
  std::vector<std::thread> thread_pool;
  thread_pool.reserve(available_cores_);
  for (size_t i = 0; i < available_cores_; ++i) {
    thread_pool.emplace_back([&] {
      std::function<void()> operation;
      while (task_queue.Pop(operation)) {
        operation();
      }
    });
  }

  std::vector<MsGridder*> gridders;
  gridders.reserve(facet_tasks_.size());
  for (const GriddingFacetTask& task : facet_tasks_) {
    const std::unique_ptr<MsGridder>& gridder = task.facet_gridder;
    gridders.emplace_back(gridder.get());
  }

  ExecuteForAllGridders(task_queue, [=](MsGridder* gridder) {
    gridder->CalculateOverallMetaData();
    gridder->StartInversion();
  });
  const size_t n_inversion_passes = gridders[0]->GetNInversionPasses();
  for (size_t pass_index = 0; pass_index < n_inversion_passes; ++pass_index) {
    ExecuteForAllGridders(task_queue, [=](MsGridder* gridder) {
      gridder->StartInversionPass(pass_index);
    });
    for (MsProviderCollection::MsData& ms_data : providers.ms_data_vector_) {
      MsGridderData shared_data(settings_);
      shared_data.CopyTaskData((*gridders[0]), solution_data_, ms_data);

      ExecuteForAllGridders(
          task_queue,
          [&](MsGridder* gridder) {
            gridder->StartMeasurementSet(providers.Count(), ms_data, false);
          },
          false);
      shared_data.StartMeasurementSet(providers.Count(), ms_data, false);
      task_queue.WaitForIdle(available_cores_);

      const aocommon::BandData band(ms_data.SelectedBand());
      const size_t n_channels = band.ChannelCount();
      const size_t n_vis_polarizations = ms_data.ms_provider->NPolarizations();

      // We need to sum constant memory usage up across all gridders as each
      // gridder has its own internal memory usage based on image size, but
      // perVisMem will always be the same as its shared across all gridders
      size_t constant_mem = 0;
      for (MsGridder* gridder : gridders) {
        constant_mem += gridder->CalculateConstantMemory();
      }
      // We incur these additional per row memory overheads with data that we
      // have to cache for later in order to apply the corrections
      size_t additional_per_vis_mem = 0;
      bool apply_corrections = gridders[0]->WillApplyCorrections();
      if (apply_corrections) {
        // For each row we have to store an antenna pair and a solution time
        // offset
        additional_per_vis_mem = sizeof(size_t) * 3;
      }

      const size_t n_max_rows_in_memory = gridders[0]->CalculateMaxRowsInMemory(
          available_memory_, constant_mem, additional_per_vis_mem, n_channels,
          apply_corrections ? n_vis_polarizations : 1);

      aocommon::UVector<double> frequencies(n_channels);
      for (size_t i = 0; i != frequencies.size(); ++i) {
        frequencies[i] = band.ChannelFrequency(i);
      }
      aocommon::UVector<bool> selected_buffer(n_channels, true);

      // Iterate over data in chunks until all visibilities
      // have been gridded.
      aocommon::Lane<InversionChunkData> task_lane(1);
      std::thread grid_chunks_thread([&] {
        GridChunks(task_lane, num_parallel_gridders, apply_corrections,
                   task_queue, frequencies, band, ms_data, n_vis_polarizations);
      });
      ReadChunksForInvert(task_lane, n_max_rows_in_memory, apply_corrections,
                          ms_data, shared_data, gridders, band,
                          n_vis_polarizations, selected_buffer.data());
      grid_chunks_thread.join();
    }
    ExecuteForAllGridders(task_queue, [=](MsGridder* gridder) {
      gridder->FinishInversionPass(pass_index);
    });
  }
  ExecuteForAllGridders(task_queue,
                        [](MsGridder* gridder) { gridder->FinishInversion(); });

  // Clean up the thread pool
  task_queue.Finish();
  for (std::thread& thread : thread_pool) {
    thread.join();
  }
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
      gridder->FinishPredictPass();
    }
    gridder->FinishPredict();
  }
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
      result.griddedVisibilityCount += gridder->GriddedVisibilityCount();
      result.visibilityWeightSum += gridder->VisibilityWeightSum();
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

void MSGridderManager::SortFacetTasks() {
  // Image size is probably an imperfect approximation of job length but should
  // on average be better than not sorting at all.
  std::sort(
      facet_tasks_.begin(), facet_tasks_.end(),
      [](const GriddingFacetTask& a, const GriddingFacetTask& b) {
        return a.facet_gridder->ImageWidth() * a.facet_gridder->ImageHeight() >
               b.facet_gridder->ImageWidth() * b.facet_gridder->ImageHeight();
      });
}

std::unique_ptr<MsGridder> MSGridderManager::ConstructGridder(
    const Resources& resources) {
  switch (settings_.gridderType) {
    case GridderType::IDG:
      return std::make_unique<IdgMsGridder>(settings_, resources,
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
    MsGridder& gridder, GriddingTask::FacetData& facet_task) {
  const schaapcommon::facets::Facet* facet = facet_task.facet.get();
  gridder.SetIsFacet(facet != nullptr);
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
