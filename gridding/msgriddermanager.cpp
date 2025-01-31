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
    aocommon::TaskQueue<std::function<void()>>& task_queue,
    const std::vector<MsGridder*>& gridders,
    MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
    MSReader& ms_reader, const aocommon::BandData band,
    const bool* selected_buffer, RowData& row_data, ChunkData& chunk_data,
    MsGridderData& shared_data) {
  switch (gain_mode) {
    case GainMode::kXX:
      return ReadChunkForInvertImplementation<GainMode::kXX>(
          apply_corrections, task_queue, gridders, ms_data, n_chunk_rows,
          ms_reader, band, selected_buffer, row_data, chunk_data, shared_data);
      break;
    case GainMode::kYY:
      return ReadChunkForInvertImplementation<GainMode::kYY>(
          apply_corrections, task_queue, gridders, ms_data, n_chunk_rows,
          ms_reader, band, selected_buffer, row_data, chunk_data, shared_data);
      break;
    case GainMode::k2VisDiagonal:
      return ReadChunkForInvertImplementation<GainMode::k2VisDiagonal>(
          apply_corrections, task_queue, gridders, ms_data, n_chunk_rows,
          ms_reader, band, selected_buffer, row_data, chunk_data, shared_data);
      break;
    case GainMode::kTrace:
      return ReadChunkForInvertImplementation<GainMode::kTrace>(
          apply_corrections, task_queue, gridders, ms_data, n_chunk_rows,
          ms_reader, band, selected_buffer, row_data, chunk_data, shared_data);
      break;
    case GainMode::kFull:
      return ReadChunkForInvertImplementation<GainMode::kFull>(
          apply_corrections, task_queue, gridders, ms_data, n_chunk_rows,
          ms_reader, band, selected_buffer, row_data, chunk_data, shared_data);
      break;
  }
  assert(false);
  return 0;
}

template <GainMode Mode>
size_t MSGridderManager::ReadChunkForInvertImplementation(
    bool apply_corrections,
    aocommon::TaskQueue<std::function<void()>>& task_queue,
    const std::vector<MsGridder*>& gridders,
    MsProviderCollection::MsData& ms_data, size_t n_chunk_rows,
    MSReader& ms_reader, const aocommon::BandData band,
    const bool* selected_buffer, RowData& row_data, ChunkData& chunk_data,
    MsGridderData& shared_data) {
  size_t n_chunk_rows_read = 0;
  MSProvider::MetaData metadata;
  std::pair<size_t, size_t>* antennas = chunk_data.antennas.data();
  double* uvw = chunk_data.uvw.data();
  std::complex<float>* visibilities = chunk_data.visibilities.data();
  while (ms_reader.CurrentRowAvailable() && n_chunk_rows_read < n_chunk_rows) {
    ms_reader.ReadMeta(metadata);
    uvw[0] = metadata.uInM;
    uvw[1] = metadata.vInM;
    uvw[2] = metadata.wInM;

    // Read and store all visibilities and weights, we need them all in
    // memory when calling 'InlineApplyWeightsAndCorrections' so that we
    // can calculate the final visibilities to return correctly
    shared_data.ReadVisibilities(ms_reader, visibilities, row_data.weights,
                                 row_data.model);
    shared_data.CalculateWeights(uvw, visibilities, band, row_data.weights,
                                 row_data.model, selected_buffer);
    if (shared_data.StoreImagingWeights())
      ms_reader.WriteImagingWeights(shared_data.scratch_image_weights_.data());

    // Sum the corrections and apply the weights.
    // We store the appropriate time_offset to be used later along with other
    // required info when we apply the corrections
    if (apply_corrections) {
      *antennas = std::make_pair(metadata.antenna1, metadata.antenna2);

      size_t time_offset;
      ExecuteForAllGridders(task_queue, [&](MsGridder* gridder) {
        time_offset = chunk_data.time_offsets.back();
        gridder->ApplyCorrections<Mode, ModifierBehaviour::kSum, true>(
            ms_data.antenna_names.size(), visibilities, band, row_data.weights,
            metadata.time, metadata.fieldId, metadata.antenna1,
            metadata.antenna2, time_offset,
            shared_data.scratch_image_weights_.data());
      });
      chunk_data.time_offsets.emplace_back(time_offset);
      ++antennas;
    }
    shared_data.ApplyWeights<Mode>(visibilities, band.ChannelCount(),
                                   row_data.weights);

    // If we aren't applying corrections then we won't have a callback that will
    // later collapse visibilities. As a result we need to collapse the
    // visibilities now. This also allows to use less memory in this scenario so
    // acts as an optimization.
    if (!apply_corrections) {
      if (ms_data.ms_provider->NPolarizations() == 2) {
        internal::CollapseData<2>(band.ChannelCount(), visibilities,
                                  shared_data.Polarization());
      } else if (ms_data.ms_provider->NPolarizations() == 4) {
        internal::CollapseData<4>(band.ChannelCount(), visibilities,
                                  shared_data.Polarization());
      }
      visibilities += band.ChannelCount();
    } else {
      visibilities +=
          band.ChannelCount() * ms_data.ms_provider->NPolarizations();
    }
    uvw += 3;

    ++n_chunk_rows_read;
    ms_reader.NextInputRow();
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

void MSGridderManager::GridChunk(
    size_t n_rows, size_t num_parallel_gridders, bool apply_corrections,
    ChunkData& chunk_data, std::vector<MsGridder*>& gridders,
    size_t gridded_visibility_count, size_t visibility_weight_sum,
    size_t max_gridded_weight, size_t total_weight, size_t n_vis_polarizations,
    aocommon::TaskQueue<std::function<void()>>& task_queue,
    const aocommon::UVector<double>& frequencies,
    const aocommon::BandData& band, MsProviderCollection::MsData& ms_data) {
  // Grid the chunk
  Logger::Info << "Gridding " + std::to_string(n_rows) + " rows for " +
                      std::to_string(gridders.size()) + " facets using " +
                      std::to_string(num_parallel_gridders) + " threads...\n";
  Logger::Info.Flush();

  ExecuteForAllGriddersWithNCores(
      task_queue, num_parallel_gridders,
      [&](MsGridder* gridder, size_t facet_index) {
        Logger::Info << "Gridding facet " + std::to_string(facet_index) + "\n";
        Logger::Info.Flush();

        gridder->gridded_visibility_count_ = gridded_visibility_count;
        gridder->visibility_weight_sum_ = visibility_weight_sum;
        gridder->max_gridded_weight_ = max_gridded_weight;
        gridder->total_weight_ = total_weight;

        gridder->GridSharedMeasurementSetChunk(
            apply_corrections, n_vis_polarizations, n_rows,
            chunk_data.uvw.data(), frequencies.data(), band,
            chunk_data.antennas.data(), chunk_data.visibilities.data(),
            apply_corrections ? chunk_data.time_offsets.data() + 1 : nullptr,
            ms_data.antenna_names.size());
        Logger::Info << "Done gridding facet " + std::to_string(facet_index) +
                            "\n";
        Logger::Info.Flush();
      });

  Logger::Info << "Finished gridding " + std::to_string(n_rows) + " rows for " +
                      std::to_string(gridders.size()) + " facets using " +
                      std::to_string(num_parallel_gridders) + " threads...\n";
  Logger::Info.Flush();
  ms_data.total_rows_processed += n_rows;
}

void MSGridderManager::ReadChunksForInvert(
    aocommon::Lane<ChunkData>& task_lane,
    aocommon::TaskQueue<std::function<void()>>& task_queue,
    size_t n_max_rows_in_memory, bool apply_corrections,
    MsProviderCollection::MsData& ms_data, MsGridderData& shared_data,
    const std::vector<MsGridder*>& gridders, const aocommon::BandData band,
    size_t n_vis_polarizations, const bool* selected_buffer) {
  // Row data
  const size_t data_size = band.ChannelCount() * n_vis_polarizations;
  aocommon::UVector<std::complex<float>> model_buffer(data_size);
  aocommon::UVector<float> weight_buffer(data_size);
  RowData row_data;
  row_data.model = model_buffer.data();
  row_data.weights = weight_buffer.data();

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
    ChunkData chunk_data(target_chunk_size, band.ChannelCount(),
                         n_vis_polarizations, apply_corrections);
    if (apply_corrections) {
      chunk_data.time_offsets.reserve(target_chunk_size + 1);
      chunk_data.time_offsets.push_back(0);
    }

    const size_t n_rows = ReadChunkForInvert(
        shared_data.GetGainMode(), apply_corrections, task_queue, gridders,
        ms_data, target_chunk_size, *ms_reader, band, selected_buffer, row_data,
        chunk_data, shared_data);

    chunk_data.gridded_visibility_count = shared_data.gridded_visibility_count_;
    chunk_data.visibility_weight_sum = shared_data.visibility_weight_sum_;
    chunk_data.max_gridded_weight = shared_data.max_gridded_weight_;
    chunk_data.total_weight = shared_data.total_weight_;
    chunk_data.n_rows = n_rows;
    Logger::Debug << "Done loading chunk " << chunk_index << ".\n";

    task_lane.write(chunk_data);
    ++chunk_index;
    target_chunk_size = n_rows_per_chunk;
  }
  Logger::Info << "All gridding rows loaded...\n";
  task_lane.write_end();
}

void MSGridderManager::GridChunks(
    aocommon::Lane<ChunkData>& task_lane, const size_t num_parallel_gridders,
    const bool apply_corrections, std::vector<MsGridder*>& gridders,
    aocommon::TaskQueue<std::function<void()>>& task_queue,
    const aocommon::UVector<double>& frequencies,
    const aocommon::BandData& band, MsProviderCollection::MsData& ms_data,
    size_t n_vis_polarizations) {
  ChunkData chunk_data;
  size_t chunk_index = 1;
  while (task_lane.read(chunk_data)) {
    Logger::Info << "Gridding chunk" << chunk_index << ".\n";
    GridChunk(chunk_data.n_rows, num_parallel_gridders, apply_corrections,
              chunk_data, gridders, chunk_data.gridded_visibility_count,
              chunk_data.visibility_weight_sum, chunk_data.max_gridded_weight,
              chunk_data.total_weight, n_vis_polarizations, task_queue,
              frequencies, band, ms_data);
    Logger::Info << "Done gridding chunk" << chunk_index << ".\n";
    ++chunk_index;
  }
  Logger::Info << "All gridding rows processed...\n";
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

  // NB! This delibritely leads to overallocation of threads
  // As this still outperforms the alternative of not overlapping the IO.
  // Future changes will implement work stealing which should
  // fix this overallocation.
  aocommon::TaskQueue<std::function<void()>> task_queue_read;
  std::vector<std::thread> thread_pool_read;
  thread_pool_read.reserve(available_cores_);
  for (size_t i = 0; i < available_cores_; ++i) {
    thread_pool_read.emplace_back([&] {
      std::function<void()> operation;
      while (task_queue_read.Pop(operation)) {
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
      aocommon::Lane<ChunkData> task_lane(1);
      std::thread grid_chunks_thread([&] {
        GridChunks(task_lane, num_parallel_gridders, apply_corrections,
                   gridders, task_queue, frequencies, band, ms_data,
                   n_vis_polarizations);
      });
      ReadChunksForInvert(task_lane, task_queue_read, n_max_rows_in_memory,
                          apply_corrections, ms_data, shared_data, gridders,
                          band, n_vis_polarizations, selected_buffer.data());
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

  // Clean up the thread pool
  task_queue_read.Finish();
  for (std::thread& thread : thread_pool_read) {
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
