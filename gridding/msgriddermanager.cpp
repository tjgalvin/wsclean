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
    MsProviderCollection::MsData& ms_data, size_t max_n_rows,
    MSReader& ms_reader, const aocommon::BandData band,
    const bool* selected_buffer, RowData& row_data, ChunkData& chunk_data,
    MsGridderData& shared_data) {
  switch (gain_mode) {
    case GainMode::kXX:
      return ReadChunkForInvertImplementation<GainMode::kXX>(
          apply_corrections, task_queue, gridders, ms_data, max_n_rows,
          ms_reader, band, selected_buffer, row_data, chunk_data, shared_data);
      break;
    case GainMode::kYY:
      return ReadChunkForInvertImplementation<GainMode::kYY>(
          apply_corrections, task_queue, gridders, ms_data, max_n_rows,
          ms_reader, band, selected_buffer, row_data, chunk_data, shared_data);
      break;
    case GainMode::k2VisDiagonal:
      return ReadChunkForInvertImplementation<GainMode::k2VisDiagonal>(
          apply_corrections, task_queue, gridders, ms_data, max_n_rows,
          ms_reader, band, selected_buffer, row_data, chunk_data, shared_data);
      break;
    case GainMode::kTrace:
      return ReadChunkForInvertImplementation<GainMode::kTrace>(
          apply_corrections, task_queue, gridders, ms_data, max_n_rows,
          ms_reader, band, selected_buffer, row_data, chunk_data, shared_data);
      break;
    case GainMode::kFull:
      return ReadChunkForInvertImplementation<GainMode::kFull>(
          apply_corrections, task_queue, gridders, ms_data, max_n_rows,
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
    MsProviderCollection::MsData& ms_data, size_t max_n_rows,
    MSReader& ms_reader, const aocommon::BandData band,
    const bool* selected_buffer, RowData& row_data, ChunkData& chunk_data,
    MsGridderData& shared_data) {
  // Clear 'per chunk' data that we have cached so that they are fresh for the
  // next chunk of data
  if (apply_corrections) {
    for (MsGridder* gridder : gridders) {
      gridder->time_offsets_.clear();
      gridder->time_offsets_[ms_data.internal_ms_index].push_back(0);
    }
  }
  size_t n_chunk_rows_read = 0;
  MSProvider::MetaData metadata;
  while (ms_reader.CurrentRowAvailable() && n_chunk_rows_read < max_n_rows) {
    ms_reader.ReadMeta(metadata);
    chunk_data.uvw[0] = metadata.uInM;
    chunk_data.uvw[1] = metadata.vInM;
    chunk_data.uvw[2] = metadata.wInM;

    // Read and store all visibilities and weights, we need them all in
    // memory when calling 'InlineApplyWeightsAndCorrections' so that we
    // can calculate the final visibilities to return correctly
    shared_data.ReadVisibilities(ms_reader, chunk_data.visibilities,
                                 row_data.weights, row_data.model);
    shared_data.CalculateWeights(chunk_data.uvw, chunk_data.visibilities, band,
                                 row_data.weights, row_data.model,
                                 selected_buffer);
    if (shared_data.StoreImagingWeights())
      ms_reader.WriteImagingWeights(shared_data.scratch_image_weights_.data());

    // Sum the corrections and apply the weights.
    // We store the appropriate time_offset to be used later along with other
    // required info when we apply the corrections
    if (apply_corrections) {
      *chunk_data.antennas =
          std::make_pair(metadata.antenna1, metadata.antenna2);

      ExecuteForAllGridders(task_queue, [&](MsGridder* gridder) {
        // TODO: Do we need time offset for every gridder or can we hold it once
        // on the shared weight manager?
        std::vector<size_t>& time_offsets =
            gridder->time_offsets_[ms_data.internal_ms_index];
        size_t time_offset = time_offsets.back();
        gridder->ApplyCorrections<Mode, ModifierBehaviour::kSum, true>(
            ms_data.antenna_names.size(), chunk_data.visibilities, band,
            row_data.weights, metadata.time, metadata.fieldId,
            metadata.antenna1, metadata.antenna2, time_offset,
            shared_data.scratch_image_weights_.data());
        time_offsets.emplace_back(time_offset);
      });
      ++chunk_data.antennas;
    }
    shared_data.ApplyWeights<Mode>(chunk_data.visibilities, band.ChannelCount(),
                                   row_data.weights);

    // If we aren't applying corrections then we won't have a callback that will
    // later collapse visibilities. As a result we need to collapse the
    // visibilities now. This also allows to use less memory in this scenario so
    // acts as an optimization.
    if (!apply_corrections) {
      if (ms_data.ms_provider->NPolarizations() == 2) {
        internal::CollapseData<2>(band.ChannelCount(), chunk_data.visibilities,
                                  shared_data.Polarization());
      } else if (ms_data.ms_provider->NPolarizations() == 4) {
        internal::CollapseData<4>(band.ChannelCount(), chunk_data.visibilities,
                                  shared_data.Polarization());
      }
      chunk_data.visibilities += band.ChannelCount();
    } else {
      chunk_data.visibilities +=
          band.ChannelCount() * ms_data.ms_provider->NPolarizations();
    }
    chunk_data.uvw += 3;

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
      const size_t n_vis_polarizations = ms_data.ms_provider->NPolarizations();
      const size_t data_size = band.ChannelCount() * n_vis_polarizations;

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
        // For each row we have to store an antenna pair and a per gridder
        // solution time offset
        additional_per_vis_mem = sizeof(size_t) * (gridders.size() + 2);
      }

      const size_t n_max_chunk_rows = gridders[0]->CalculateMaxRowsInMemory(
          available_memory_, constant_mem, additional_per_vis_mem,
          band.ChannelCount(), apply_corrections ? 1 : n_vis_polarizations);

      if (apply_corrections) {
        for (MsGridder* gridder : gridders) {
          gridder->time_offsets_[ms_data.internal_ms_index].reserve(
              n_max_chunk_rows + 1);
        }
      }

      aocommon::UVector<double> frequencies(band.ChannelCount());
      for (size_t i = 0; i != frequencies.size(); ++i)
        frequencies[i] = band.ChannelFrequency(i);

      aocommon::UVector<bool> selected_buffer(band.ChannelCount(), true);

      // Row data
      aocommon::UVector<std::complex<float>> model_buffer(data_size);
      aocommon::UVector<float> weight_buffer(data_size);
      RowData row_data;
      row_data.model = model_buffer.data();
      row_data.weights = weight_buffer.data();

      // Chunk data
      aocommon::UVector<std::pair<size_t, size_t>> antennas(
          apply_corrections ? n_max_chunk_rows : 0);
      aocommon::UVector<double> uvw_buffer(n_max_chunk_rows * 3);
      // If we don't apply corrections then we collapse the visibilities when
      // storing them to save memory. In order to be able to collapse the
      // already copied data in place we have to slightly overallocate the
      // buffer by one uncollapsed row.
      const size_t visibility_size =
          apply_corrections
              ? (n_max_chunk_rows * data_size)
              : ((n_max_chunk_rows * band.ChannelCount()) + data_size);
      aocommon::UVector<std::complex<float>> visibilities(visibility_size);
      ChunkData chunk_data;

      // Iterate over chunks until all data has been gridded
      std::unique_ptr<MSReader> ms_reader = ms_data.ms_provider->MakeReader();
      while (ms_reader->CurrentRowAvailable()) {
        Logger::Debug << "Max " << n_max_chunk_rows << " rows fit in memory.\n";
        Logger::Info << "Loading data in memory...\n";

        // Read / fill the chunk
        chunk_data.antennas = antennas.data();
        chunk_data.uvw = uvw_buffer.data();
        chunk_data.visibilities = visibilities.data();
        size_t n_rows = ReadChunkForInvert(
            shared_data.GetGainMode(), apply_corrections, task_queue, gridders,
            ms_data, n_max_chunk_rows, *ms_reader, band, selected_buffer.data(),
            row_data, chunk_data, shared_data);

        // Grid the chunk
        Logger::Info << "Gridding " + std::to_string(n_rows) + " rows for " +
                            std::to_string(gridders.size()) + " facets using " +
                            std::to_string(num_parallel_gridders) +
                            " threads...\n";
        ExecuteForAllGriddersWithNCores(
            task_queue, num_parallel_gridders,
            [&](MsGridder* gridder, size_t facet_index) {
              Logger::Info << "Gridding facet " + std::to_string(facet_index) +
                                  "\n";

              gridder->gridded_visibility_count_ =
                  shared_data.gridded_visibility_count_;
              gridder->visibility_weight_sum_ =
                  shared_data.visibility_weight_sum_;
              gridder->max_gridded_weight_ = shared_data.max_gridded_weight_;
              gridder->total_weight_ = shared_data.total_weight_;

              gridder->GridSharedMeasurementSetChunk(
                  apply_corrections, shared_data.n_vis_polarizations_, n_rows,
                  uvw_buffer.data(), frequencies.data(), band, antennas.data(),
                  visibilities.data(),
                  apply_corrections
                      ? gridder->time_offsets_[ms_data.internal_ms_index]
                                .data() +
                            1
                      : nullptr,
                  ms_data.antenna_names.size());
              Logger::Info << "Done gridding facet " +
                                  std::to_string(facet_index) + "\n";
            });
        ms_data.total_rows_processed += n_rows;
      }
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
