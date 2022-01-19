#include "griddingtaskmanager.h"

#include "mpischeduler.h"
#include "threadedscheduler.h"

#include "../main/settings.h"

#include "../gridding/msgridderbase.h"
#include "../gridding/wsmsgridder.h"
#include "../gridding/directmsgridder.h"

#include "../idg/idgmsgridder.h"

#include <schaapcommon/facets/facet.h>

#ifdef HAVE_WGRIDDER
#include "../wgridder/wgriddingmsgridder.h"
#endif

GriddingTaskManager::GriddingTaskManager(const class Settings& settings)
    : _settings(settings) {}

GriddingTaskManager::~GriddingTaskManager() {}

std::unique_ptr<GriddingTaskManager> GriddingTaskManager::Make(
    const class Settings& settings) {
  if (settings.useMPI) {
#ifdef HAVE_MPI
    return std::unique_ptr<GriddingTaskManager>(new MPIScheduler(settings));
#else
    throw std::runtime_error("MPI not available");
#endif
  } else if (settings.parallelGridding == 1) {
    return std::unique_ptr<GriddingTaskManager>(
        new GriddingTaskManager(settings));
  } else {
    return std::unique_ptr<GriddingTaskManager>(
        new ThreadedScheduler(settings));
  }
}

void GriddingTaskManager::Run(
    GriddingTask&& task, std::function<void(GriddingResult&)> finishCallback) {
  GriddingResult result = RunDirect(std::move(task));
  finishCallback(result);
}

GriddingResult GriddingTaskManager::RunDirect(GriddingTask&& task) {
  std::unique_ptr<MSGridderBase> gridder(makeGridder());
  return runDirect(std::move(task), *gridder);
}

GriddingResult GriddingTaskManager::runDirect(GriddingTask&& task,
                                              MSGridderBase& gridder) {
  gridder.ClearMeasurementSetList();
  std::vector<std::unique_ptr<MSProvider>> msProviders;
  for (auto& p : task.msList) {
    msProviders.emplace_back(p->GetProvider());
    gridder.AddMeasurementSet(msProviders.back().get(), p->Selection());
  }

  gridder.SetFacetGroupIndex(task.facetGroupIndex);
  gridder.SetAdditivePredict(task.facet != nullptr);
  if (task.facet != nullptr) {
    gridder.SetFacetIndex(task.facetIndex);
    gridder.SetImageWidth(task.facet->GetUntrimmedBoundingBox().Width());
    gridder.SetImageHeight(task.facet->GetUntrimmedBoundingBox().Height());
    gridder.SetTrimSize(task.facet->GetTrimmedBoundingBox().Width(),
                        task.facet->GetTrimmedBoundingBox().Height());
    gridder.SetFacetDirectionRA(task.facet->RA());
    gridder.SetFacetDirectionDec(task.facet->Dec());
  } else {
    gridder.SetImageWidth(_settings.paddedImageWidth);
    gridder.SetImageHeight(_settings.paddedImageHeight);
    gridder.SetTrimSize(_settings.trimmedImageWidth,
                        _settings.trimmedImageHeight);
  }
  gridder.SetImagePadding(_settings.imagePadding);
  gridder.SetPhaseCentreDec(task.observationInfo.phaseCentreDec);
  gridder.SetPhaseCentreRA(task.observationInfo.phaseCentreRA);
  gridder.SetPhaseCentreDM(task.observationInfo.shiftM);
  gridder.SetPhaseCentreDL(task.observationInfo.shiftL);
  gridder.SetPolarization(task.polarization);
  gridder.SetIsComplex(task.polarization == aocommon::Polarization::XY ||
                       task.polarization == aocommon::Polarization::YX);
  gridder.SetIsFirstIteration(task.verbose);
  if (task.cache)
    gridder.SetMetaDataCache(std::move(task.cache));
  else
    gridder.SetMetaDataCache(
        std::unique_ptr<MetaDataCache>(new MetaDataCache()));
  gridder.SetImageWeights(task.imageWeights.get());
  if (task.operation == GriddingTask::Invert) {
    gridder.SetDoImagePSF(task.imagePSF);
    gridder.SetDoSubtractModel(task.subtractModel);
    gridder.SetStoreImagingWeights(task.storeImagingWeights);
    gridder.Invert();
  } else {
    gridder.SetWriterLockManager(this);
    gridder.Predict(std::move(task.modelImages));
  }

  GriddingResult result;
  result.images = gridder.ResultImages();
  result.startTime = gridder.StartTime();
  result.beamSize = gridder.BeamSize();
  result.imageWeight = gridder.ImageWeight();
  result.normalizationFactor = gridder.NormalizationFactor();
  result.actualWGridSize = gridder.ActualWGridSize();
  result.griddedVisibilityCount = gridder.GriddedVisibilityCount();
  result.effectiveGriddedVisibilityCount =
      gridder.EffectiveGriddedVisibilityCount();
  result.visibilityWeightSum = gridder.VisibilityWeightSum();
  result.cache = gridder.AcquireMetaDataCache();
  return result;
}

std::unique_ptr<MSGridderBase> GriddingTaskManager::constructGridder() const {
  if (_settings.useIDG) {
    return std::unique_ptr<MSGridderBase>(new IdgMsGridder(_settings));
  } else if (_settings.useWGridder) {
#ifdef HAVE_WGRIDDER
    return std::unique_ptr<MSGridderBase>(new WGriddingMSGridder(_settings));
#else
    throw std::runtime_error(
        "WGridder cannot be used: WGridder requires a C++17 compiler, which "
        "was not found during compilation. Update your compiler and recompile "
        "wsclean.");
#endif
  } else if (_settings.directFT) {
    switch (_settings.directFTPrecision) {
      case DirectFTPrecision::Float:
        return std::unique_ptr<MSGridderBase>(
            new DirectMSGridder<float>(_settings));
        break;
      default:
      case DirectFTPrecision::Double:
        return std::unique_ptr<MSGridderBase>(
            new DirectMSGridder<double>(_settings));
        break;
      case DirectFTPrecision::LongDouble:
        return std::unique_ptr<MSGridderBase>(
            new DirectMSGridder<long double>(_settings));
        break;
    }
  } else
    return std::unique_ptr<MSGridderBase>(new WSMSGridder(_settings));
}

std::unique_ptr<MSGridderBase> GriddingTaskManager::makeGridder() const {
  std::unique_ptr<MSGridderBase> gridder(constructGridder());
  gridder->SetGridMode(_settings.gridMode);
  return gridder;
}
