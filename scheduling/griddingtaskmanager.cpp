#include "griddingtaskmanager.h"

#include "mpischeduler.h"
#include "threadedscheduler.h"

#include "../main/settings.h"

#include "../gridding/wsmsgridder.h"
#include "../gridding/directmsgridder.h"

#include "../idg/idgmsgridder.h"

#ifdef HAVE_WGRIDDER
#include "../wgridder/wgriddingmsgridder.h"
#endif

GriddingTaskManager::GriddingTaskManager(const class Settings& settings)
    : _settings(settings) {}

GriddingTaskManager::~GriddingTaskManager() {}

std::unique_ptr<GriddingTaskManager> GriddingTaskManager::Make(
    const class Settings& settings, bool useDirectScheduler) {
  if (settings.useMPI && !useDirectScheduler) {
#ifdef HAVE_MPI
    return std::unique_ptr<GriddingTaskManager>(new MPIScheduler(settings));
#else
    throw std::runtime_error("MPI not available");
#endif
  } else if (settings.parallelGridding == 1 || useDirectScheduler) {
    return std::unique_ptr<GriddingTaskManager>(
        new GriddingTaskManager(settings));
  } else {
    return std::unique_ptr<GriddingTaskManager>(
        new ThreadedScheduler(settings));
  }
}

std::unique_ptr<MSGridderBase> GriddingTaskManager::createGridder() const {
  if (_settings.useIDG) {
    return std::unique_ptr<MSGridderBase>(new IdgMsGridder(_settings));
  } else if (_settings.useWGridder) {
#ifdef HAVE_WGRIDDER
    return std::unique_ptr<MSGridderBase>(new WGriddingMSGridder(
        _settings.threadCount, _settings.memFraction, _settings.absMemLimit,
        _settings.wgridderAccuracy));
#else
    throw std::runtime_error(
        "WGridder cannot be used: WGridder requires a C++17 compiler, which "
        "was not found during compilation. Update your compiler and recompile "
        "wsclean.");
#endif
  } else if (_settings.directFT) {
    switch (_settings.directFTPrecision) {
      case DirectFTPrecision::Half:
        throw std::runtime_error("Half precision is not implemented");
        // return std::unique_ptr<MSGridderBase>(new
        // DirectMSGridder<half_float::half>(&_imageAllocator,
        // _settings.threadCount));
        break;
      case DirectFTPrecision::Float:
        return std::unique_ptr<MSGridderBase>(
            new DirectMSGridder<float>(_settings.threadCount));
        break;
      default:
      case DirectFTPrecision::Double:
        return std::unique_ptr<MSGridderBase>(
            new DirectMSGridder<double>(_settings.threadCount));
        break;
      case DirectFTPrecision::LongDouble:
        return std::unique_ptr<MSGridderBase>(
            new DirectMSGridder<long double>(_settings.threadCount));
        break;
    }
  } else
    return std::unique_ptr<MSGridderBase>(new WSMSGridder(
        _settings.threadCount, _settings.memFraction, _settings.absMemLimit));
}

void GriddingTaskManager::Run(
    GriddingTask&& task, std::function<void(GriddingResult&)> finishCallback) {
  if (!_gridder) {
    _gridder = createGridder();
    prepareGridder(*_gridder);
  }

  GriddingResult result = runDirect(std::move(task), *_gridder);
  finishCallback(result);
}

GriddingResult GriddingTaskManager::RunDirect(GriddingTask&& task) {
  if (!_gridder) {
    _gridder = createGridder();
    prepareGridder(*_gridder);
  }

  return runDirect(std::move(task), *_gridder);
}

GriddingResult GriddingTaskManager::runDirect(GriddingTask&& task,
                                              MSGridderBase& gridder) {
  gridder.ClearMeasurementSetList();
  std::vector<std::unique_ptr<MSProvider>> msProviders;
  for (auto& p : task.msList) {
    msProviders.emplace_back(p->GetProvider());
    gridder.AddMeasurementSet(msProviders.back().get(), p->Selection());
  }
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
    gridder.SetAddToModel(task.addToModel);
    if (task.polarization == aocommon::Polarization::XY ||
        task.polarization == aocommon::Polarization::YX)
      gridder.Predict(std::move(task.modelImageReal),
                      std::move(task.modelImageImaginary));
    else
      gridder.Predict(std::move(task.modelImageReal));
  }

  GriddingResult result;
  result.imageRealResult = gridder.ImageRealResult();
  if (gridder.IsComplex())
    result.imageImaginaryResult = gridder.ImageImaginaryResult();
  result.startTime = gridder.StartTime();
  result.beamSize = gridder.BeamSize();
  result.imageWeight = gridder.ImageWeight();
  result.normalizationFactor = gridder.NormalizationFactor();
  result.actualWGridSize = gridder.ActualWGridSize();
  result.griddedVisibilityCount = gridder.GriddedVisibilityCount();
  result.effectiveGriddedVisibilityCount =
      gridder.EffectiveGriddedVisibilityCount();
  result.visibilityWeightSum = gridder.VisibilityWeightSum();
  result.actualInversionWidth = gridder.ActualInversionWidth();
  result.actualInversionHeight = gridder.ActualInversionHeight();
  result.cache = gridder.AcquireMetaDataCache();
  return result;
}

void GriddingTaskManager::prepareGridder(MSGridderBase& gridder) {
  gridder.SetGridMode(_settings.gridMode);
  gridder.SetImageWidth(_settings.paddedImageWidth);
  gridder.SetImageHeight(_settings.paddedImageHeight);
  gridder.SetTrimSize(_settings.trimmedImageWidth,
                      _settings.trimmedImageHeight);
  gridder.SetNWSize(_settings.widthForNWCalculation,
                    _settings.heightForNWCalculation);
  gridder.SetNWFactor(_settings.nWLayersFactor);
  gridder.SetPixelSizeX(_settings.pixelScaleX);
  gridder.SetPixelSizeY(_settings.pixelScaleY);
  if (_settings.nWLayers != 0)
    gridder.SetWGridSize(_settings.nWLayers);
  else
    gridder.SetNoWGridSize();
  gridder.SetAntialiasingKernelSize(_settings.antialiasingKernelSize);
  gridder.SetOverSamplingFactor(_settings.overSamplingFactor);
  gridder.SetDataColumnName(_settings.dataColumnName);
  gridder.SetWeighting(_settings.weightMode);
  gridder.SetWLimit(_settings.wLimit / 100.0);
  gridder.SetSmallInversion(_settings.smallInversion);
  gridder.SetVisibilityWeightingMode(_settings.visibilityWeightingMode);
}
