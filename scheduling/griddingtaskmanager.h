#ifndef GRIDDING_OPERATION_H
#define GRIDDING_OPERATION_H

#include "griddingtask.h"
#include "griddingresult.h"

#include "../gridding/measurementsetgridder.h"
#include "../gridding/msgridderbase.h"

#include "../structures/imageweights.h"
#include "../structures/observationinfo.h"
#include "../structures/msselection.h"

#include <aocommon/polarization.h>

#include "../msproviders/msdatadescription.h"

#include <aocommon/lane.h>

#include <cstring>
#include <functional>
#include <vector>

class GriddingTaskManager {
 public:
  virtual ~GriddingTaskManager();

  virtual void Run(GriddingTask& task,
                   std::function<void(GriddingResult&)> finishCallback);

  GriddingResult RunDirect(GriddingTask& task);

  virtual void Finish(){};

  // MSGridderBase* Gridder() { return _gridder.get(); }

  static std::unique_ptr<GriddingTaskManager> Make(
      const class Settings& settings, const struct ObservationInfo& obsInfo,
      bool useDirectScheduler = false);

 protected:
  const class Settings& _settings;

  GriddingTaskManager(const class Settings& settings,
                      const struct ObservationInfo& obsInfo);

  std::unique_ptr<MSGridderBase> createGridder() const;
  void prepareGridder(MSGridderBase& gridder);
  GriddingResult runDirect(GriddingTask& task, MSGridderBase& gridder);

 private:
  std::unique_ptr<MSGridderBase> _gridder;
  const struct ObservationInfo _obsInfo;
};

#endif
