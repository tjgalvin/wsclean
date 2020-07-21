#ifndef GRIDDING_OPERATION_H
#define GRIDDING_OPERATION_H

#include "griddingtask.h"
#include "griddingresult.h"

#include "../wsclean/observationinfo.h"
#include "../wsclean/measurementsetgridder.h"
#include "../wsclean/msgridderbase.h"

#include "../imageweights.h"
#include <aocommon/polarization.h>
#include "../msselection.h"

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
      const class WSCleanSettings& settings, bool useDirectScheduler = false);

 protected:
  const class WSCleanSettings& _settings;

  GriddingTaskManager(const class WSCleanSettings& settings);

  std::unique_ptr<MSGridderBase> createGridder() const;
  void prepareGridder(MSGridderBase& gridder);
  GriddingResult runDirect(GriddingTask& task, MSGridderBase& gridder);

 private:
  std::unique_ptr<MSGridderBase> _gridder;
};

#endif
