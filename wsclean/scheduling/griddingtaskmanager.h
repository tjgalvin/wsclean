#ifndef GRIDDING_OPERATION_H
#define GRIDDING_OPERATION_H

#include <cstring>
#include <functional>
#include <vector>

#include "griddingtask.h"
#include "griddingresult.h"

#include "../wsclean/observationinfo.h"
#include "../wsclean/measurementsetgridder.h"
#include "../wsclean/msgridderbase.h"

#include "../lane.h"
#include "../imageweights.h"
#include "../polarization.h"
#include "../msselection.h"

#include "../msproviders/msdatadescription.h"

class GriddingTaskManager
{
public:
	virtual ~GriddingTaskManager();
	
	virtual void Run(GriddingTask& task, std::function<void(GriddingResult&)> finishCallback);
	
	virtual void Finish() { };
	
	//MSGridderBase* Gridder() { return _gridder.get(); }
	
	static std::unique_ptr<GriddingTaskManager> Make(const class WSCleanSettings& settings, bool useDirectScheduler = false);
	GriddingResult RunDirect(GriddingTask& task)
	{
		if(!_gridder) {
			_gridder = createGridder();
			prepareGridder(*_gridder);
		}
	
		return runDirect(task, *_gridder);
	}
	
protected:
	const class WSCleanSettings& _settings;
	
	GriddingTaskManager(const class WSCleanSettings& settings);
	
	std::unique_ptr<MSGridderBase> createGridder() const;
	void prepareGridder (MSGridderBase& gridder);
	GriddingResult runDirect(GriddingTask& task, MSGridderBase& gridder);
	
private:
	std::unique_ptr<MSGridderBase> _gridder;
};

#endif
