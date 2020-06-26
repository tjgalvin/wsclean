#ifndef GRIDDING_TASK_H
#define GRIDDING_TASK_H

#include "../polarization.h"
#include "../image.h"
#include "../imageweights.h"

#include "../msproviders/msdatadescription.h"

#include "metadatacache.h"

#include <istream>
#include <ostream>

class GriddingTask
{
public:
	enum Operation { Invert, Predict } operation;
	bool imagePSF;
	bool subtractModel;
	PolarizationEnum polarization;
	bool verbose;
	std::unique_ptr<MetaDataCache> cache;
	bool storeImagingWeights;
	
	std::shared_ptr<ImageWeights> imageWeights;
	std::vector<std::unique_ptr<MSDataDescription>> msList;
	
	// For prediction
	bool addToModel;
	Image modelImageReal;
	Image modelImageImaginary;
	
	void Serialize(std::ostream& stream) const;
	void Unserialize(std::istream& stream);
};

#endif
