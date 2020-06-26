#ifndef META_DATA_CACHE_H
#define META_DATA_CACHE_H

#include <memory>
#include <vector>

#include "../idg/averagebeam.h"

struct MetaDataCache
{
	struct Entry {
		double minW, maxW, maxWWithFlags, maxBaselineUVW, maxBaselineInM, integrationTime;
	};
	std::vector<Entry> msDataVector;
	std::unique_ptr<class AverageBeam> averageBeam;
	
	void Serialize(std::ostream& stream) const;
	void Unserialize(std::istream& stream);
};


#endif
