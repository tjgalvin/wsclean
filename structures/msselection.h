#ifndef MS_SELECTION
#define MS_SELECTION

#include "imagingtableentry.h"

#include <schaapcommon/reordering/msselection.h>
#include <aocommon/multibanddata.h>

namespace wsclean {

/**
 * Change this selection object so that its datadescid and channel range
 * correspond with the given entry. If the specified bands are not necessary
 * for this entry, the msselection is not changed and the function returns
 * false.
 */
bool SelectMsChannels(schaapcommon::reordering::MSSelection& selection,
                      const aocommon::MultiBandData& msBands, size_t dataDescId,
                      const ImagingTableEntry& entry);

}  // namespace wsclean

#endif
