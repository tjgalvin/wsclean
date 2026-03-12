#ifndef WSCLEAN_IO_READ_IN_MEMORY_H_
#define WSCLEAN_IO_READ_IN_MEMORY_H_

#include "../structures/inmemorypart.h"

#include "../interface/inmemoryms.h"

#include <string>

namespace schaapcommon::reordering {
class MSSelection;
}  // namespace schaapcommon::reordering

namespace wsclean {

InMemoryMs ReadMsInMemory(const std::string& measurement_set, const std::string& data_column_name,
                          const schaapcommon::reordering::MSSelection& selection,
                          bool sort_in_time);

/**
 * Removes rows that are not selected. If @p selection is left unset, then
 * the function returns immediately. The channel range selection is not
 * applied, because this is done during the reordering step.
 */
void ApplyRowSelectionToMs(InMemoryMs& ms, const schaapcommon::reordering::MSSelection& selection);

}  // namespace wsclean

#endif
