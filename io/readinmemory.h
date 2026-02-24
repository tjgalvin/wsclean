#ifndef WSCLEAN_IO_READ_IN_MEMORY_H_
#define WSCLEAN_IO_READ_IN_MEMORY_H_

#include "../structures/inmemorypart.h"

#include "../interface/inmemoryms.h"

#include <string>

namespace schaapcommon::reordering {
class MSSelection;
}  // namespace schaapcommon::reordering

namespace wsclean {

InMemoryMs ReadInMemory(const std::string& measurement_set,
                        const std::string& data_column_name,
                        const schaapcommon::reordering::MSSelection& selection);

}

#endif
