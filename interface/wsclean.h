#ifndef WSCLEAN_INTERFACE_WSCLEAN_H_
#define WSCLEAN_INTERFACE_WSCLEAN_H_

/**
 * This is the public header file for WSClean, that other programs
 * can use to call WSClean operations from C++ applications.
 */

#include <string>

// The following list of dependencies (including transitive ones) should be kept
// minimal: no files outside the interface/ directory should be included,
// because they are not installed.
#include "inmemoryms.h"

namespace wsclean {

void Image(const std::string& command_line_parameters, InMemoryMs&& ms);

}  // namespace wsclean

#endif
