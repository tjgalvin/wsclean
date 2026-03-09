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

/**
 * Perform a WSClean imaging run on in-memory data. The @p
 * command_line_parameters should specify the settings of the run, and should be
 * similar to a normal wsclean command line parameter list. The parameters
 * should not include the executable ("wsclean") and should not include
 * measurement set filenames.
 *
 * By default, WSClean will produce logging output on the command line. If this
 * is undesirable, the option "-quiet" can be added to surpress this output.
 */
void Image(const std::string& command_line_parameters,
           std::vector<InMemoryMs>&& ms_list);

}  // namespace wsclean

#endif
