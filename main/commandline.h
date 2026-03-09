#ifndef WSCLEAN_COMMAND_LINE_H
#define WSCLEAN_COMMAND_LINE_H

#include <string>
#include <cstring>

#include "settings.h"

namespace wsclean {

class WSClean;

class CommandLine {
 public:
  /**
   * Initialize wsclean from the given command line options
   * @returns @c true when parameters indicate wsclean should
   * be run, @c false e.g. when this is a dry run.
   */
  static bool Parse(WSClean& wsclean, int argc, const char* argv[],
                    bool is_slave, bool use_in_memory_data) {
    bool fullRun = ParseWithoutValidation(wsclean, argc, argv, is_slave,
                                          use_in_memory_data);
    if (fullRun) Validate(wsclean);
    return fullRun;
  }

  /**
   * Initialize wsclean from the given command line options, but
   * do only limitted validation of the parameters. This is useful
   * when the logging needs to be set up before further validation
   * is done, as is used in distributed/wsclean-mp.cpp.
   * Otherwise similar to @ref Parse().
   */
  static bool ParseWithoutValidation(WSClean& wsclean, int argc,
                                     const char* argv[], bool is_slave,
                                     bool use_in_memory_data);
  /**
   * Finish a call to ParseWithoutValidation().
   */
  static void Validate(WSClean& wsclean);

  static void Run(WSClean& wsclean);
};

}  // namespace wsclean

#endif
