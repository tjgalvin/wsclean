#ifndef WSCLEAN_COMMAND_LINE_H
#define WSCLEAN_COMMAND_LINE_H

#include <string>
#include <cstring>

#include "settings.h"

class CommandLine {
 public:
  /**
   * Initialize wsclean from the given command line options
   * @returns @c true when parameters indicate wsclean should
   * be run, @c false e.g. when this is a dry run.
   */
  static bool Parse(class WSClean& wsclean, int argc, const char* argv[],
                    bool isSlave) {
    bool fullRun = ParseWithoutValidation(wsclean, argc, argv, isSlave);
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
  static bool ParseWithoutValidation(class WSClean& wsclean, int argc,
                                     const char* argv[], bool isSlave);
  /**
   * Finish a call to ParseWithoutValidation().
   */
  static void Validate(class WSClean& wsclean);

  static void Run(class WSClean& wsclean);

 private:
  static void deprecated(bool isSlave, const std::string& param,
                         const std::string& replacement);
  static void printHeader();
  static void printHelp();
  static std::vector<std::string> parseStringList(const char* param);
  static size_t parse_size_t(const char* param, const char* name);
  static double parse_double(const char* param, double lowerLimit,
                             const char* name, bool inclusive = true);
  static double parse_double(const char* param, const char* name);
};

#endif
