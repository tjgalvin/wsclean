#include "commandline.h"
#include "wsclean.h"

#include <aocommon/checkblas.h>
#include <aocommon/logger.h>

#include <exception>
#include <iostream>

using aocommon::Logger;

int main(int argc, char* argv[]) {
  try {
    check_openblas_multithreading();
    WSClean wsclean;
    if (CommandLine::Parse(wsclean, argc, const_cast<const char**>(argv),
                           false))
      CommandLine::Run(wsclean);
    return 0;
  } catch (std::exception& e) {
    Logger::Error << "+ + + + + + + + + + + + + + + + + + +\n"
                  << "+ An exception occured:\n";
    std::istringstream iss(e.what());
    for (std::string line; std::getline(iss, line);) {
      Logger::Error << "+ >>> " << line << "\n";
    }
    Logger::Error << "+ + + + + + + + + + + + + + + + + + +\n";
    return -1;
  }
}
