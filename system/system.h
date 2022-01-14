#ifndef MSIOSYSTEM_H
#define MSIOSYSTEM_H

#include <casacore/casa/OS/HostInfo.h>

#include <stdio.h>
#include <unistd.h>
#include <sched.h>
#include <string.h>

#include <stdexcept>
#include <string>

#include <aocommon/threadpool.h>
class System {
 public:
  static long TotalMemory() { return casacore::HostInfo::memoryTotal() * 1024; }

  static unsigned ProcessorCount() { return aocommon::ThreadPool::NCPUs(); }

  static std::string StrError(int errnum) {
    // Because strerror_r() has different return values on different platforms,
    // two overloads of handle_strerror are used to make this compile and work
    // in either case of int or char*.
    char buffer[1024];
    char* ret = handle_strreturn(strerror_r(errnum, buffer, 1024));
    if (ret == nullptr)
      return std::string(buffer);
    else
      return std::string(ret);
  }

  static std::string FindPythonFilePath(const std::string& filename);

 private:
  static char* handle_strreturn(int value) {
    if (value != 0) throw std::runtime_error("strerror_r() reported an error");
    return nullptr;
  }
  static char* handle_strreturn(char* value) { return value; }
};

#endif  // MSIOSYSTEM_H
