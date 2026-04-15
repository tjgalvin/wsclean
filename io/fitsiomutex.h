#ifndef WSCLEAN_IO_FITSIOMUTEX_H_
#define WSCLEAN_IO_FITSIOMUTEX_H_

#include <mutex>

namespace wsclean::io {

/**
 * Process-wide guard for FITS I/O.
 *
 * CachedImageSet and WSCFitsWriter perform FITS reads/writes from worker
 * threads. WSClean already initializes cfitsio once during startup to avoid an
 * initialization race; this mutex keeps the actual FitsReader/FitsWriter calls
 * serialized as well.
 */
inline std::mutex& FitsIOMutex() {
  static std::mutex mutex;
  return mutex;
}

}  // namespace wsclean::io

#endif  // WSCLEAN_IO_FITSIOMUTEX_H_
