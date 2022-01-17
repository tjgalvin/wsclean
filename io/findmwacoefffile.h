#ifndef WSCLEAN_MWA_FINDCOEFFFILE_H_
#define WSCLEAN_MWA_FINDCOEFFFILE_H_

#include "../system/pythonfilepath.h"

#include <boost/filesystem.hpp>

#include <string>

#ifndef DEFAULT_H5_FILE
#define DEFAULT_H5_FILE "mwa_full_embedded_element_pattern.h5"
#endif

#ifndef DEFAULT_H5_FILE_PATH
#define DEFAULT_H5_FILE_PATH "mwapy/data/"
#endif

namespace wsclean {
namespace mwa {

/**
 * @brief Search for MWA h5 coefficient file
 * (mwa_full_embedded_element_pattern.h5) on the suggested path \param
 * search_path.
 *
 * @param search_path Search path
 * @return std::string Path to file
 */
inline static std::string FindCoeffFile(const std::string& search_path) {
  std::string h5_path;
  if (search_path.empty()) {
    std::string h5_test_path = DEFAULT_H5_FILE_PATH;
    h5_test_path += DEFAULT_H5_FILE;
    h5_path = wsclean::system::FindPythonFilePath(h5_test_path);
  } else {
    boost::filesystem::path p =
        boost::filesystem::path(search_path) / DEFAULT_H5_FILE;
    if (!boost::filesystem::exists(p))
      throw std::runtime_error(
          "Manually specified MWA directory did not contain H5 beam file: '" +
          p.string() + "' not found.");
    h5_path = p.string();
  }
  return h5_path;
}
}  // namespace mwa
}  // namespace wsclean

#endif  // WSCLEAN_MWA_FINDCOEFFFILE_H_
