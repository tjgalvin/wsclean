#include "pythonfilepath.h"

#include "../io/logger.h"

#include <boost/filesystem.hpp>

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <random>

namespace wsclean {
namespace system {
std::string FindPythonFilePath(const std::string& filename) {
  if (boost::filesystem::exists(filename)) return filename;
  Logger::Debug << "Searching " << filename << "... ";
  Logger::Debug.Flush();
  std::random_device rndDev;
  std::mt19937 gen(rndDev());
  std::stringstream filenameStr;
  filenameStr << "/tmp/ao-python-path-list" << gen() << ".tmp";
  boost::filesystem::path tempPath =
      filenameStr.str();  // boost::filesystem::unique_path();
  const std::string tempFilename = tempPath.string();  // optional
  std::string command =
      std::string(
          "echo \"from __future__ import print_function\nimport sys\nfor a in "
          "sys.path:\n  print(a)\"|python>") +
      tempFilename;
  int status = std::system(command.c_str());
  if (status != 0) {
    // retry with python3
    command = std::string(
                  "echo \"from __future__ import print_function\nimport "
                  "sys\nfor a in "
                  "sys.path:\n  print(a)\"|python3>") +
              tempFilename;
    status = std::system(command.c_str());
  }
  if (status != 0)
    throw std::runtime_error(
        "std::system() returned non-zero error code: might be out of memory, "
        "or "
        "python might not be working properly.\nCommand was:\n" +
        command);
  std::ifstream searchPathsFile(tempFilename.c_str());
  if (!searchPathsFile.good())
    throw std::runtime_error(("Error in findPythonFilePath: system call did "
                              "not create expected temporary file " +
                              tempFilename)
                                 .c_str());
  while (searchPathsFile.good()) {
    std::string prefixPath;
    std::getline(searchPathsFile, prefixPath);
    boost::filesystem::path searchPath(prefixPath);
    searchPath /= filename;

    bool pathExists = false;
    try {
      pathExists = boost::filesystem::exists(searchPath);
    } catch (...) {
    }
    if (pathExists) {
      const std::string result = searchPath.string();
      Logger::Debug << result << '\n';
      searchPathsFile.close();
      boost::filesystem::remove(tempPath);
      return result;
    }
  }
  searchPathsFile.clear();
  searchPathsFile.seekg(0, std::ios::beg);
  std::string err(std::string("Could not find Python file ") + filename +
                  ". Paths searched:\n");
  while (searchPathsFile.good()) {
    std::string prefixPath;
    std::getline(searchPathsFile, prefixPath);
    err += prefixPath + '\n';
  }
  searchPathsFile.close();
  boost::filesystem::remove(tempPath);
  throw std::runtime_error(err);
}
}  // namespace system
}  // namespace wsclean