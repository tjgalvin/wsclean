#include "wsclean.h"

#include <boost/program_options/parsers.hpp>

#include "../main/commandline.h"
#include "../main/wsclean.h"

namespace wsclean {

void Image(const std::string& command_line_parameters, InMemoryMs&& ms) {
  const std::vector<std::string> parameters =
      boost::program_options::split_unix(command_line_parameters);
  std::vector<const char*> argv;
  argv.reserve(parameters.size());
  for (const std::string& p : parameters) argv.emplace_back(p.c_str());

  WSClean wsclean;
  CommandLine::Parse(wsclean, parameters.size(), argv.data(), false);

  // TODO in memory data must be reordered (reordering not yet implemented).
}

}  // namespace wsclean
