#include "wsclean.h"

#include <boost/program_options/parsers.hpp>

#include "../main/commandline.h"
#include "../main/wsclean.h"

namespace wsclean {

void Image(const std::string& command_line_parameters, std::vector<InMemoryMs>&& ms_list) {
  const std::vector<std::string> parameters =
      boost::program_options::split_unix(command_line_parameters);
  std::vector<const char*> argv;
  argv.reserve(parameters.size() + 1);
  // The first parameter is the executable, but is unused:
  argv.emplace_back();
  for (const std::string& p : parameters) argv.emplace_back(p.c_str());

  WSClean wsclean;
  constexpr bool kSlave = false;
  constexpr bool kInMemoryData = true;
  CommandLine::Parse(wsclean, argv.size(), argv.data(), kSlave, kInMemoryData);

  std::vector<std::unique_ptr<InMemoryMs>> in_memory_ms_list;
  in_memory_ms_list.reserve(ms_list.size());
  for (InMemoryMs& ms : ms_list) {
    in_memory_ms_list.emplace_back(std::make_unique<InMemoryMs>(std::move(ms)));
  }
  ms_list.clear();

  wsclean.SetInMemoryData(std::move(in_memory_ms_list));
  wsclean.RunClean();
}

}  // namespace wsclean
