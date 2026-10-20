#include "worker.h"

#include "taskmessage.h"

#include "../main/commandline.h"
#include "../main/wsclean.h"

#include <aocommon/checkblas.h>
#include <aocommon/logger.h>

#include <exception>
#include <iostream>

#include <mpi.h>

using aocommon::Logger;

using wsclean::CommandLine;
using wsclean::Settings;
using wsclean::TaskMessage;
using wsclean::Worker;
using wsclean::WSClean;

namespace {
void SetMpiSettings(Settings& settings, size_t n_nodes) {
  settings.nMpiNodes = n_nodes;

  if (settings.channelToNode.empty()) {
    // Create a default channel-to-node mapping.
    settings.channelToNode.reserve(settings.channelsOut);
    if (!settings.masterDoesWork) --n_nodes;
    for (size_t channelIndex = 0; channelIndex < settings.channelsOut; ++channelIndex) {
      // Use a round-robin distribution, since later channels have higher
      // frequencies and gridding is therefore more expensive.
      size_t node_index = channelIndex % n_nodes;
      if (!settings.masterDoesWork) ++node_index;
      settings.channelToNode.push_back(node_index);
    }
  }
}
}  // namespace

int main(int argc, char* argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  if (provided != MPI_THREAD_MULTIPLE) {
    std::cout << "This MPI implementation does not support multiple threads.\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int result = 0;
  WSClean wsclean;
  int world_size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  bool main = (rank == 0);

  char hostname[256];
  gethostname(hostname, sizeof(hostname));
  std::cout << "Node " << rank << ", PID " << getpid() << " on " << hostname << "\n";

  // During parsing of parameters, we don't want all processes to report
  // bad parameters. This variable is used to keep track if full errors
  // should be reported
  bool shortException = false;
  try {
    bool parseResult = false;
    shortException = !main;
    constexpr bool kInMemory = false;
    parseResult = CommandLine::ParseWithoutValidation(wsclean, argc, const_cast<const char**>(argv),
                                                      !main, kInMemory);
    shortException = !main && !Logger::IsVerbose();
    check_openblas_multithreading();
    if (parseResult) {
      Settings& settings = wsclean.GetSettings();
      SetMpiSettings(settings, world_size);
      CommandLine::Validate(wsclean);
      shortException = false;
      if (main) {
        CommandLine::Run(wsclean);
        const TaskMessage message(TaskMessage::Type::kFinish, 0);
        aocommon::SerialOStream msgStream;
        message.Serialize(msgStream);
        for (int i = 1; i != world_size; ++i) {
          MPI_Send(msgStream.data(), msgStream.size(), MPI_BYTE, i, 0, MPI_COMM_WORLD);
        }
      } else {
        if (settings.mode == Settings::PredictMode) {
          wsclean.InitializePredictImageSettings();
        }
        Worker worker(settings);
        worker.Run();
      }
    }
    Logger::Error << "Process " << rank << " finished.\n";
  } catch (std::exception& e) {
    if (shortException)
      Logger::Error << "Process " << rank << " stopped because of an exception.\n";
    else {
      Logger::Error << "+ + + + + + + + + + + + + + + + + + +\n"
                    << "+ An exception occured in process " << rank << ":\n";
      std::istringstream iss(e.what());
      for (std::string line; std::getline(iss, line);) {
        Logger::Error << "+ >>> " << line << "\n";
      }
      Logger::Error << "+ + + + + + + + + + + + + + + + + + +\n";
    }
    result = -1;
  }
  MPI_Finalize();
  return result;
}
