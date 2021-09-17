#include "slave.h"

#include "taskmessage.h"

#include "../main/commandline.h"
#include "../main/wsclean.h"

#include "../io/logger.h"
#include <aocommon/checkblas.h>

#include <exception>
#include <iostream>

#include <mpi.h>

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
  bool master = (rank == 0);

  char hostname[256];
  gethostname(hostname, sizeof(hostname));
  std::cout << "Node " << rank << ", PID " << getpid() << " on " << hostname
            << "\n";

  // During parsing of parameters, we don't want all processes to report
  // bad parameters. This variable is used to keep track if full errors
  // should be reported
  bool shortException = false;
  try {
    bool parseResult = false;
    shortException = !master;
    parseResult = CommandLine::ParseWithoutValidation(
        wsclean, argc, const_cast<const char**>(argv), !master);
    shortException = !master && !Logger::IsVerbose();
    check_openblas_multithreading();
    if (parseResult) {
      CommandLine::Validate(wsclean);
      Settings& settings = wsclean.GetSettings();
      settings.useMPI = true;
      shortException = false;
      if (master) {
        CommandLine::Run(wsclean);
        TaskMessage message;
        message.type = TaskMessage::Type::kFinish;
        message.bodySize = 0;
        aocommon::SerialOStream msgStream;
        message.Serialize(msgStream);
        for (int i = 1; i != world_size; ++i) {
          MPI_Send(msgStream.data(), msgStream.size(), MPI_BYTE, i, 0,
                   MPI_COMM_WORLD);
        }
      } else {
        Slave slave(settings);
        slave.Run();
      }
    }
    Logger::Error << "Process " << rank << " finished.\n";
  } catch (std::exception& e) {
    if (shortException)
      Logger::Error << "Process " << rank
                    << " stopped because of an exception.\n";
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
