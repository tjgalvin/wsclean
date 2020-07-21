#include "slave.h"
#include "taskmessage.h"

#include "../scheduling/griddingtask.h"
#include "../scheduling/griddingtaskmanager.h"

#include "../wsclean/logger.h"

#include <mpi.h>

void Slave::Run() {
  TaskMessage message;
  do {
    MPI_Status status;
    MPI_Recv(&message, sizeof(TaskMessage), MPI_BYTE, 0, 0, MPI_COMM_WORLD,
             &status);

    switch (message.type) {
      case TaskMessage::GriddingRequest:
        grid(message.bodySize);
        break;
      default:
        break;
    }

  } while (message.type != TaskMessage::Finish);
  Logger::Info << "Worker node received exit message.\n";
}

void Slave::grid(size_t bodySize) {
  MPI_Status status;
  aocommon::UVector<unsigned char> buffer(bodySize);
  MPI_Recv(buffer.data(), bodySize, MPI_BYTE, 0, 0, MPI_COMM_WORLD, &status);
  SerialIStream stream(std::move(buffer));

  GriddingTask task;
  task.Unserialize(stream);
  std::unique_ptr<GriddingTaskManager> scheduler =
      GriddingTaskManager::Make(_settings, true);
  Logger::Info << "Worker node is starting gridding.\n";
  GriddingResult result = scheduler->RunDirect(task);
  Logger::Info << "Worker node is done gridding.\n";

  SerialOStream resStream;
  result.Serialize(resStream);

  TaskMessage message;
  message.type = TaskMessage::GriddingResult;
  message.bodySize = resStream.size();
  MPI_Send(&message, sizeof(TaskMessage), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
  MPI_Send(resStream.data(), resStream.size(), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
}
