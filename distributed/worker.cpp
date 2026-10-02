#include "worker.h"

#include <mpi.h>

#include <aocommon/io/serialistream.h>
#include <aocommon/io/serialostream.h>
#include <aocommon/logger.h>

#include "../scheduling/griddingtask.h"

#include "mpibig.h"
#include "taskmessage.h"

namespace wsclean {

namespace {
constexpr int kMainNode = 0;
constexpr int kTag = 0;
}  // namespace

void Worker::Run() {
  TaskMessage message;
  const size_t maximum_message_size = scheduler_.GetSettings().maxMpiMessageSize;

  do {
    MPI_Status status;

    aocommon::UVector<unsigned char> buffer(TaskMessage::kSerializedSize);
    MPI_Recv(buffer.data(), TaskMessage::kSerializedSize, MPI_BYTE, kMainNode, kTag, MPI_COMM_WORLD,
             &status);
    aocommon::SerialIStream stream(std::move(buffer));
    message.Unserialize(stream);

    switch (message.type) {
      case TaskMessage::Type::kStart:
        scheduler_.Start(message.n_writer_groups, message.n_parallel_tasks);
        break;

      case TaskMessage::Type::kFinish:
        // The do..while loop will exit, too.
        break;

      case TaskMessage::Type::kGriddingRequest: {
        buffer.resize(message.body_size);
        MPI_Recv_Big(buffer.data(), message.body_size, kMainNode, kTag, MPI_COMM_WORLD, &status,
                     maximum_message_size);
        aocommon::SerialIStream stream(std::move(buffer));
        stream.UInt64();  // skip the nr of packages

        GriddingTask task;
        task.Unserialize(stream);
        scheduler_.Run(std::move(task), [](GriddingResult&) {});
        break;
      }

      default:
        aocommon::Logger::Warn << "wsclean-mp MPI worker received an unknown message type!\n";
        break;
    }

  } while (message.type != TaskMessage::Type::kFinish);
  aocommon::Logger::Info << "Worker node " << scheduler_.Rank() << " received exit message.\n";
}

}  // namespace wsclean
