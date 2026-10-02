#include "mpiworkerscheduler.h"

#include <mpi.h>

#include <aocommon/io/serialostream.h>
#include <aocommon/io/serialistream.h>
#include <aocommon/logger.h>

#include "../distributed/taskmessage.h"
#include "../distributed/mpibig.h"

namespace wsclean {

namespace {
constexpr int kMainNode = 0;
constexpr int kTag = 0;
}  // namespace

MpiWorkerScheduler::MpiWorkerScheduler(const Settings& settings)
    : GriddingTaskManager{settings},
      rank_{-1},
      local_scheduler_{settings.parallelGridding, settings} {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  local_scheduler_.SetWriterLockManager(*this);
}

void MpiWorkerScheduler::Run(
    GriddingTask&& task, [[maybe_unused]] std::function<void(GriddingResult&)> ignored_callback) {
  aocommon::Logger::Info << "Worker node " << rank_ << " is starting gridding task "
                         << task.unique_id << ".\n";
  local_scheduler_.Run(std::move(task), [this](GriddingResult& result) {
    aocommon::Logger::Info << "Worker node " << rank_ << " has finished gridding task "
                           << result.unique_id << ".\n";

    aocommon::SerialOStream resStream;
    resStream.UInt64(0);  // reserve nr of packages for MPI_Send_Big
    result.Serialize(resStream);

    const TaskMessage message(TaskMessage::Type::kGriddingResult, resStream.size());
    aocommon::SerialOStream msgStream;
    message.Serialize(msgStream);
    assert(msgStream.size() == TaskMessage::kSerializedSize);

    std::lock_guard<std::mutex> lock(mutex_);
    MPI_Send(msgStream.data(), msgStream.size(), MPI_BYTE, kMainNode, kTag, MPI_COMM_WORLD);
    MPI_Send_Big(resStream.data(), resStream.size(), kMainNode, kTag, MPI_COMM_WORLD,
                 GetSettings().maxMpiMessageSize);
  });
}

}  // namespace wsclean
