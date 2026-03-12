#include "mpibig.h"

#include <aocommon/logger.h>

#include <algorithm>
#include <cstdint>

using aocommon::Logger;

namespace wsclean {

int MPI_Send_Big(unsigned char* buf, size_t count, int dest, int tag, MPI_Comm comm,
                 size_t maximum_message_size) {
  size_t n_packages = (count + maximum_message_size - 1) / maximum_message_size;

  *reinterpret_cast<uint64_t*>(buf) = n_packages;

  Logger::Debug << "Sending " << n_packages << " packages...\n";
  for (size_t i = 0; i != n_packages - 1; ++i) {
    const unsigned char* part_buffer = buf + i * maximum_message_size;
    int return_value = MPI_Send(part_buffer, maximum_message_size, MPI_BYTE, dest, tag, comm);
    if (return_value != MPI_SUCCESS) return return_value;
    Logger::Debug << "Package " << (i + 1) << " sent.\n";
  }

  const unsigned char* part_buffer = buf + (n_packages - 1) * maximum_message_size;
  size_t part_count = count % maximum_message_size;
  int return_value = MPI_Send(part_buffer, part_count, MPI_BYTE, dest, tag, comm);
  Logger::Debug << "Package " << n_packages << " sent.\n";
  return return_value;
}

int MPI_Recv_Big(unsigned char* buf, size_t count, int source, int tag, MPI_Comm comm,
                 MPI_Status* status, size_t maximum_message_size) {
  int first_size = std::min(maximum_message_size, count);
  int return_value = MPI_Recv(buf, first_size, MPI_BYTE, source, tag, comm, status);
  if (return_value != MPI_SUCCESS) return return_value;

  size_t n_packages = *reinterpret_cast<uint64_t*>(buf);
  buf += first_size;
  count -= size_t(first_size);

  Logger::Debug << "Received package 1/" << n_packages << ".\n";
  for (size_t i = 1; i != n_packages; ++i) {
    int part_size = std::min(maximum_message_size, count);
    return_value = MPI_Recv(buf, part_size, MPI_BYTE, source, tag, comm, status);
    if (return_value != MPI_SUCCESS) return return_value;

    buf += part_size;
    count -= size_t(part_size);
    Logger::Debug << "Received package " << (i + 1) << "/" << n_packages << ".\n";
  }
  return MPI_SUCCESS;
}

}  // namespace wsclean
