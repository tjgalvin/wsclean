#include "mpibig.h"

#include <aocommon/logger.h>

#include <algorithm>
#include <cstdint>
#include <limits>

using aocommon::Logger;

int MPI_Send_Big(unsigned char* buf, size_t count, int dest, int tag,
                 MPI_Comm comm) {
  size_t nPackages = (count + std::numeric_limits<int>::max() - 1) /
                     std::numeric_limits<int>::max();

  *reinterpret_cast<uint64_t*>(buf) = nPackages;

  Logger::Debug << "Sending " << nPackages << " packages...\n";
  for (size_t i = 0; i != nPackages - 1; ++i) {
    const unsigned char* partBuffer = buf + i * std::numeric_limits<int>::max();
    int returnValue = MPI_Send(partBuffer, std::numeric_limits<int>::max(),
                               MPI_BYTE, dest, tag, comm);
    if (returnValue != MPI_SUCCESS) return returnValue;
    Logger::Debug << "Package " << (i + 1) << " sent.\n";
  }

  const unsigned char* partBuffer =
      buf + (nPackages - 1) * std::numeric_limits<int>::max();
  size_t partCount = count % std::numeric_limits<int>::max();
  int returnValue = MPI_Send(partBuffer, partCount, MPI_BYTE, dest, tag, comm);
  Logger::Debug << "Package " << nPackages << " sent.\n";
  return returnValue;
}

int MPI_Recv_Big(unsigned char* buf, size_t count, int source, int tag,
                 MPI_Comm comm, MPI_Status* status) {
  int firstSize = std::min<size_t>(std::numeric_limits<int>::max(), count);
  int returnValue =
      MPI_Recv(buf, firstSize, MPI_BYTE, source, tag, comm, status);
  if (returnValue != MPI_SUCCESS) return returnValue;

  size_t nPackages = *reinterpret_cast<uint64_t*>(buf);
  buf += firstSize;
  count -= size_t(firstSize);

  Logger::Debug << "Received package 1/" << nPackages << ".\n";
  for (size_t i = 1; i != nPackages; ++i) {
    int partSize = std::min<size_t>(std::numeric_limits<int>::max(), count);
    returnValue = MPI_Recv(buf, partSize, MPI_BYTE, source, tag, comm, status);
    if (returnValue != MPI_SUCCESS) return returnValue;

    buf += partSize;
    count -= size_t(partSize);
    Logger::Debug << "Received package " << (i + 1) << "/" << nPackages
                  << ".\n";
  }
  return MPI_SUCCESS;
}
