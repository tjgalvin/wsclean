#ifndef MPI_BIG_H
#define MPI_BIG_H

#include <mpi.h>

/**
 * Send a big message. If count does not fit in an integer, the message is split
 * into multiple messages. The first 8 bytes of the buf are changed and used to
 * communicate the number of messages. Therefore the first 8 bytes HAVE to be
 * reserved. The parameters are equivallent to MPI_Send, except it will always
 * use the type MPI_Byte.
 * @returns MPI_SUCCESS when successful.
 */
int MPI_Send_Big(unsigned char* buf, size_t count, int dest, int tag,
                 MPI_Comm comm);

/** Receives a big message.
 * The first 8 bytes of the buffer will receive the number of
 * packages. The count parameter should match the one given to @ref
 * MPI_Send_Big().
 * @returns MPI_SUCCESS when successful.
 */
int MPI_Recv_Big(unsigned char* buf, size_t count, int source, int tag,
                 MPI_Comm comm, MPI_Status* status);

#endif
