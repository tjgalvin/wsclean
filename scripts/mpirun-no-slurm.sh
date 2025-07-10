#!/bin/bash
# Ensure that OpenMPI does not detect a SLURM environment so it just runs
# without doing anything smart.
unset SLURM_JOBID

# Changing the call by the following commented line can
# be useful for debugging:
# exec mpirun -np 2 xterm -e gdb -ex run --args "${@:3}"

exec mpirun "$@"


