#!/bin/sh
# Ensure that OpenMPI does not detect a SLURM environment so it just runs
# without doing anything smart.
unset SLURM_JOBID
exec mpirun "$@"