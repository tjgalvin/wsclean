module purge

module load spack/20250523

module load gcc
module load cmake
module load boost
module load casacore
module load cfitsio
module load fftw
module load hdf5
module load openblas
module load openmpi
module load gsl

module load dp3
module load everybeam
module load idg

# Use a venv instead of spack for python packages
source /var/software/spack-extras/20250523/wsclean-ci-env/bin/activate
