module purge

module load spack/20251124

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
module load py-pybind11

module load dp3
module load everybeam
module load idg

# Use a venv instead of spack for python packages
source /var/software/spack-extras/20250523/wsclean-ci-env/bin/activate
export EVERYBEAM_DATADIR=/var/software/spack-extras/20250109/opt/spack/linux-rocky8-zen2/gcc-13.1.0/everybeam-0.6.1-wm3lduhu6psslnrqr6jyc3ru5l4vvnou/share/everybeam/
