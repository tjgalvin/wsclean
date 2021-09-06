# Extend MODULEPATH, so that openmpi modulefiles can be found
if [[ ${MODULEPATH} != *"/cm/shared/modulefiles"* ]]; then
  export MODULEPATH=$MODULEPATH:/cm/shared/modulefiles
fi

module load git/2.30.0
module load cmake/3.16.2
module load openmpi/gcc/64/4.0.2

# Everybeam loads casacore, which loads many more modules.
module load everybeam/LATEST-gcc-8.3.0
module load idg/LATEST-gcc-8.3.0-CUDA