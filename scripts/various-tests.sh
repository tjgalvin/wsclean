#!/bin/sh

if [[ "$1" == "" ]] ; then
    echo Syntax: various-tests.sh \<ms\> [optional] \<mwa_coefficients\>
    exit 1
fi

export MWA_MS="$1"
export MWA_COEFFS_PATH="$2"
# Let pytest know not to remove MWA_MS and MWA_COEFFS_PATH
# after the last test is finished
export SKIP_CLEANUP="true"

dims="-size 1024 1024 -scale 1amin"
rectdims="-size 1536 1024 -scale 1amin"
facetfile_2facets="../tests/data/ds9_2facets.reg"
facetfile_nfacets="../tests/data/ds9_nfacets.reg"

# workdir should point to the directory where testconfig.py is located
# this typically is the build directory
workdir="../../build"

cd $workdir

# Execute tests
python3 -B -m pytest -v source/tests.py
