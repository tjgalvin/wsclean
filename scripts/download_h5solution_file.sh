#!/bin/sh

# Script for downloading a H5 solution file

set -e

mkdir -p test_data
cd test_data/

H5_SOLUTION_FILE=mock_soltab_2pol.h5

if [ ! -f ${H5_SOLUTION_FILE} ] ; then
    wget -N -q www.astron.nl/citt/ci_data/wsclean/${H5_SOLUTION_FILE}
fi
