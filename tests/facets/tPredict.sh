#!/bin/bash

# Runs a predict on a provided model image with and without faceting
# Script is assumed to be run from within the build/ directory

set -e

# Make sure that wsclean executable from build directory is used
export PATH=${PWD}:$PATH
cd ${PWD}/test_data

MWA_MOCK_ARCHIVE=MWA_ARCHIVE.tar.bz2
MWA_MOCK_MS=MWA_MOCK.ms
# MS used for full image and facet based predict respectively
MWA_MOCK_FULL=MWA_MOCK_FULL.ms
MWA_MOCK_FACET=MWA_MOCK_FACET.ms

# Download a mock measurement set from the citt server
if [ ! -f ${MWA_MOCK_MS}/table.f1 ]; then

    if [ ! -f "$MWA_MOCK_ARCHIVE" ]; then
	wget -q www.astron.nl/citt/EveryBeam/MWA-single-timeslot.tar.bz2 -O $MWA_MOCK_ARCHIVE
    fi

    mkdir -p $MWA_MOCK_MS

    tar -xf $MWA_MOCK_ARCHIVE  -C $MWA_MOCK_MS --strip-components=1
    rm $MWA_MOCK_ARCHIVE
fi

cp -r $MWA_MOCK_MS $MWA_MOCK_FULL
cp -r $MWA_MOCK_MS $MWA_MOCK_FACET

# Predict settings
common="-predict -use-wgridder -name point-source"
facetfile="ds9_facet.reg"

echo "===== Predicting full image ====="
wsclean ${common} ${MWA_MOCK_FULL}

echo "===== Predicting facet ====="
wsclean ${common} -facet-regions ${facetfile} ${MWA_MOCK_FACET}

# taql requires casacore(-tools)
taql "select from MWA_MOCK_FULL.ms t1, MWA_MOCK_FACET.ms t2 where not all(near(t1.MODEL_DATA,t2.MODEL_DATA,3e-3))" > taql.out

# Create expected taql output.
echo "    select result of 0 rows" > taql.ref
diff taql.out taql.ref  ||  exit 1