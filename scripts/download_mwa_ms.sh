#!/bin/sh
#
# Author: Jakob Maljaars
# Email: jakob.maljaars_@_stcorp.nl

# Script for downloading a mock MWA Measurement set

set -e

# Download into a 'test_data' directory inside the current directory.
mkdir -p test_data
cd test_data/

MWA_MOCK_ARCHIVE=MWA_ARCHIVE.tar.bz2
MWA_MOCK_MS=MWA_MOCK.ms

if [ ! -f ${MWA_MOCK_MS}/table.f1 ]; then

    if [ ! -f "$MWA_MOCK_ARCHIVE" ]; then
	wget -q www.astron.nl/citt/EveryBeam/MWA-single-timeslot.tar.bz2 -O $MWA_MOCK_ARCHIVE
    fi

    mkdir -p $MWA_MOCK_MS

    tar -xf $MWA_MOCK_ARCHIVE  -C $MWA_MOCK_MS --strip-components=1
    rm $MWA_MOCK_ARCHIVE
fi
