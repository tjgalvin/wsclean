add_executable(
  unittest
  test.cpp
  testbaselinedependentaveraging.cpp
  testcommandline.cpp
  testfitsdateobstime.cpp
  testcachedimageset.cpp
  testparsetreader.cpp
  testprimarybeamimageset.cpp
  testserialization.cpp
  idg/taveragebeam.cpp
  io/tcachedimageset.cpp
  math/trenderer.cpp
  msproviders/tnoisemsrowprovider.cpp
  msproviders/tbdamsrowproviderdata.cpp
  msproviders/tbdamsrowprovider.cpp
  msproviders/tmsrowproviderbase.cpp
  structures/testimagingtable.cpp
  structures/tfacetutil.cpp
  structures/toutputchannelinfo.cpp
  structures/tprimarybeamimageset.cpp
  structures/tresources.cpp
  system/tmappedfile.cpp
  ${WSCLEANFILES})

add_definitions(
  -DFACET_DEFINITION_FILE_2FACETS="${CMAKE_CURRENT_SOURCE_DIR}/data/ds9_2facets.reg"
  -DFACET_DEFINITION_FILE_4FACETS="${CMAKE_CURRENT_SOURCE_DIR}/data/ds9_4facets.reg"
)

add_custom_target(
  download_mocks
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_mwa_ms.sh
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_h5solution_file.sh
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

set(WSCLEAN_DATA_URL https://support.astron.nl/software/ci_data/wsclean)

ExternalProject_Add(
  jvla_test_ms
  URL ${WSCLEAN_DATA_URL}/JVLA-MultiBand-S1_C5-minimal.ms.tar.bz2
  URL_HASH
    SHA256=7c8d41b5ff59c8736b1223e6b855a96e410f27fb4be05d179c8292fdb78cdc7e
  SOURCE_DIR ${CMAKE_BINARY_DIR}/test_data/JVLA-MultiBand-S1_C5-minimal.ms
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND "")

# A dataset generated using the following command:
# DP3 msin=MWA-1052736496-averaged.ms msin.ntimes=4 msin.nchan=1 \
#   msin.baseline="1,2,126,217" msout=MWA_BDA_MOCK.ms msout.overwrite=True \
#   steps=[filter,bdaaverager] filter.blrange=[1000,10000] \
#   bdaaverager.timebase=5000
ExternalProject_Add(
  mwa_bda_mock_ms
  URL ${WSCLEAN_DATA_URL}/MWA_BDA_MOCK.ms.tar.bz2
  URL_HASH
    SHA256=cc249aa62771d14d163e2a5d90fbb6110c57c13bba07b762acad875a64a15f41
  SOURCE_DIR ${CMAKE_BINARY_DIR}/test_data/MWA_BDA_MOCK.ms
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND "")

# This is a 9h long SKA-mid simulated measurement set.
# It contains a point source off the phase center.
# These conditions ensure a strong variablity of the PSF in different directions.
# The dataset is heavily averaged to reduce its size.
ExternalProject_Add(
  ska_low_ms
  URL ${WSCLEAN_DATA_URL}/ska-low-sim-averaged.tar.gz
  URL_HASH
    SHA256=79f723f2165456d715e3019893f1a6701ef92cf9e106c4ea251c95be567215c5
  SOURCE_DIR ${CMAKE_BINARY_DIR}/test_data/ska-low-sim-averaged.ms
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND "")

set(TEST_LIBRARIES pybind11::embed ${ALL_LIBRARIES}
                   ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

target_link_libraries(unittest PRIVATE ${TEST_LIBRARIES})

add_custom_target(
  unitchecks
  DEPENDS unittest download_mocks jvla_test_ms mwa_bda_mock_ms ska_low_ms
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMAND unittest)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/test_data/)

add_subdirectory(python)

add_custom_target(check DEPENDS unitchecks;integrationchecks)

add_custom_target(
  coverage
  DEPENDS check
  COMMAND gcovr -r .. -e '.*/tests/.*' -e '.*/CompilerIdCXX/.*')

# Collect all tests that should run on das5 in a single target
add_custom_target(nightlychecks DEPENDS integrationchecks;facetchecks)
