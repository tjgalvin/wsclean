add_executable(
  unittest
  test.cpp
  tbaselinedependentaveraging.cpp
  tfitsdateobstime.cpp
  tserialization.cpp
  gridding/taveragecorrection.cpp
  gridding/tvisibilitymodifier.cpp
  idg/taveragebeam.cpp
  io/tcachedimageset.cpp
  io/tparsetreader.cpp
  main/tcommandline.cpp
  main/tmshelper.cpp
  main/tsettings.cpp
  math/timageoperations.cpp
  math/trenderer.cpp
  math/ttophatconvolution.cpp
  model/tbbsmodel.cpp
  msproviders/tnoisemsrowprovider.cpp
  msproviders/tbdamsrowprovider.cpp
  msproviders/tbdamsrowproviderdata.cpp
  msproviders/tmsrowproviderbase.cpp
  scheduling/tgriddingtaskfactory.cpp
  structures/timagingtable.cpp
  structures/tfacetutil.cpp
  structures/toutputchannelinfo.cpp
  structures/tprimarybeam.cpp
  structures/tprimarybeamimageset.cpp
  structures/tresources.cpp
  system/tmappedfile.cpp)

target_compile_definitions(
  unittest
  PRIVATE
    FACET_DEFINITION_FILE_2FACETS="${CMAKE_CURRENT_SOURCE_DIR}/data/ds9_2facets.reg"
    FACET_DEFINITION_FILE_4FACETS="${CMAKE_CURRENT_SOURCE_DIR}/data/ds9_4facets.reg"
)

target_link_libraries(
  unittest
  PRIVATE wsclean-lib ${CASACORE_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY}
          ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

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

ExternalProject_Add(
  gmrt_ms
  URL https://support.astron.nl/software/ci_data/EveryBeam/GMRT_MOCK.tar.bz2
  URL_HASH
    SHA256=5286b6f19f226ed97a8198438edcc3c5b09481376ac134135dd595107dc39fc3
  SOURCE_DIR ${CMAKE_BINARY_DIR}/test_data/GMRT_MOCK.ms
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND "")

ExternalProject_Add(
  lofar_3c196_ms
  URL ${WSCLEAN_DATA_URL}/LOFAR_3C196.ms.tar.bz2
  URL_HASH
    SHA256=f834f4eefb2ceffc09194791744c53ebc5154e947a5437a8ea42d259ccdd5fc3
  SOURCE_DIR ${CMAKE_BINARY_DIR}/test_data/LOFAR_3C196.ms
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND "")

add_custom_target(
  unitchecks
  DEPENDS unittest download_mocks jvla_test_ms mwa_bda_mock_ms ska_low_ms
          gmrt_ms
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMAND unittest)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/test_data/)

add_subdirectory(python)

add_custom_target(check DEPENDS unitchecks;integrationchecks)

add_custom_target(
  coverage
  DEPENDS check
  COMMAND gcovr -r ${CMAKE_SOURCE_DIR} -e '.*/tests/.*' -e '.*/external/.*' -e
          '.*/CompilerIdCXX/.*' ${CMAKE_BINARY_DIR})

# Collect all tests that should run on das5 in a single target
add_custom_target(nightlychecks DEPENDS integrationchecks;facetchecks)
