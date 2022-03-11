add_executable(
  unittest
  test.cpp
  testbaselinedependentaveraging.cpp
  testclean.cpp
  testcomponentlist.cpp
  testcommandline.cpp
  testfitsdateobstime.cpp
  testcachedimageset.cpp
  testparsetreader.cpp
  testpeakfinder.cpp
  testprimarybeamimageset.cpp
  testserialization.cpp
  deconvolution/testdeconvolutiontable.cpp
  deconvolution/testimageset.cpp
  idg/taveragebeam.cpp
  math/tdijkstrasplitter.cpp
  math/tgaussianfitter.cpp
  math/tpolynomialfitter.cpp
  math/tnlplfitter.cpp
  math/tpolynomialchannelfitter.cpp
  msproviders/tnoisemsrowprovider.cpp
  msproviders/tbdamsrowproviderdata.cpp
  msproviders/tbdamsrowprovider.cpp
  msproviders/tmsrowproviderbase.cpp
  structures/testimagingtable.cpp
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

# A dataset generated using the following command:
# DP3 msin=MWA-1052736496-averaged.ms msin.ntimes=4 msin.nchan=1 \
#   msin.baseline="1,2,126,217" msout=MWA_BDA_MOCK.ms msout.overwrite=True \
#   steps=[filter,bdaaverager] filter.blrange=[1000,10000] \
#   bdaaverager.timebase=5000
ExternalProject_Add(
  mwa_bda_mock_ms
  URL http://astron.nl/citt/ci_data/wsclean/MWA_BDA_MOCK.ms.tar.bz2
  URL_HASH
    SHA256=cc249aa62771d14d163e2a5d90fbb6110c57c13bba07b762acad875a64a15f41
  SOURCE_DIR ${CMAKE_BINARY_DIR}/test_data/MWA_BDA_MOCK.ms
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND "")

set(TEST_LIBRARIES pybind11::embed ${ALL_LIBRARIES}
                   ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

target_link_libraries(unittest PRIVATE ${TEST_LIBRARIES})

add_custom_target(
  unitchecks
  DEPENDS unittest download_mocks mwa_bda_mock_ms
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
