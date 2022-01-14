add_executable(
  unittest EXCLUDE_FROM_ALL
  test.cpp
  testbaselinedependentaveraging.cpp
  testclean.cpp
  testcomponentlist.cpp
  testcommandline.cpp
  testfitsdateobstime.cpp
  testimageset.cpp
  testcachedimageset.cpp
  testparsetreader.cpp
  testpeakfinder.cpp
  testserialization.cpp
  math/tdijkstrasplitter.cpp
  math/tgaussianfitter.cpp
  math/tpolynomialfitter.cpp
  math/tnlplfitter.cpp
  math/tpolynomialchannelfitter.cpp
  msproviders/tnoisemsrowprovider.cpp
  msproviders/tbdamsrowproviderdata.cpp
  msproviders/tbdamsrowprovider.cpp
  msproviders/tmsrowproviderbase.cpp
  structures/testimage.cpp
  structures/testimagingtable.cpp
  ${WSCLEANFILES})

add_custom_target(
  download_mocks
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_mwa_ms.sh
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
  checkunit
  DEPENDS unittest download_mocks mwa_bda_mock_ms
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMAND unittest)

if(EVERYBEAM_LIB AND IDGAPI_LIBRARIES)
  set(DOWNLOAD_COEFF true)
else()
  set(DOWNLOAD_COEFF false)
endif()

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/test_data/)

add_custom_target(
  checkintegration
  COMMAND
    python3 -B -m pytest --exitfirst -v
    --junitxml=${CMAKE_BINARY_DIR}/test_integration.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/test_basics.py
    ${CMAKE_CURRENT_SOURCE_DIR}/test_veladeconvolution.py)

add_custom_target(check DEPENDS checkunit;checkintegration)

add_custom_target(
  coverage
  DEPENDS check
  COMMAND gcovr -r .. -e '.*/tests/.*' -e '.*/CompilerIdCXX/.*')

# The 'source' symbolic link simplifies running the tests manually inside
# ${CMAKE_CURRENT_BINARY_DIR}:
# It allows using 'python3 -m pytest [OPTIONS] source/tests.py' instead of
# 'python3 -m pytest [OPTIONS] ../wsclean/tests/systemtests/tests.py'.
execute_process(
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_BINARY_DIR}/source)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/testconfig.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/testconfig.py)

add_custom_target(
  checkfacets
  COMMAND
    python3 -B -m pytest --exitfirst -v
    --junitxml=${CMAKE_BINARY_DIR}/test_facets.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/test_facets.py)

# Test a number of wsclean commands, which can be considered system level tests3
# FIXME: maybe other name for target, e.g. checkwscleancommands?
add_custom_target(
  checkcommandcatalogue
  COMMAND
    python3 -B -m pytest --exitfirst -v
    --junitxml=${CMAKE_BINARY_DIR}/test_command_catalogue.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/test_command_catalogue.py)

# Collect all tests that should run on das5 in a single target
add_custom_target(checkintegration-nightly DEPENDS checkintegration;checkfacets)
