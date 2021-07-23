add_executable(
  unittest EXCLUDE_FROM_ALL
  test.cpp
  testbaselinedependentaveraging.cpp
  testclean.cpp
  testcomponentlist.cpp
  testfitsdateobstime.cpp
  testfluxdensity.cpp
  testimageset.cpp
  testcachedimageset.cpp
  testmultibanddata.cpp
  testparsetreader.cpp
  testpeakfinder.cpp
  testserialization.cpp
  math/tdijkstrasplitter.cpp
  math/tgaussianfitter.cpp
  math/tpolynomialfitter.cpp
  math/tnlplfitter.cpp
  math/tpolynomialchannelfitter.cpp
  msproviders/tnoisemsrowprovider.cpp
  structures/testimage.cpp
  structures/testimagingtable.cpp
  ${WSCLEANFILES})

add_executable(
  integrationtest EXCLUDE_FROM_ALL
  test.cpp integration/tbasics.cpp integration/tstitchfacets.cpp
  integration/tveladeconvolution.cpp ${WSCLEANFILES})

set(TEST_LIBRARIES pybind11::embed ${ALL_LIBRARIES}
                   ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

target_link_libraries(unittest PRIVATE ${TEST_LIBRARIES})
target_link_libraries(integrationtest PRIVATE ${TEST_LIBRARIES})

add_custom_target(
  checkunit
  DEPENDS unittest
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
  DEPENDS integrationtest
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_mwa_ms.sh
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_mwa_coeff.sh ${DOWNLOAD_COEFF}
  COMMAND
    ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/data/ds9_2facets.reg
    ${CMAKE_BINARY_DIR}/test_data/ds9_2facets.reg
  COMMAND integrationtest -l unit_scope)

add_custom_target(check DEPENDS checkunit;checkintegration)

add_custom_target(
  coverage
  DEPENDS check
  COMMAND gcovr -r .. -e '.*/tests/.*' -e '.*/CompilerIdCXX/.*')

# The 'source' symbolic link simplifies running the tests manually inside
# ${CMAKE_BINARY_DIR}: It allows using 'python3 -m pytest [OPTIONS] source/tests.py' instead of
# 'python3 -m pytest [OPTIONS] ../wsclean/tests/systemtests/tests.py'. (Using 'RESULT', fatal
# errors won't occur on systems without symlink support.)
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                        ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/source)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/testconfig.py.in
               ${CMAKE_BINARY_DIR}/testconfig.py)

add_custom_target(
  checkfacets
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMAND
    python3 -m pytest --exitfirst -v
    --junitxml=${CMAKE_BINARY_DIR}/test_facets.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/test_facets.py)

# Test a number of wsclean commands, which can be considered system level tests3
# FIXME: maybe other name for target, e.g. checkwscleancommands?
add_custom_target(
  checkcommandcatalogue
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMAND
    python3 -m pytest --exitfirst -v
    --junitxml=${CMAKE_BINARY_DIR}/test_command_catalogue.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/test_command_catalogue.py)

# Collect all tests that should run on das5 in a single target
add_custom_target(checkintegration-nightly DEPENDS checkintegration;checkfacets)
