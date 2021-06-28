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
  COMMAND integrationtest)

add_custom_target(check DEPENDS checkunit;checkintegration)

add_custom_target(
  coverage
  DEPENDS check
  COMMAND gcovr -r .. -e '.*/tests/.*' -e '.*/CompilerIdCXX/.*')

# Test the facet based and non-facet based predict
add_custom_target(
  checkpredict
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/tests/facets/resources/ds9_facet.reg
    ${CMAKE_BINARY_DIR}/test_data/ds9_facet.reg
  COMMAND
    ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/tests/facets/resources/point-source-model.fits
    ${CMAKE_BINARY_DIR}/test_data/point-source-model.fits
  COMMAND ${CMAKE_SOURCE_DIR}/tests/facets/tPredict.sh)

# Test the facet based imaging
add_custom_target(
  checkfacetdegridding
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/tests/facets/resources/ds9_facet.reg
    ${CMAKE_BINARY_DIR}/test_data/ds9_facet.reg
  COMMAND
    ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/tests/facets/resources/point-source-model.fits
    ${CMAKE_BINARY_DIR}/test_data/point-source-model.fits
  COMMAND ${CMAKE_SOURCE_DIR}/tests/facets/tFacetDegridding.sh)

# Collect all tests that should run on das5 in a single target
add_custom_target(checkintegration-nightly
                  DEPENDS checkintegration;checkpredict;checkfacetdegridding)

# The 'source' symbolic link simplifies running the tests manually inside
# ${CMAKE_BINARY_DIR}: It allows using 'python3 -m pytest [OPTIONS] source/tests.py' instead of
# 'python3 -m pytest [OPTIONS] ../wsclean/tests/systemtests/tests.py'. (Using 'RESULT', fatal
# errors won't occur on systems without symlink support.)
execute_process(
  COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/systemtests ${CMAKE_BINARY_DIR}/source)

# Test a number of wsclean commands, which can be considered system level tests3
# FIXME: maybe other name for target, e.g. checkwscleancommands?
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/systemtests/testconfig.py.in
               ${CMAKE_BINARY_DIR}/testconfig.py)
add_custom_target(
  checksystemtests
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMAND
    python3 -m pytest --exitfirst -v
    --junitxml=${CMAKE_BINARY_DIR}/system_tests.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/systemtests/tests.py)
