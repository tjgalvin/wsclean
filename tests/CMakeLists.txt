add_executable(unittest EXCLUDE_FROM_ALL
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
  structures/testimage.cpp
  structures/testimagingtable.cpp
  ${WSCLEANFILES}
)

add_executable(integrationtest EXCLUDE_FROM_ALL
  test.cpp
  integration/tbasics.cpp
  integration/tveladeconvolution.cpp
  ${WSCLEANFILES}
)

set(TEST_LIBRARIES
  pybind11::embed
  ${ALL_LIBRARIES}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

target_link_libraries(unittest PRIVATE ${TEST_LIBRARIES})
target_link_libraries(integrationtest PRIVATE ${TEST_LIBRARIES})

add_custom_target(checkunit
  DEPENDS unittest
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMAND unittest
)

add_custom_target(checkintegration
  DEPENDS integrationtest
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_mwa_ms.sh
  COMMAND integrationtest
)

add_custom_target(check
  DEPENDS checkunit;checkintegration
)

add_custom_target(coverage
  DEPENDS check
  COMMAND gcovr -r .. -e '.*/tests/.*' -e '.*/CompilerIdCXX/.*'
)
