cmake_minimum_required(VERSION 2.8.8)

project(wsclean)

include(CMakeVersionInfo.txt)

#if(NOT CMAKE_BUILD_TYPE)
#  set(CMAKE_BUILD_TYPE "Release")
#endif()

# Using pybind11 requires using -fvisibility=hidden.
# See https://pybind11.readthedocs.io/en/stable/faq.html
add_compile_options(-O3 -Wall -DNDEBUG -fvisibility=hidden)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(STATUS "Debug build selected: setting linking flag --no-undefined")
  string(APPEND CMAKE_SHARED_LINKER_FLAGS " -Wl,--no-undefined")
endif()

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++17" COMPILER_HAS_CXX17)

option(DISABLE_CXX17 "Do not use C++17 for compilation" OFF)
option(PORTABLE "Generate portable code" OFF)
option(BUILD_PACKAGES "Build Debian packages" OFF)

if(COMPILER_HAS_CXX17 AND NOT DISABLE_CXX17)
  add_compile_options(-std=c++17)
  add_definitions(-DHAVE_WGRIDDER)
  set(WGRIDDER_FILES wgridder/wgriddingmsgridder.cpp
    external/wgridder/ducc0/infra/threading.cc
    wgridder/wgriddinggridder_simple.cpp)
  include_directories(external/wgridder)
else()
  add_compile_options(-std=c++11)
  message(WARNING "Selected compiler does not support the C++17 standard. This will disable the 'wgridder' gridder, which requires C++17.")
  set(WGRIDDER_FILES)
endif()

if(NOT PORTABLE)
  add_compile_options(-march=native)
endif(NOT PORTABLE)

# Find and include git submodules
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive --checkout
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

# Include aocommon
include_directories("${CMAKE_SOURCE_DIR}/external/aocommon/include/")

# Schaapcommon
add_subdirectory("${CMAKE_SOURCE_DIR}/external/schaapcommon")
include_directories("${CMAKE_SOURCE_DIR}/external/schaapcommon/include")

# Casacore has a separate CMake file in this directory
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

find_package(HDF5 COMPONENTS C CXX REQUIRED)

SET(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(Casacore REQUIRED COMPONENTS casa ms tables measures fits)
find_package(CFITSIO REQUIRED)

find_library(FFTW3_LIB fftw3 REQUIRED)
find_library(FFTW3_THREADS_LIB fftw3_threads REQUIRED)
find_library(FFTW3F_LIB fftw3f REQUIRED)
find_library(FFTW3F_THREADS_LIB fftw3f_threads REQUIRED)
find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h)

find_package(PythonLibs 3 REQUIRED)
find_package(PythonInterp REQUIRED)
message(STATUS "Using python version ${PYTHON_VERSION_STRING}")

# Include pybind11
add_subdirectory("${CMAKE_SOURCE_DIR}/external/pybind11")
include_directories(${pybind11_INCLUDE_DIR})

#Prevent accidentally finding old BoostConfig.cmake file from casapy
set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost COMPONENTS date_time filesystem system program_options REQUIRED)
find_library(PTHREAD_LIB pthread REQUIRED)
find_library(DL_LIB dl REQUIRED)

# Once we bump minimum CMake version to >= 3.2, we can use `find_package(GSL REQUIRED)`.
find_library(GSL_LIB NAMES gsl)
find_path(GSL_INCLUDE_DIR NAMES gsl/gsl_version.h)
find_library(GSL_CBLAS_LIB NAMES gslcblas)
if(NOT GSL_LIB OR NOT GSL_INCLUDE_DIR OR NOT GSL_CBLAS_LIB)
  message(FATAL_ERROR "GSL not found, but required to build WSClean!")
endif()

find_package(MPI)

if(MPI_FOUND)
  add_compile_options(${MPI_COMPILE_FLAGS})
  add_definitions(-DHAVE_MPI)
  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
  set(MPI_CPP_FILES scheduling/mpischeduler.cpp distributed/mpibig.cpp)
else(MPI_FOUND)
  message(WARNING "MPI not found, multi-processing executable wsclean-mp will not be build.")
  set(MPI_LIBRARIES "")
  set(MPI_CPP_FILES)
endif(MPI_FOUND)

find_path(EVERYBEAM_INCLUDE_DIR NAMES EveryBeam/load.h)
find_library(EVERYBEAM_LIB everybeam)

if(EVERYBEAM_INCLUDE_DIR AND EVERYBEAM_LIB)
  include_directories(${EVERYBEAM_INCLUDE_DIR})
  add_definitions(-DHAVE_EVERYBEAM)
  message(STATUS "EveryBeam library found.")
else()
  set(EVERYBEAM_LIB "")
  message("EveryBeam library not found: Functionality for calculating beams for telescopes such as LOFAR will not be available.")
endif(EVERYBEAM_INCLUDE_DIR AND EVERYBEAM_LIB)

find_package(IDGAPI NO_MODULE QUIET)

if(IDGAPI_LIBRARIES AND IDGAPI_INCLUDE_DIRS)
  set(IDG_FILES idg/idgmsgridder.cpp)
  include_directories(${IDGAPI_INCLUDE_DIRS})
  add_definitions(-DHAVE_IDG)
  message("Image domain gridder libraries found.")
else(IDGAPI_LIBRARIES AND IDGAPI_INCLUDE_DIRS)
  set(IDGAPI_LIBRARIES)
  set(IDG_FILES)
  message("Image domain gridder libraries NOT found. Experimental gridder will not be available.")
endif(IDGAPI_LIBRARIES AND IDGAPI_INCLUDE_DIRS)

include_directories(${CASACORE_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${CFITSIO_INCLUDE_DIR})
include_directories(${FFTW3_INCLUDE_DIR})
include_directories(${GSL_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${PYTHON_INCLUDE_DIRS})

include(CheckCXXSourceCompiles)

# GSL is required for WSClean, so always available
add_definitions(-DHAVE_GSL)

# The following stuff will set the "rpath" correctly, so that
# LD_LIBRARY_PATH doesn't have to be set.

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

configure_file("${PROJECT_SOURCE_DIR}/wscversion.h.in" "${PROJECT_BINARY_DIR}/wscversion.h")
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Add CPack directory if user wants to generate Debian packages
if(BUILD_PACKAGES)
  add_subdirectory(CPack)
endif()

add_library(wsclean-object OBJECT
  main/commandline.cpp
  main/progressbar.cpp
  main/stopwatch.cpp
  main/wsclean.cpp
  main/settings.cpp
  deconvolution/casamaskreader.cpp
  deconvolution/componentlist.cpp
  deconvolution/deconvolution.cpp
  deconvolution/deconvolutionalgorithm.cpp
  deconvolution/genericclean.cpp
  deconvolution/imageset.cpp
  deconvolution/moresane.cpp
  deconvolution/paralleldeconvolution.cpp
  deconvolution/peakfinder.cpp
  deconvolution/pythondeconvolution.cpp
  deconvolution/simpleclean.cpp
  deconvolution/spectralfitter.cpp
  deconvolution/subminorloop.cpp
  gridding/directmsgridder.cpp
  gridding/msgridderbase.cpp
  gridding/wsmsgridder.cpp
  gridding/wstackinggridder.cpp
  idg/averagebeam.cpp
  interface/wscleaninterface.cpp
  io/facetreader.cpp
  io/logger.cpp
  io/parsetreader.cpp
  io/wscfitswriter.cpp
  iuwt/imageanalysis.cpp
  iuwt/iuwtdecomposition.cpp
  iuwt/iuwtdeconvolutionalgorithm.cpp
  iuwt/iuwtmask.cpp
  math/fftconvolver.cpp
  math/fftresampler.cpp
  math/gaussianfitter.cpp
  math/imageoperations.cpp
  math/modelrenderer.cpp
  math/nlplfitter.cpp
  math/polynomialchannelfitter.cpp
  math/polynomialfitter.cpp
  math/rmsimage.cpp
  model/model.cpp
  msproviders/averagingmsrowprovider.cpp
  msproviders/contiguousms.cpp
  msproviders/msdatadescription.cpp
  msproviders/directmsrowprovider.cpp
  msproviders/msprovider.cpp
  msproviders/msrowprovider.cpp
  msproviders/partitionedms.cpp
  msproviders/synchronizedms.cpp
  multiscale/multiscalealgorithm.cpp
  multiscale/multiscaletransforms.cpp
  multiscale/threadeddeconvolutiontools.cpp
  scheduling/griddingtaskmanager.cpp
  scheduling/griddingresult.cpp
  scheduling/griddingtask.cpp
  scheduling/metadatacache.cpp
  scheduling/threadedscheduler.cpp
  structures/image.cpp
  structures/imageweights.cpp
  structures/imagingtable.cpp
  structures/imagingtableentry.cpp
  structures/msselection.cpp
  structures/multibanddata.cpp
  structures/observationinfo.cpp
  structures/primarybeam.cpp
  system/fftwmanager.cpp
  system/system.cpp
  ${IDG_FILES} ${MPI_CPP_FILES} ${WGRIDDER_FILES})


# A number of files perform the 'core' high-performance floating point
# operations. In these files, NaNs are avoided and thus -ffast-math is
# allowed. Note that visibilities can be NaN hence this can not be turned
# on for all files.
set_source_files_properties(
	deconvolution/subminorloop.cpp
	deconvolution/genericclean.cpp
	deconvolution/imageset.cpp
	deconvolution/simpleclean.cpp
	deconvolution/spectralfitter.cpp
	multiscale/multiscalealgorithm.cpp
	multiscale/multiscaletransforms.cpp
	multiscale/threadeddeconvolutiontools.cpp
	wsclean/wstackinggridder.cpp
	wgridder/wgriddinggridder_simple.cpp
	PROPERTIES COMPILE_FLAGS -ffast-math)

set_property(TARGET wsclean-object PROPERTY POSITION_INDEPENDENT_CODE 1)

set(WSCLEANFILES $<TARGET_OBJECTS:wsclean-object>)

set(ALL_LIBRARIES ${CASACORE_LIBRARIES} ${FFTW3_LIB} ${FFTW3_THREADS_LIB} ${FFTW3F_LIB} ${FFTW3F_THREADS_LIB} ${Boost_DATE_TIME_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${CFITSIO_LIBRARY} ${GSL_LIB} ${GSL_CBLAS_LIB} ${PTHREAD_LIB} ${IDGAPI_LIBRARIES} ${HDF5_LIBRARIES} ${MPI_LIBRARIES} ${EVERYBEAM_LIB} schaapcommon)

add_library(wsclean-lib STATIC ${WSCLEANFILES})
target_link_libraries(wsclean-lib PRIVATE pybind11::embed ${ALL_LIBRARIES})
set_target_properties(wsclean-lib PROPERTIES OUTPUT_NAME wsclean)
set_target_properties(wsclean-lib PROPERTIES SOVERSION ${WSCLEAN_VERSION_SO})

add_library(wsclean-shared SHARED ${WSCLEANFILES})
target_link_libraries(wsclean-shared PRIVATE pybind11::embed ${ALL_LIBRARIES})
set_target_properties(wsclean-shared PROPERTIES SOVERSION ${WSCLEAN_VERSION_SO})

add_executable(wsclean main/main.cpp)
target_link_libraries(wsclean wsclean-lib)

if(MPI_FOUND)
  add_executable(wsclean-mp
    distributed/wsclean-mp.cpp
    distributed/slave.cpp)
  target_link_libraries(wsclean-mp wsclean-lib)
endif(MPI_FOUND)

add_executable(chgcentre chgcentre/main.cpp chgcentre/progressbar.cpp structures/multibanddata.cpp)
target_link_libraries(chgcentre ${CASACORE_LIBRARIES} ${GSL_LIB} ${GSL_CBLAS_LIB} ${MPI_LIBRARIES})

add_executable(wsuvbinning EXCLUDE_FROM_ALL gridding/examples/wsuvbinning.cpp ${WSCLEANFILES})
target_link_libraries(wsuvbinning ${ALL_LIBRARIES})

add_executable(wspredictionexample EXCLUDE_FROM_ALL gridding/examples/wspredictionexample.cpp gridding/wstackinggridder.cpp system/fftwmanager.cpp)
target_link_libraries(wspredictionexample ${ALL_LIBRARIES})

add_executable(mscaleexample EXCLUDE_FROM_ALL multiscale/mscaleexample.cpp math/fftconvolver.cpp multiscale/multiscalealgorithm.cpp multiscale/multiscaletransforms.cpp deconvolution/simpleclean.cpp)
target_link_libraries(mscaleexample ${ALL_LIBRARIES})

install(TARGETS wsclean DESTINATION bin)
install(TARGETS wsclean-lib DESTINATION lib)
install(TARGETS chgcentre DESTINATION bin)
install(FILES interface/wscleaninterface.h DESTINATION include)

if(MPI_FOUND)
  install(TARGETS wsclean-mp DESTINATION bin)
endif(MPI_FOUND)

add_custom_target(sphinxdoc
  make html
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
  COMMENT "Generating documentation with Sphinx" VERBATIM)

# add target to generate API documentation with Doxygen
find_package(Doxygen)

if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doxygendoc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM)
  add_custom_target(doc DEPENDS doxygendoc sphinxdoc)
else(DOXYGEN_FOUND)
  message(STATUS "Doxygen not found: API documentation can not compiled.")
  add_custom_target(doc DEPENDS sphinxdoc)
endif(DOXYGEN_FOUND)

# Boost 1.59 introduced BOOST_TEST, which several tests use.
find_package(Boost 1.59.0 COMPONENTS unit_test_framework)
if(Boost_FOUND)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
  add_subdirectory(tests)
else()
  message("Boost testing framework not found (not required for wsclean: only required for running tests).")
endif()

get_directory_property(MAIN_COMPILE_OPTIONS COMPILE_OPTIONS)
string(REPLACE ";" " " MAIN_COMPILE_OPTIONS "${MAIN_COMPILE_OPTIONS}")
message(STATUS "Flags passed to C++ compiler: ${MAIN_COMPILE_OPTIONS}")
