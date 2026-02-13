cmake_minimum_required(VERSION 3.23)

project(wsclean LANGUAGES CXX)

include(ExternalProject)
include(CMakeVersionInfo.txt)

add_compile_options(
  -O3
  -Wall
  -Wnon-virtual-dtor
  -Wzero-as-null-pointer-constant
  -Wduplicated-branches
  -Wundef
  -Wvla
  -Wpointer-arith
  -Wextra
  -Wno-unused-parameter
  -Wno-psabi)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(STATUS "Debug build selected: setting linking flag --no-undefined")
  string(APPEND CMAKE_SHARED_LINKER_FLAGS " -Wl,--no-undefined")
else()
  add_compile_options(-DNDEBUG)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS OFF)

option(PORTABLE "Build portable binaries (with slightly decreased performance)"
       OFF)
option(BUILD_PACKAGES "Build Debian packages" OFF)
option(ENABLE_WTOWERS
       "Enable compiling w-towers gridder (requires ska-sdp-func library)" OFF)

set(RADLER_CXX_FLAGS
    " "
    CACHE STRING "Extra C++ flags used for compiling Radler")

set(ExternalSubmoduleDirectories aocommon radler schaapcommon)

find_package(pybind11 CONFIG REQUIRED)

foreach(ExternalSubmodule ${ExternalSubmoduleDirectories})
  if(NOT EXISTS ${CMAKE_SOURCE_DIR}/external/${ExternalSubmodule})
    message(
      FATAL_ERROR
        "The external submodule '${ExternalSubmodule}' is missing in the external/ subdirectory. "
        "This is likely the result of downloading a git tarball without submodules. "
        "This is not supported: git tarballs do not provide the required versioning "
        "information for the submodules. Please perform a git clone of this repository."
    )
  endif()
endforeach()

# Find and include git submodules
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  # Update submodules as needed
  option(GIT_SUBMODULE "Check submodules during build" ON)
  if(GIT_SUBMODULE)
    message(STATUS "Submodule update")
    execute_process(
      COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive --checkout
              --depth 1
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(
        FATAL_ERROR
          "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules"
      )
    endif()
  endif()
endif()

# User may optionally set `TARGET_CPU` if `PORTABLE=OFF`
include(external/aocommon/CMake/SetTargetCPU.cmake)

# Include aocommon.
include_directories("${CMAKE_SOURCE_DIR}/external/aocommon/include/")

# Schaapcommon
set(SCHAAPCOMMON_MODULES ducc0 facets fitters h5parm math reordering)
add_subdirectory("${CMAKE_SOURCE_DIR}/external/schaapcommon")
include_directories("${CMAKE_SOURCE_DIR}/external/schaapcommon/include")
include_directories("${CMAKE_SOURCE_DIR}/external/schaapcommon/external")

# Casacore has a separate CMake file in this directory
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

# W-towers requires ska-sdp-func
if(ENABLE_WTOWERS)
  find_package(SKA_SDP_FUNC NAMES ska-sdp-func)
  if(SKA_SDP_FUNC_FOUND)
    if(EXISTS
       "${SKA_SDP_FUNC_INCLUDE_DIR}/ska-sdp-func/grid_data/sdp_grid_wstack_wtower.h"
    )
      message(
        STATUS
          "Enabling w-towers gridder; Using ska-sdp-func library found in ${SKA_SDP_FUNC_INCLUDE_DIR}"
      )
      add_definitions(-DBUILD_WTOWERS)
      include_directories(BEFORE SYSTEM ${SKA_SDP_FUNC_INCLUDE_DIR})
    else()
      message(
        FATAL_ERROR
          "Can't enable w-towers gridder: unsuitable ska-sdp-func library found in ${SKA_SDP_FUNC_INCLUDE_DIR}"
      )
    endif()
  else(SKA_SDP_FUNC_FOUND)
    message(
      FATAL_ERROR
        "Can't enable w-towers gridder: ska-sdp-func library not found.")
  endif(SKA_SDP_FUNC_FOUND)
else(ENABLE_WTOWERS)
  message(STATUS "Disabling w-towers gridder.")
endif(ENABLE_WTOWERS)

find_package(
  HDF5
  COMPONENTS CXX
  REQUIRED)
add_definitions(${HDF5_DEFINITIONS})

set(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(Casacore 3.6 REQUIRED COMPONENTS casa ms tables measures fits)

find_package(CFITSIO REQUIRED)
# Chgcentre uses a lapack function. EveryBeam provides lapack, but because
# EveryBeam is not required, chgcentre will fail to link when EveryBeam is not
# available, if LAPACK would not be explicitly required.
find_package(LAPACK REQUIRED)

find_library(FFTW3_LIB fftw3 REQUIRED HINTS ENV FFTW3_LIB)
find_library(FFTW3_THREADS_LIB fftw3_threads REQUIRED HINTS ENV FFTW3_LIB)
find_library(FFTW3F_LIB fftw3f REQUIRED HINTS ENV FFTW3_LIB)
find_library(FFTW3F_THREADS_LIB fftw3f_threads REQUIRED HINTS ENV FFTW3_LIB)
find_path(
  FFTW3_INCLUDE_DIR
  NAMES fftw3.h
  HINTS ENV FFTW3_INCLUDE)

find_package(
  Python3
  COMPONENTS Interpreter Development
  REQUIRED)

#Prevent accidentally finding old BoostConfig.cmake file from casapy
set(Boost_NO_BOOST_CMAKE ON)
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.30")
  # Disable warning; both old and new policies are fine
  cmake_policy(SET CMP0167 NEW)
endif()
#For providing compatibility with Boost 1.89+, the Boost COMPONENTS do not
#include date_time and system, which were removed in Boost 1.89. They have been
#header-only since Boost 1.73 / 1.69, so the minimum Boost version is 1.73.
find_package(Boost 1.73 REQUIRED COMPONENTS filesystem program_options)
find_package(Threads REQUIRED)
find_library(DL_LIB dl REQUIRED)

# Once we bump minimum CMake version to >= 3.2, we can use `find_package(GSL REQUIRED)`.
find_library(GSL_LIB NAMES gsl)
find_path(GSL_INCLUDE_DIR NAMES gsl/gsl_version.h)
find_library(GSL_CBLAS_LIB NAMES gslcblas)
if(NOT GSL_LIB
   OR NOT GSL_INCLUDE_DIR
   OR NOT GSL_CBLAS_LIB)
  message(FATAL_ERROR "GSL not found, but required to build WSClean!")
endif()

find_package(MPI)

if(MPI_FOUND)
  # FindMPI in CMake >= 3.10 provides MPI_CXX_COMPILE_OPTIONS, a list that can
  # be fed into add_compile_options(). In older versions of CMake, FindMPI
  # provides MPI_CXX_COMPILE_FLAGS, a white-space separated string that can be
  # fed into add_definitions().
  if(DEFINED MPI_CXX_COMPILE_OPTIONS)
    add_compile_options(${MPI_CXX_COMPILE_OPTIONS})
  else()
    add_definitions(${MPI_CXX_COMPILE_FLAGS})
  endif()
  add_definitions(-DHAVE_MPI)
  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
  set(MPI_CPP_FILES scheduling/mpischeduler.cpp
                    scheduling/mpiworkerscheduler.cpp distributed/mpibig.cpp)
else(MPI_FOUND)
  message(
    WARNING
      "MPI not found, multi-processing executable wsclean-mp will not be build."
  )
  set(MPI_LIBRARIES "")
  set(MPI_CPP_FILES)
endif(MPI_FOUND)

find_package(EveryBeam NO_MODULE QUIET)
if(${EVERYBEAM_FOUND})
  if(${EVERYBEAM_VERSION} VERSION_LESS "0.7.4" OR ${EVERYBEAM_VERSION}
                                                  VERSION_GREATER_EQUAL "0.9.0")
    message(
      FATAL_ERROR
        "EveryBeam version ${EVERYBEAM_VERSION} was found, but it is not compatible with this version of WSClean. This WSClean version requires EveryBeam version 0.7.x (x >= 4) or 0.8"
    )
  else()
    include_directories(${EVERYBEAM_INCLUDE_DIRS})
    add_definitions(-DHAVE_EVERYBEAM)
    message(STATUS "EveryBeam library found.")
  endif()
else(${EVERYBEAM_FOUND})
  message(
    STATUS
      "EveryBeam library not found: Functionality for calculating beams for telescopes such as LOFAR will not be available."
  )
  set(EVERYBEAM_LIB "")
endif(${EVERYBEAM_FOUND})

find_package(IDGAPI NO_MODULE QUIET)

if(IDGAPI_FOUND)
  if(${IDGAPI_VERSION} VERSION_LESS "0.8.1")
    message(STATUS "IDG version is incompatible")
    message(
      FATAL_ERROR
        "WSClean needs IDG version 0.8.x - with x >= 1 - but found version ${IDGAPI_VERSION}"
    )
  else()
    set(IDG_FILES idg/idgcommon.cpp idg/idgmsgridder.cpp
                  idg/facetidgmsgridder.cpp)
    include_directories(${IDGAPI_INCLUDE_DIRS})
    add_definitions(-DHAVE_IDG)
    message(STATUS "Image domain gridder libraries found.")
  endif()
else()
  set(IDGAPI_LIBRARIES)
  set(IDG_FILES)
  message(
    STATUS
      "Image domain gridder libraries NOT found. This gridder will not be available."
  )
endif()

include_directories(SYSTEM ${CASACORE_INCLUDE_DIRS})
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
include_directories(SYSTEM ${CFITSIO_INCLUDE_DIR})
include_directories(SYSTEM ${FFTW3_INCLUDE_DIR})
include_directories(SYSTEM ${GSL_INCLUDE_DIR})
include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})
include_directories(SYSTEM ${Python3_INCLUDE_DIRS})

include(CheckCXXSourceCompiles)

# GSL is required for WSClean, so always available
add_definitions(-DHAVE_GSL)

# Radler
set(RADLER_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/external/radler)
set(RADLER_LIB
    ${RADLER_INSTALL_PREFIX}/lib/libradler${CMAKE_STATIC_LIBRARY_SUFFIX})

# By default, CMake incorrectly forwards a prefix path with multiple items.
# Using an alternate LIST_SEPARATOR works around this behavior. See
# https://stackoverflow.com/questions/45414507/pass-a-list-of-prefix-paths-to-externalproject-add-in-cmake-args
# Use ':', since the CMAKE_PREFIX_PATH environment variable already uses it.
string(REPLACE ";" ":" ALT_CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}")
ExternalProject_Add(
  radler-module
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/radler
  BINARY_DIR ${CMAKE_BINARY_DIR}/external/radler/build
  BUILD_ALWAYS TRUE
  BUILD_BYPRODUCTS ${RADLER_LIB}
  LIST_SEPARATOR :
  CMAKE_ARGS
    -DCOMPILE_AS_EXTERNAL_PROJECT=On
    -DPORTABLE=${PORTABLE}
    -DTARGET_CPU=${TARGET_CPU}
    -DCMAKE_INCLUDE_PATH=${CMAKE_INCLUDE_PATH}
    -DCMAKE_LIBRARY_PATH=${CMAKE_LIBRARY_PATH}
    # The casacore root is set to the location where casacore was found.
    # Note that CASACORE_ROOT_DIR isn't set/updated by the find script, so we
    # shouldn't set Radler's CASACORE_ROOT_DIR to WSClean's one.
    -DCASACORE_ROOT_DIR=${CASACORE_INCLUDE_DIR}/..
    -DAOCOMMON_INCLUDE_DIR=${CMAKE_SOURCE_DIR}/external/aocommon/include/
    -DSCHAAPCOMMON_SOURCE_DIR=${CMAKE_SOURCE_DIR}/external/schaapcommon/
    -DBUILD_TESTING=OFF
    -DCMAKE_PREFIX_PATH=${ALT_CMAKE_PREFIX_PATH}:$ENV{CMAKE_PREFIX_PATH}
    -DCMAKE_INSTALL_PREFIX=${RADLER_INSTALL_PREFIX}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_CXX_FLAGS=${RADLER_CXX_FLAGS})

include_directories(${RADLER_INSTALL_PREFIX}/include)

# The following section will set the "rpath" correctly, so that
# LD_LIBRARY_PATH doesn't have to be set.

# Include GNUInstallDirs for CMAKE_INSTALL_FULL_LIBDIR
include(GNUInstallDirs)
# Use, i.e. don't skip the full RPATH for the build tree.
set(CMAKE_SKIP_BUILD_RPATH FALSE)
# When building, don't use the install RPATH already
# (but later on when installing).
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
# Add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# The RPATH to be used when installing, but only if it's not a system directory.
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
     "${CMAKE_INSTALL_FULL_LIBDIR}" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
endif("${isSystemDir}" STREQUAL "-1")

configure_file("${PROJECT_SOURCE_DIR}/wscversion.h.in"
               "${PROJECT_BINARY_DIR}/wscversion.h")
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Add CPack directory if user wants to generate Debian packages
if(BUILD_PACKAGES)
  add_subdirectory(CPack)
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # GCC 8.x requires linking with stdc++fs for the filesystem library
  # https://gcc.gnu.org/onlinedocs/gcc-9.1.0/libstdc++/manual/manual/status.html#status.iso.2017
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0)
    message(
      FATAL_ERROR "The GCC version is too old, upgrade to GCC 8.0 or newer")
  elseif(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    link_libraries(stdc++fs)
  endif()
endif()

add_library(
  wsclean-object OBJECT
  # These wgridder files are given first because they take a long time to compile. In multi-threaded compilation,
  # this has the advantage that it allows better parallelization.
  # Instantiation is split over multiple files to help make it more parallel and reduce memory requirements
  wgridder/wgriddingmsgridder.cpp
  wgridder/wgridder_float.cpp
  wgridder/wgridder_double.cpp
  main/commandline.cpp
  main/imageweightinitializer.cpp
  main/mshelper.cpp
  main/progressbar.cpp
  main/settings.cpp
  main/stopwatch.cpp
  main/wsclean.cpp
  gridding/averagecorrection.cpp
  gridding/directmsgridder.cpp
  gridding/h5solutiondata.cpp
  gridding/msgridder.cpp
  gridding/msgridderdata.cpp
  gridding/msgriddermanager.cpp
  gridding/msprovidercollection.cpp
  gridding/visibilitymodifier.cpp
  gridding/wsmsgridder.cpp
  gridding/wstackinggridder.cpp
  idg/averagebeam.cpp
  io/componentlistwriter.cpp
  io/facetreader.cpp
  io/parsetreader.cpp
  io/wscfitswriter.cpp
  math/imageoperations.cpp
  math/renderer.cpp
  math/subpixelrenderer.cpp
  math/fourierdomainrenderer.cpp
  math/tophatconvolution.cpp
  model/model.cpp
  msproviders/contiguousms.cpp
  msproviders/inmemoryprovider.cpp
  msproviders/msdatadescription.cpp
  msproviders/msprovider.cpp
  msproviders/operations.cpp
  msproviders/reorderedmsprovider.cpp
  msproviders/synchronizedms.cpp
  msproviders/timestepbuffer.cpp
  msproviders/msreaders/contiguousmsreader.cpp
  msproviders/msreaders/reorderedmsreader.cpp
  msproviders/msreaders/timestepbufferreader.cpp
  msproviders/rowproviders/averagingmsrowprovider.cpp
  msproviders/rowproviders/bdamsrowprovider.cpp
  msproviders/rowproviders/directmsrowprovider.cpp
  msproviders/rowproviders/mscolumns.cpp
  msproviders/rowproviders/msrowprovider.cpp
  msproviders/rowproviders/msrowproviderbase.cpp
  msproviders/rowproviders/msweightcolumn.cpp
  msproviders/rowproviders/noisemsrowprovider.cpp
  scheduling/griddingtask.cpp
  scheduling/griddingtaskfactory.cpp
  scheduling/griddingtaskmanager.cpp
  scheduling/griddingresult.cpp
  scheduling/metadatacache.cpp
  scheduling/threadedscheduler.cpp
  structures/facetutil.cpp
  structures/imageweights.cpp
  structures/imagingtable.cpp
  structures/imagingtableentry.cpp
  structures/inmemoryreordering.cpp
  structures/msselection.cpp
  structures/observationinfo.cpp
  structures/primarybeam.cpp
  structures/resources.cpp
  system/pythonfilepath.cpp
  ${IDG_FILES}
  ${MPI_CPP_FILES}
  $<$<BOOL:${ENABLE_WTOWERS}>:wtowers/wtowersmsgridder.cpp>)

# A number of files perform the 'core' high-performance floating point
# operations. In these files, NaNs are avoided and thus -ffast-math is
# allowed. Note that visibilities can be NaN hence this can not be turned
# on for all files.
set_source_files_properties(
  wsclean/directmsgridder.cpp wsclean/wstackinggridder.cpp
  wgridder/wgridder_float.cpp wgridder/wgridder_double.cpp
  gridding/msgriddermanager.cpp PROPERTIES COMPILE_FLAGS "-ffast-math")

# Add radler-module dependency to wsclean-object to ensure all
# wsclean targets can link against ${RADLER_LIB}
add_dependencies(wsclean-object radler-module)

set_property(TARGET wsclean-object PROPERTY POSITION_INDEPENDENT_CODE 1)

target_link_libraries(
  wsclean-object PRIVATE schaapcommon
                         $<$<BOOL:${ENABLE_WTOWERS}>:${SKA_SDP_FUNC_LIBRARY}>)

set(EXTERNAL_LIBRARIES
    ${CASACORE_LIBRARIES}
    ${FFTW3_LIB}
    ${FFTW3_THREADS_LIB}
    ${FFTW3F_LIB}
    ${FFTW3F_THREADS_LIB}
    Boost::filesystem
    Boost::program_options
    ${CFITSIO_LIBRARY}
    ${GSL_LIB}
    ${GSL_CBLAS_LIB}
    Threads::Threads
    ${IDGAPI_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${MPI_LIBRARIES}
    ${EVERYBEAM_LIB}
    schaapcommon
    pybind11::embed)

# Perform the BLAS check from aocommon
include("${CMAKE_CURRENT_LIST_DIR}/external/aocommon/CMake/CheckBLAS.cmake")
check_blas(LIBRARIES ${EXTERNAL_LIBRARIES})
# Radler lib cannot be checked on BLAS consistenty as it doesn't exist
# at cmake time. It needs to be first in the link order because Radler
# itself also depends on some of the libraries in the EXTERNAL_LIBRARIES list.
list(PREPEND EXTERNAL_LIBRARIES ${RADLER_LIB})

add_library(wsclean-lib STATIC)
target_link_libraries(wsclean-lib PRIVATE ${EXTERNAL_LIBRARIES} wsclean-object)
set_target_properties(wsclean-lib PROPERTIES OUTPUT_NAME wsclean)

add_library(wsclean-shared SHARED interface/wsclean.cpp)
target_link_libraries(wsclean-shared PRIVATE ${EXTERNAL_LIBRARIES}
                                             wsclean-object)
set_target_properties(wsclean-shared PROPERTIES OUTPUT_NAME wsclean)
# Define the headers. Because we use aocommon, we also have to install some files
# from aocommon. Only dependent aocommon files are installed.
target_sources(
  wsclean-shared
  INTERFACE FILE_SET
            wsclean_headers
            TYPE
            HEADERS
            BASE_DIRS
            "${CMAKE_CURRENT_SOURCE_DIR}/interface"
            FILES
            interface/wsclean.h
            interface/inmemoryms.h)
target_sources(
  wsclean-shared
  INTERFACE FILE_SET
            aocommon_headers
            TYPE
            HEADERS
            BASE_DIRS
            "${CMAKE_CURRENT_SOURCE_DIR}/external/aocommon/include"
            FILES
            external/aocommon/include/aocommon/banddata.h
            external/aocommon/include/aocommon/constants.h
            external/aocommon/include/aocommon/multibanddata.h
            external/aocommon/include/aocommon/optionalnumber.h
            external/aocommon/include/aocommon/polarization.h
            external/aocommon/include/aocommon/uvector.h
            external/aocommon/include/aocommon/vectormap.h
            external/aocommon/include/aocommon/io/serialistream.h
            external/aocommon/include/aocommon/io/serialostream.h)

# Make schaapcommon a PUBLIC dependency, since some WSClean headers use
# schaapcommon headers. An application/test linking to wsclean-lib or
# wsclean-shared thus needs schaapcommon as well.
target_link_libraries(wsclean-lib PUBLIC schaapcommon)
target_link_libraries(wsclean-shared PUBLIC schaapcommon)

add_executable(wsclean main/main.cpp)
target_link_libraries(wsclean wsclean-lib)

if(MPI_FOUND)
  add_executable(wsclean-mp distributed/wsclean-mp.cpp distributed/worker.cpp)
  target_link_libraries(wsclean-mp wsclean-lib)
endif(MPI_FOUND)

add_executable(chgcentre chgcentre/main.cpp main/progressbar.cpp)
target_link_libraries(chgcentre ${CASACORE_LIBRARIES} ${GSL_LIB}
                      ${GSL_CBLAS_LIB} ${MPI_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable(
  wspredictionexample EXCLUDE_FROM_ALL gridding/examples/wspredictionexample.cpp
                                       gridding/wstackinggridder.cpp)
target_link_libraries(wspredictionexample ${EXTERNAL_LIBRARIES})

install(TARGETS wsclean DESTINATION bin)
install(
  TARGETS wsclean-shared
  LIBRARY DESTINATION lib
          FILE_SET wsclean_headers
          DESTINATION include/wsclean
          FILE_SET aocommon_headers
          DESTINATION include)
install(TARGETS chgcentre DESTINATION bin)
if(MPI_FOUND)
  install(TARGETS wsclean-mp DESTINATION bin)
endif(MPI_FOUND)

add_custom_target(
  sphinxdoc
  make html
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
  COMMENT "Generating documentation with Sphinx"
  VERBATIM)

# add target to generate API documentation with Doxygen
find_package(Doxygen)

if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
                 ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(
    doxygendoc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)
  add_custom_target(doc DEPENDS doxygendoc sphinxdoc)
else(DOXYGEN_FOUND)
  message(STATUS "Doxygen not found: API documentation can not compiled.")
  add_custom_target(doc DEPENDS sphinxdoc)
endif(DOXYGEN_FOUND)

if(BUILD_TESTING)
  # Boost 1.59 introduced BOOST_TEST, which several tests use.
  find_package(Boost 1.59.0 COMPONENTS unit_test_framework)
  if(Boost_FOUND)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
    add_subdirectory(tests)
  else()
    message(
      STATUS
        "Boost testing framework not found (not required for wsclean: only required for running tests)."
    )
  endif()
endif()

get_directory_property(MAIN_COMPILE_OPTIONS COMPILE_OPTIONS)
string(REPLACE ";" " " MAIN_COMPILE_OPTIONS "${MAIN_COMPILE_OPTIONS}")
message(STATUS "Flags passed to C++ compiler: ${MAIN_COMPILE_OPTIONS}")
