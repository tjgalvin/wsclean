# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Ensure that we use /usr as install prefix when packaging.
set(CPACK_PACKAGING_INSTALL_PREFIX "/usr")

# We use git to determine the version string from the latest tag.
find_package(Git)
if(Git_FOUND)
  execute_process(
    COMMAND bash -c "${GIT_EXECUTABLE} describe --tags --dirty"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE VERSION_STRING
    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  # Strip leading non-numeric characters, like "v"
  string(REGEX REPLACE "^[^[0-9]]*" "" VERSION_STRING ${VERSION_STRING})
endif()
if(NOT VERSION_STRING)
  set(VERSION_STRING "0.0")
  message(WARNING "[CPack] Could not find tag to create version string, "
                  "using ${VERSION_STRING} instead.")
endif()

# Determine target CPU
include(${CMAKE_CURRENT_SOURCE_DIR}/DetermineTargetCPU.cmake)

# Set package name: encode target CPU, if defined, in package name.
if(DEFINED IDENTIFIED_TARGET_CPU)
  set(CPACK_PACKAGE_NAME "wsclean-${IDENTIFIED_TARGET_CPU}")
else()
  set(CPACK_PACKAGE_NAME "wsclean")
endif()

set(CPACK_PACKAGE_VENDOR "ASTRON")
set(CPACK_PACKAGE_VERSION "${VERSION_STRING}")

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "WSClean (w-stacking clean)")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://gitlab.com/aroffringa/wsclean")

set(CPACK_GENERATOR "DEB")
set(CPACK_SOURCE_GENERATOR "TGZ")

set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
if(DEFINED IDENTIFIED_TARGET_CPU)
  set(CPACK_DEBIAN_PACKAGE_DEPENDS
      "everybeam (>=0.1), idg-${IDENTIFIED_TARGET_CPU}-api (>=0.8)")
else()
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "everybeam (>=0.1), idg-api (>=0.8)")
endif()
set(CPACK_DEBIAN_PACKAGE_DESCRIPTION
    "WSClean (w-stacking clean) is a fast generic widefield imager.
 It implements several gridding algorithms and offers fully-automated
 multi-scale multi-frequency deconvolution.
 Manual: https://wsclean.readthedocs.io/")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "deb-packages@astron.nl")
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)

# Determine list of conflicting package names
if(DEFINED IDENTIFIED_TARGET_CPU)
  set(_conflicts wsclean)
else()
  set(_conflicts)
endif()
foreach(_cpu ${KNOWN_TARGET_CPUS})
  if(NOT "${_cpu}" STREQUAL "${IDENTIFIED_TARGET_CPU}")
    if("${_conflicts}" STREQUAL "")
      set(_conflicts "wsclean-${_cpu}")
    else()
      set(_conflicts "${_conflicts}, wsclean-${_cpu}")
    endif()
  endif()
endforeach()
set(CPACK_DEBIAN_PACKAGE_CONFLICTS "${_conflicts}")

include(CPack)

message(STATUS "Package name: ${CPACK_PACKAGE_NAME}")
message(STATUS "Package version: ${CPACK_PACKAGE_VERSION}")
