cmake_minimum_required(VERSION 3.16)
project(MolSurf
  VERSION 2.0.0
  DESCRIPTION "Volumetric image visualization and processing"
  LANGUAGES C CXX)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Make CGAL's exact constructions stable with GCC. Other compilers do not
# require this and would warn or fail on the flag.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  add_compile_options(-frounding-math)
endif()

# Honor environment-provided flags (CPPFLAGS/LDFLAGS) without trampling
# anything the user passed on the command line.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    Debug Release RelWithDebInfo MinSizeRel)
endif()

option(USING_MAKE_ENV_VARS_FOR_INSTALL
  "Use the LDFLAGS environment variable as a directory to look up libs for a bundle" ON)
mark_as_advanced(USING_MAKE_ENV_VARS_FOR_INSTALL)
if(USING_MAKE_ENV_VARS_FOR_INSTALL)
  set(ENV_LIB_DIR "$ENV{ENV}/lib")
endif()
if(NOT BUNDLE_LIB_SEARCH_DIRS)
  set(BUNDLE_LIB_SEARCH_DIRS "${ENV_LIB_DIR}")
endif()
message(STATUS "BUNDLE_LIB_SEARCH_DIRS: ${BUNDLE_LIB_SEARCH_DIRS}")

# CMake module path: only our own modules and whatever the user/toolchain
# provides. No more reaching into ${PROJECT_BINARY_DIR} for vendored libs.
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")

# Dependency discovery macros. Each Setup<Foo>.cmake calls find_package
# against the system / package manager (apt, brew, vcpkg). NFFT may fall
# back to FetchContent because it has no first-class package on every OS.
include(SetupBoost)
include(SetupFFTW)
include(SetupNFFT)
include(SetupLapack)
include(SetupCGAL)

#output 
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
mark_as_advanced(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)


# Optional Subversion revision tag (legacy; mostly for backward-compat banners).
option(USING_SUBVERSION "Use the svn revision for version info" OFF)
set(CVC_SVN_REVISION "0")
if(USING_SUBVERSION)
  find_package(Subversion)
  if(Subversion_FOUND)
    Subversion_WC_INFO(${CMAKE_SOURCE_DIR} MolSurf)
    set(CVC_SVN_REVISION ${MolSurf_WC_REVISION})
    message(STATUS "Building Revision: ${CVC_SVN_REVISION}")
  endif()
endif()

# Version is set by project() above; alias the components for downstream
# code / CPack metadata that still expects MolSurf_VERSION_*.
set(MolSurf_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(MolSurf_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(MolSurf_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(MolSurf_VERSION       ${PROJECT_VERSION})
message(STATUS "MolSurf version: ${MolSurf_VERSION}")

# OSX Bundle settings
IF(APPLE)
  SET(CMAKE_BUNDLE_NAME
    "MolSurf-${MolSurf_VERSION_MAJOR}.${MolSurf_VERSION_MINOR}.${MolSurf_VERSION_PATCH}")
  SET(CMAKE_BUNDLE_LOCATION "${CMAKE_INSTALL_PREFIX}")
  # make sure CMAKE_INSTALL_PREFIX ends in /
  STRING(LENGTH "${CMAKE_INSTALL_PREFIX}" LEN)
  MATH(EXPR LEN "${LEN} -1" )
  STRING(SUBSTRING "${CMAKE_INSTALL_PREFIX}" ${LEN} 1 ENDCH)
  IF(NOT "${ENDCH}" STREQUAL "/")
    SET(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/")
  ENDIF(NOT "${ENDCH}" STREQUAL "/")
  SET(CMAKE_INSTALL_PREFIX
    "${CMAKE_INSTALL_PREFIX}${CMAKE_BUNDLE_NAME}.app/Contents")
ENDIF(APPLE)



# ************* CPack configuration for MolSurf *************
# Include the packaging system for distributing binaries
if(WIN32)
  set(CPACK_GENERATOR "ZIP")
elseif(APPLE)
  set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY 1)
  set(CPACK_GENERATOR "ZIP")
else(WIN32)
  set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY 1)
  set(CPACK_GENERATOR "TGZ")
endif(WIN32)
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Volumetric image visualization and processing")
SET(CPACK_PACKAGE_VENDOR "Computational Visualization Center")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/README")
SET(CPACK_PACKAGE_VERSION_MAJOR "${MolSurf_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${MolSurf_VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${MolSurf_VERSION_PATCH}")
SET(CPACK_PACKAGE_VERSION "${MolSurf_VERSION}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "MolSurf ${MolSurf_VERSION_MAJOR}.${MolSurf_VERSION_MINOR}")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "MolSurf-${MolSurf_VERSION}")
SET(CPACK_PACKAGE_NAME "${CMAKE_PROJECT_NAME}")
# Installers for 32- vs. 64-bit CMake:
#  - Root install directory (displayed to end user at installer-run time)
#  - "NSIS package/display name" (text used in the installer GUI)
#  - Registry key used to store info about the installation
IF(CMAKE_CL_64)
  SET(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES64")
  SET(CPACK_NSIS_PACKAGE_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} (Win64)")
  SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${CPACK_PACKAGE_NAME} ${CPACK_PACKAGE_VERSION} (Win64)")
ELSE()
  SET(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES")
  SET(CPACK_NSIS_PACKAGE_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY}")
  SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} MolSurf")
  SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\cvcweb.ices.utexas.edu")
  SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\cvcweb.ices.utexas.edu")
  SET(CPACK_NSIS_CONTACT "transfix@ices.utexas.edu")
  SET(CPACK_NSIS_MODIFY_PATH ON)
  SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${CPACK_PACKAGE_NAME} ${CPACK_PACKAGE_VERSION}")
ENDIF()
IF(NOT DEFINED CPACK_SYSTEM_NAME)
  # make sure package is not Cygwin-unknown, for Cygwin just
  # cygwin is good for the system name
  IF("${CMAKE_SYSTEM_NAME}" STREQUAL "CYGWIN")
    SET(CPACK_SYSTEM_NAME Cygwin)
  ELSE("${CMAKE_SYSTEM_NAME}" STREQUAL "CYGWIN")
    SET(CPACK_SYSTEM_NAME ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  ENDIF("${CMAKE_SYSTEM_NAME}" STREQUAL "CYGWIN")
ENDIF(NOT DEFINED CPACK_SYSTEM_NAME)
IF(${CPACK_SYSTEM_NAME} MATCHES Windows)
  IF(CMAKE_CL_64)
    SET(CPACK_SYSTEM_NAME win64-x64)
  ELSE(CMAKE_CL_64)
    SET(CPACK_SYSTEM_NAME win32-x86)
  ENDIF(CMAKE_CL_64)
ENDIF(${CPACK_SYSTEM_NAME} MATCHES Windows)
IF(NOT DEFINED CPACK_PACKAGE_FILE_NAME)
  # if the CPACK_PACKAGE_FILE_NAME is not defined by the cache
  # default to source package - system, on cygwin system is not
  # needed
  IF(CYGWIN)
    SET(CPACK_PACKAGE_FILE_NAME "${CPACK_SOURCE_PACKAGE_FILE_NAME}")
  ELSE(CYGWIN)
    SET(CPACK_PACKAGE_FILE_NAME
      "${CPACK_SOURCE_PACKAGE_FILE_NAME}-${CPACK_SYSTEM_NAME}")
  ENDIF(CYGWIN)
ENDIF(NOT DEFINED CPACK_PACKAGE_FILE_NAME)


SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Alex Rand") #required
SET(CPACK_PACKAGE_CONTACT "arand@ices.utexas.edu")
SET(CPACK_STRIP_FILES ON)
SET(CPACK_SOURCE_STRIP_FILES "")
SET(CPACK_PACKAGE_EXECUTABLES "MolSurf" "MolSurf")
SET(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_SOURCE_DIR}/README")
include( CPack )
include(InstallRequiredSystemLibraries)

option(PRE_BUILD "Build vendored third-party dependencies (REMOVED)" OFF)
if(PRE_BUILD)
  message(FATAL_ERROR
    "PRE_BUILD=ON is no longer supported. The vendored third-party/ directory "
    "has been removed; install dependencies through your system package "
    "manager (apt / brew / vcpkg). See INSTALL for details.")
endif()

add_subdirectory(src)



