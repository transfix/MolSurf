cmake_minimum_required (VERSION 2.6)
project (MolSurf)

# arand: either of the lines below eliminate an annoying
#        CGAL problem, 7-20-2010
#ADD_DEFINITIONS(-DCGAL_DISABLE_ROUNDING_MATH_CHECK)
ADD_DEFINITIONS(-frounding-math)

#set (DISABLE_BUILD_COARSEGRAIN TRUE)

# Global Variables
if(CMAKE_HOST_UNIX)
  # standard setup
  set (CMAKE_CXX_FLAGS "$ENV{CPPFLAGS} -O3")
  set (CMAKE_LD_FLAGS "$ENV{LDFLAGS}")

  # Debugging and profiling friendly flags:
  #set (CMAKE_CXX_FLAGS "$ENV{CPPFLAGS} -pg")
  #set (CMAKE_LD_FLAGS "$ENV{LDFLAGS} -pg")

  #set (CMAKE_CXX_FLAGS "-I/usr/include -O2")
  #set (CMAKE_LD_FLAGS "-L/usr/lib")
endif(CMAKE_HOST_UNIX)

SET(CGAL_DIR
  ${PROJECT_BINARY_DIR}/lib/CGAL)


option(USING_MAKE_ENV_VARS_FOR_INSTALL "If this is enabled, use the LDFLAGS environment variable as a directory to look up libs for a bundle" ON)
mark_as_advanced(USING_MAKE_ENV_VARS_FOR_INSTALL)
if(USING_MAKE_ENV_VARS_FOR_INSTALL)
  set(ENV_LIB_DIR "$ENV{ENV}/lib")
endif(USING_MAKE_ENV_VARS_FOR_INSTALL)

if(NOT BUNDLE_LIB_SEARCH_DIRS)
  set(BUNDLE_LIB_SEARCH_DIRS "${ENV_LIB_DIR}")
endif(NOT BUNDLE_LIB_SEARCH_DIRS)
message(STATUS "BUNDLE_LIB_SEARCH_DIRS: ${BUNDLE_LIB_SEARCH_DIRS}")


# CMake Modules
SET(CMAKE_MODULE_PATH
  ${PROJECT_BINARY_DIR}/lib/CGAL
  ${PROJECT_BINARY_DIR}/lib/cmake/lapack-3.4.2
  ${CMAKE_SOURCE_DIR}/CMake
  ${CMAKE_MODULE_PATH})


# include macro for setting up Boost
include( SetupBoost )
include( SetupFFTW )
include( SetupNFFT )
include( SetupLapack )
include( SetupCGAL )

#output 
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
mark_as_advanced(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)


# Set the revision number from svn.
option(USING_SUBVERSION "If you're using subversion, MolSurf can use the revision number for version info" OFF)
set(CVC_SVN_REVISION "0")
if(USING_SUBVERSION)
  find_package(Subversion)
  if(Subversion_FOUND)
    Subversion_WC_INFO(${CMAKE_SOURCE_DIR} MolSurf)
    set(CVC_SVN_REVISION ${MolSurf_WC_REVISION})
    message(STATUS "Building Revision: ${CVC_SVN_REVISION}")
  endif(Subversion_FOUND)
endif(USING_SUBVERSION)


# The MolSurf version number.
SET(MolSurf_VERSION_MAJOR 1)
SET(MolSurf_VERSION_MINOR 0)
SET(MolSurf_VERSION_PATCH ${CVC_SVN_REVISION})

# The full version string
SET(MolSurf_VERSION ${MolSurf_VERSION_MAJOR}.${MolSurf_VERSION_MINOR}.${MolSurf_VERSION_PATCH})
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

option(PRE_BUILD "Builds dependencies" ON)
if(PRE_BUILD)
  ADD_SUBDIRECTORY (third-party)
  MESSAGE(STATUS "*******************************************************")
  MESSAGE(STATUS "** If your expected target is not available, consider")
  MESSAGE(STATUS "** changing the PRE_BUILD option using ccmake.")
  MESSAGE(STATUS "*******************************************************")
else(PRE_BUILD)
  include_directories(${PROJECT_BINARY_DIR}/include)
  link_directories(${PROJECT_BINARY_DIR}/lib)
  ADD_SUBDIRECTORY (src)
  MESSAGE(STATUS "*******************************************************")
  MESSAGE(STATUS "** If needed dependencies are not available, consider")
  MESSAGE(STATUS "** changing the PRE_BUILD option using ccmake.")
  MESSAGE(STATUS "*******************************************************")
endif(PRE_BUILD)



