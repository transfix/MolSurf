cmake_minimum_required(VERSION 3.16)
project(MolSurf
  VERSION 2.0.0
  DESCRIPTION "Volumetric image visualization and processing"
  LANGUAGES C CXX)

# C++ standard. C++20 is required because the libcvc dependency (pulled in
# below via FetchContent) uses std::format, which is a C++20 library feature.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Prefer the upstream BoostConfig.cmake shipped with modern Boost (>= 1.70).
# CMake 3.30 removed the legacy FindBoost module; Homebrew's Boost 1.90 only
# ships BoostConfig. Setting CMP0167 to NEW makes find_package(Boost) consume
# BoostConfig directly instead of falling back to the legacy finder.
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

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

# ── libcvc-deps prebuilt dependency bundle ───────────────────────────
# libcvc-deps (https://github.com/transfix/libcvc-deps) ships a single
# CMAKE_PREFIX_PATH-shaped tree containing Boost, FFTW, NFFT3, CGAL,
# GMP/MPFR, LAPACK/BLAS, and the rest of the libcvc/volrover3 dependency
# graph for Linux / macOS / Windows × Debug / Release. Pulling it down
# here means every Setup<Foo>.cmake (and the libcvc find below) sees a
# matching set of dependencies without per-platform apt/brew/vcpkg
# scaffolding in CI. Set -DLIBCVC_DEPS_SKIP=ON to fall back to the host
# package manager (useful for local dev when the bundle is already
# installed elsewhere).
set(LIBCVC_DEPS_VERSION "1.0.2" CACHE STRING
    "libcvc-deps release version to download")
option(LIBCVC_DEPS_SKIP
    "Do not download libcvc-deps; rely on system/CMAKE_PREFIX_PATH deps" OFF)

if(NOT LIBCVC_DEPS_SKIP)
  set(_dep_plat "")
  if(WIN32)
    set(_dep_plat "windows-x86_64")
    set(_dep_ext  "zip")
  elseif(APPLE)
    set(_dep_plat "macos-arm64")
    set(_dep_ext  "zip")
  elseif(UNIX)
    set(_dep_plat "linux-x86_64")
    set(_dep_ext  "tar.gz")
  endif()
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(_dep_cfg "debug")
  else()
    set(_dep_cfg "release")
  endif()
  if(_dep_plat)
    set(_dep_stem    "libcvc-deps-${LIBCVC_DEPS_VERSION}-${_dep_plat}-${_dep_cfg}")
    set(_dep_archive "${_dep_stem}.${_dep_ext}")
    set(_dep_url     "https://github.com/transfix/libcvc-deps/releases/download/v${LIBCVC_DEPS_VERSION}/${_dep_archive}")
    set(_dep_stage   "${CMAKE_BINARY_DIR}/_deps/libcvc-deps")
    set(_dep_root    "${_dep_stage}/${_dep_stem}")
    set(_dep_marker  "${_dep_root}/include")
    if(NOT EXISTS "${_dep_marker}")
      file(MAKE_DIRECTORY "${_dep_stage}")
      message(STATUS "Downloading libcvc-deps: ${_dep_url}")
      file(DOWNLOAD "${_dep_url}" "${_dep_stage}/${_dep_archive}"
           STATUS _dep_dl SHOW_PROGRESS TLS_VERIFY ON)
      list(GET _dep_dl 0 _dep_dl_code)
      if(_dep_dl_code EQUAL 0)
        execute_process(
          COMMAND ${CMAKE_COMMAND} -E tar xf "${_dep_archive}"
          WORKING_DIRECTORY "${_dep_stage}"
          RESULT_VARIABLE _dep_extract)
        if(NOT _dep_extract EQUAL 0)
          message(WARNING "Failed to extract ${_dep_archive} (exit ${_dep_extract})")
        endif()
      else()
        message(WARNING "libcvc-deps download failed: ${_dep_dl}")
      endif()
    endif()
    if(EXISTS "${_dep_marker}")
      list(PREPEND CMAKE_PREFIX_PATH "${_dep_root}")
      set(libcvc-deps_ROOT "${_dep_root}" CACHE PATH "" FORCE)
      message(STATUS "Using libcvc-deps ${LIBCVC_DEPS_VERSION} from ${_dep_root}")
    else()
      message(WARNING
        "libcvc-deps bundle not available; falling back to system packages on CMAKE_PREFIX_PATH")
    endif()
  endif()
endif()

# Dependency discovery macros. Each Setup<Foo>.cmake calls find_package
# against the libcvc-deps prefix prepended above (or the host package
# manager when LIBCVC_DEPS_SKIP=ON).
include(SetupBoost)
include(SetupFFTW)
include(SetupNFFT)
include(SetupLapack)
include(SetupCGAL)

# ── libcvc dependency ────────────────────────────────────────────────
# Resolution order, falling through on failure:
#   1. find_package(cvc CONFIG) against a pre-installed libcvc on
#      CMAKE_PREFIX_PATH (fastest, no download or build).
#   2. Download the prebuilt release tarball matching this platform /
#      build type and find_package against the extracted tree.
#   3. FetchContent_Declare(GIT_REPOSITORY ...) and build from source.
#
# MolSurf only needs the cvc::cvc target, which is shipped by every
# published libcvc release archive, so strategies 1 and 2 normally
# succeed and the source build (which pulls Boost/HDF5/CGAL/FFTW
# rediscovery) is avoided. Override with -DLIBCVC_FORCE_SOURCE=ON if
# you need to test against an unreleased libcvc commit.
set(LIBCVC_TAG     "v3.2.1" CACHE STRING "libcvc release tag to consume")
set(LIBCVC_VERSION "3.2.0"  CACHE STRING "libcvc package version to look up")
set(LIBCVC_GIT_SHA "beb23324b25208dba62d8bfde29ea8058ac5bcd7"
    CACHE STRING "libcvc commit SHA used when falling through to source build")
option(LIBCVC_FORCE_SOURCE "Skip prebuilt libcvc and always build from source" OFF)

include(FetchContent)

# Cache options consumed by libcvc's own CMakeLists when we fall through
# to the source build. Harmless when find_package succeeds.
set(CVC_BUILD_TESTS     OFF CACHE BOOL "" FORCE)
set(CVC_ENABLE_CUDA     OFF CACHE BOOL "" FORCE)
set(CVC_ENABLE_COVERAGE OFF CACHE BOOL "" FORCE)
set(CVC_BUILD_VOLROVER3 OFF CACHE BOOL "" FORCE)

set(_libcvc_usable FALSE)

# 1. Pre-installed package
if(NOT LIBCVC_FORCE_SOURCE)
  find_package(cvc ${LIBCVC_VERSION} CONFIG QUIET)
  if(cvc_FOUND)
    message(STATUS "Using pre-installed libcvc ${cvc_VERSION}")
    set(_libcvc_usable TRUE)
  endif()
endif()

# 2. Prebuilt release tarball
if(NOT _libcvc_usable AND NOT LIBCVC_FORCE_SOURCE)
  if(WIN32)
    set(_libcvc_plat "windows-x86_64")
    set(_libcvc_ext  "zip")
  elseif(APPLE)
    set(_libcvc_plat "macos-arm64")
    set(_libcvc_ext  "zip")
  elseif(UNIX)
    set(_libcvc_plat "linux-x86_64")
    set(_libcvc_ext  "tar.gz")
  endif()
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(_libcvc_cfg "debug")
  else()
    set(_libcvc_cfg "release")
  endif()
  if(DEFINED _libcvc_plat)
    set(_libcvc_archive "libcvc-${LIBCVC_VERSION}-${_libcvc_plat}-${_libcvc_cfg}.${_libcvc_ext}")
    set(_libcvc_url     "https://github.com/transfix/libcvc/releases/download/${LIBCVC_TAG}/${_libcvc_archive}")
    set(_libcvc_stage   "${CMAKE_BINARY_DIR}/_deps/libcvc-prebuilt")
    set(_libcvc_cfgfile "${_libcvc_stage}/lib/cmake/cvc/cvcConfig.cmake")
    if(NOT EXISTS "${_libcvc_cfgfile}")
      message(STATUS "Downloading prebuilt libcvc: ${_libcvc_url}")
      file(MAKE_DIRECTORY "${_libcvc_stage}")
      file(DOWNLOAD "${_libcvc_url}" "${_libcvc_stage}/${_libcvc_archive}"
           STATUS _libcvc_dl TLS_VERIFY ON)
      list(GET _libcvc_dl 0 _libcvc_dl_code)
      if(_libcvc_dl_code EQUAL 0)
        execute_process(
          COMMAND ${CMAKE_COMMAND} -E tar xzf "${_libcvc_archive}"
          WORKING_DIRECTORY "${_libcvc_stage}"
          RESULT_VARIABLE _libcvc_extract)
        if(NOT _libcvc_extract EQUAL 0)
          message(STATUS "Failed to extract ${_libcvc_archive} (exit ${_libcvc_extract})")
        endif()
      else()
        message(STATUS "Prebuilt libcvc download failed: ${_libcvc_dl}")
      endif()
    endif()
    if(EXISTS "${_libcvc_cfgfile}")
      list(PREPEND CMAKE_PREFIX_PATH "${_libcvc_stage}")
      find_package(cvc ${LIBCVC_VERSION} CONFIG QUIET)
      if(cvc_FOUND)
        message(STATUS "Using prebuilt libcvc ${cvc_VERSION} from ${_libcvc_url}")
        set(_libcvc_usable TRUE)
      endif()
    endif()
  endif()
endif()

# 3. Source build (fallback)
if(NOT _libcvc_usable)
  message(STATUS "Building libcvc ${LIBCVC_TAG} from source")
  FetchContent_Declare(
    libcvc
    GIT_REPOSITORY https://github.com/transfix/libcvc.git
    GIT_TAG        ${LIBCVC_GIT_SHA}
  )
  FetchContent_MakeAvailable(libcvc)
endif()

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



