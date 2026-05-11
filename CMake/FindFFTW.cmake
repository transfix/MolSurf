# FindFFTW.cmake
#
# Locate a system FFTW3 installation (double + threads) and define imported
# targets compatible with the package-manager layouts on Linux, macOS, and
# Windows (apt's libfftw3-dev, Homebrew's `fftw`, vcpkg's `fftw3`).
#
# Provides:
#   - Imported target FFTW3::fftw3        (double precision)
#   - Imported target FFTW3::fftw3_threads (double-precision threaded; optional)
#   - Imported target FFTW3::fftw3f       (single precision; optional)
#   - Variables FFTW_INCLUDE, FFTW_LIB    (legacy aliases for existing call sites)

include(FindPackageHandleStandardArgs)

find_path(FFTW_INCLUDE_DIR
  NAMES fftw3.h
  PATH_SUFFIXES fftw
  DOC "FFTW3 include directory")

find_library(FFTW_DOUBLE_LIB
  NAMES fftw3 libfftw3-3 libfftw3
  DOC "FFTW3 double-precision library")

find_library(FFTW_DOUBLE_THREADS_LIB
  NAMES fftw3_threads libfftw3_threads
  DOC "FFTW3 double-precision threaded library")

find_library(FFTW_SINGLE_LIB
  NAMES fftw3f libfftw3f-3 libfftw3f
  DOC "FFTW3 single-precision library")

find_package_handle_standard_args(FFTW
  REQUIRED_VARS FFTW_INCLUDE_DIR FFTW_DOUBLE_LIB)

mark_as_advanced(FFTW_INCLUDE_DIR FFTW_DOUBLE_LIB FFTW_DOUBLE_THREADS_LIB FFTW_SINGLE_LIB)

if(FFTW_FOUND)
  if(NOT TARGET FFTW3::fftw3)
    add_library(FFTW3::fftw3 UNKNOWN IMPORTED)
    set_target_properties(FFTW3::fftw3 PROPERTIES
      IMPORTED_LOCATION "${FFTW_DOUBLE_LIB}"
      INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIR}")
  endif()
  if(FFTW_DOUBLE_THREADS_LIB AND NOT TARGET FFTW3::fftw3_threads)
    add_library(FFTW3::fftw3_threads UNKNOWN IMPORTED)
    set_target_properties(FFTW3::fftw3_threads PROPERTIES
      IMPORTED_LOCATION "${FFTW_DOUBLE_THREADS_LIB}"
      INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIR}"
      INTERFACE_LINK_LIBRARIES FFTW3::fftw3)
  endif()
  if(FFTW_SINGLE_LIB AND NOT TARGET FFTW3::fftw3f)
    add_library(FFTW3::fftw3f UNKNOWN IMPORTED)
    set_target_properties(FFTW3::fftw3f PROPERTIES
      IMPORTED_LOCATION "${FFTW_SINGLE_LIB}"
      INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIR}")
  endif()

  # Legacy aliases for call sites that still reference FFTW_INCLUDE / FFTW_LIB.
  set(FFTW_INCLUDE "${FFTW_INCLUDE_DIR}")
  set(FFTW_LIB     "${FFTW_DOUBLE_LIB}")
  if(FFTW_DOUBLE_THREADS_LIB)
    list(APPEND FFTW_LIB "${FFTW_DOUBLE_THREADS_LIB}")
  endif()
endif()
