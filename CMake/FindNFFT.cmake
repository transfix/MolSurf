# FindNFFT.cmake
#
# Locate NFFT3. Available natively on Debian/Ubuntu (libnfft3-dev) but
# absent from vcpkg and core Homebrew, so we may have to fall back to a
# FetchContent build (handled in SetupNFFT.cmake).

include(FindPackageHandleStandardArgs)

# Try pkg-config first when available (covers most Linux installs).
find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
  pkg_check_modules(PC_NFFT QUIET nfft3)
endif()

find_path(NFFT_INCLUDE_DIR
  NAMES nfft3.h
  HINTS ${PC_NFFT_INCLUDE_DIRS}
  PATH_SUFFIXES nfft nfft3
  DOC "NFFT3 include directory")

find_library(NFFT_LIBRARY
  NAMES nfft3 libnfft3
  HINTS ${PC_NFFT_LIBRARY_DIRS}
  DOC "NFFT3 library")

find_package_handle_standard_args(NFFT
  REQUIRED_VARS NFFT_INCLUDE_DIR NFFT_LIBRARY)

mark_as_advanced(NFFT_INCLUDE_DIR NFFT_LIBRARY)

if(NFFT_FOUND AND NOT TARGET NFFT::nfft3)
  add_library(NFFT::nfft3 UNKNOWN IMPORTED)
  set_target_properties(NFFT::nfft3 PROPERTIES
    IMPORTED_LOCATION "${NFFT_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${NFFT_INCLUDE_DIR}")

  # Legacy aliases for old call sites.
  set(NFFT_INCLUDE "${NFFT_INCLUDE_DIR}")
  set(NFFT_LIB     "${NFFT_LIBRARY}")
endif()
