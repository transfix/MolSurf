# SetupNFFT.cmake
#
# Discover NFFT3 and expose it as the imported target `NFFT::nfft3`,
# plus the helper macro `SetupNFFT(<target>)` which wires NFFT and
# its FFTW dependency onto a consumer target.
#
# MolSurf now consumes NFFT 3.5.x from the libcvc-deps prebuilt
# bundle that the top-level CMakeLists.txt downloads on first
# configure (Linux uses the apt `libnfft3-dev` package inside the
# bundle; macOS and Windows use binaries built from source against
# the bundled FFTW3). The discovery here is therefore intentionally
# minimal: a CMake config-file lookup, a Find module lookup, and a
# manual find_path / find_library fallback.

if(NOT MolSurf_NFFT_Resolved)
  # 1) CMake package config (some NFFT distributions ship one).
  find_package(NFFT CONFIG QUIET)

  # 2) Find module (Find<Pkg>.cmake somewhere on CMAKE_MODULE_PATH).
  if(NOT NFFT_FOUND)
    find_package(NFFT QUIET)
  endif()

  # 3) Manual probe against CMAKE_PREFIX_PATH (libcvc-deps lives there).
  if(NOT TARGET NFFT::nfft3 AND NOT NFFT_FOUND)
    find_path(NFFT_INCLUDE_DIR
      NAMES nfft3.h
      PATH_SUFFIXES include)
    find_library(NFFT_LIBRARY
      NAMES nfft3 libnfft3
      PATH_SUFFIXES lib lib64)
    if(NFFT_INCLUDE_DIR AND NFFT_LIBRARY)
      add_library(NFFT::nfft3 UNKNOWN IMPORTED)
      set_target_properties(NFFT::nfft3 PROPERTIES
        IMPORTED_LOCATION "${NFFT_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${NFFT_INCLUDE_DIR}")
      set(NFFT_FOUND TRUE)
      message(STATUS "NFFT: found header ${NFFT_INCLUDE_DIR}/nfft3.h and library ${NFFT_LIBRARY}")
    endif()
  endif()

  if(NOT TARGET NFFT::nfft3 AND NOT NFFT_FOUND)
    message(FATAL_ERROR
      "NFFT3 was not found. MolSurf normally pulls NFFT3 from the "
      "libcvc-deps prebuilt bundle downloaded by CMakeLists.txt. If "
      "that download was disabled (-DLIBCVC_DEPS_SKIP=ON) install "
      "NFFT3 system-wide (apt: libnfft3-dev) or pass "
      "-DNFFT_INCLUDE_DIR=<dir> -DNFFT_LIBRARY=<lib> on the CMake "
      "command line.")
  endif()

  set(MolSurf_NFFT_Resolved TRUE CACHE INTERNAL "")
endif()

macro(SetupNFFT TargetName)
  if(TARGET NFFT::nfft3)
    target_link_libraries(${TargetName} PUBLIC NFFT::nfft3)
  elseif(NFFT_LIBRARY)
    target_link_libraries(${TargetName} PUBLIC "${NFFT_LIBRARY}")
    target_include_directories(${TargetName} PUBLIC "${NFFT_INCLUDE_DIR}")
  endif()
  # NFFT calls into FFTW so always pull FFTW in alongside.
  SetupFFTW(${TargetName})
endmacro()
