# SetupNFFT.cmake
#
# Discover NFFT3, falling back to a FetchContent autotools build when the
# system has FFTW3 but no NFFT3 package (typical on Windows + vcpkg and on
# macOS + Homebrew). NFFT depends on FFTW so SetupFFTW must run first.
#
# Cache options:
#   MOLSURF_NFFT_FORCE_FETCHCONTENT  : skip find_package; always build from source.
#   MOLSURF_NFFT_VERSION             : tarball version to fetch (default 3.5.3).
#   MOLSURF_NFFT_URL                 : override the source tarball URL.

option(MOLSURF_NFFT_FORCE_FETCHCONTENT
  "Always build NFFT from source via FetchContent instead of using a system install" OFF)
set(MOLSURF_NFFT_VERSION "3.5.3" CACHE STRING "NFFT source version to fetch when not found")
set(MOLSURF_NFFT_URL "" CACHE STRING "Override URL for the NFFT source tarball")

if(NOT MolSurf_NFFT_Resolved)
  if(NOT MOLSURF_NFFT_FORCE_FETCHCONTENT)
    find_package(NFFT QUIET)
  endif()

  if(NOT NFFT_FOUND)
    # System NFFT not available -- build a private copy. Autotools-based, so
    # we drive it through ExternalProject and re-expose it via an imported
    # INTERFACE target named NFFT::nfft3.
    if(CMAKE_SYSTEM_NAME STREQUAL "Windows" AND NOT CMAKE_HOST_UNIX)
      message(FATAL_ERROR
        "NFFT3 was not found and the FetchContent build requires an autotools "
        "toolchain (configure + make), which is not available on a plain "
        "Windows host. Install NFFT3 (e.g. through MSYS2 with "
        "`pacman -S mingw-w64-x86_64-nfft`) or configure with -DMOLSURF_NFFT_FORCE_FETCHCONTENT=OFF "
        "and supply -DNFFT_INCLUDE_DIR / -DNFFT_LIBRARY pointing at an existing build.")
    endif()

    include(ExternalProject)
    include(FetchContent)

    if(NOT MOLSURF_NFFT_URL)
      set(MOLSURF_NFFT_URL
        "https://www-user.tu-chemnitz.de/~potts/nfft/download/nfft-${MOLSURF_NFFT_VERSION}.tar.gz")
    endif()

    set(_nfft_prefix "${CMAKE_BINARY_DIR}/_deps/nfft")
    set(_nfft_install "${_nfft_prefix}/install")

    # FFTW must already be discovered (SetupFFTW is included before us).
    if(NOT FFTW_INCLUDE_DIR OR NOT FFTW_DOUBLE_LIB)
      message(FATAL_ERROR
        "NFFT FetchContent build requires FFTW3 to be already located. "
        "Include SetupFFTW before SetupNFFT in the top-level CMakeLists.txt.")
    endif()
    get_filename_component(_fftw_lib_dir "${FFTW_DOUBLE_LIB}" DIRECTORY)

    ExternalProject_Add(nfft_external
      PREFIX            "${_nfft_prefix}"
      URL               "${MOLSURF_NFFT_URL}"
      CONFIGURE_COMMAND <SOURCE_DIR>/configure
                          --prefix=${_nfft_install}
                          --enable-static
                          --disable-shared
                          --with-fftw3-includedir=${FFTW_INCLUDE_DIR}
                          --with-fftw3-libdir=${_fftw_lib_dir}
      BUILD_COMMAND     ${CMAKE_MAKE_PROGRAM}
      INSTALL_COMMAND   ${CMAKE_MAKE_PROGRAM} install
      BUILD_BYPRODUCTS  "${_nfft_install}/lib/libnfft3.a"
      LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON)

    file(MAKE_DIRECTORY "${_nfft_install}/include")
    add_library(NFFT::nfft3 STATIC IMPORTED GLOBAL)
    set_target_properties(NFFT::nfft3 PROPERTIES
      IMPORTED_LOCATION "${_nfft_install}/lib/libnfft3.a"
      INTERFACE_INCLUDE_DIRECTORIES "${_nfft_install}/include")
    add_dependencies(NFFT::nfft3 nfft_external)

    set(NFFT_INCLUDE "${_nfft_install}/include")
    set(NFFT_LIB     "${_nfft_install}/lib/libnfft3.a")
    set(NFFT_FOUND TRUE)
    message(STATUS "NFFT: not found via find_package -- building ${MOLSURF_NFFT_VERSION} from source")
  endif()

  set(MolSurf_NFFT_Resolved TRUE CACHE INTERNAL "")
endif()

macro(SetupNFFT TargetName)
  if(TARGET NFFT::nfft3)
    target_link_libraries(${TargetName} PUBLIC NFFT::nfft3)
  endif()
  # Modern NFFT releases (>=3.3) and the libnfft3-dev package no longer
  # install the legacy `nfft3util.h` header that NFFT's bundled fastsum
  # application -- from which MolSurf's FastSummation sources are derived --
  # still includes. We ship a vendored copy of the 3.2-era header under
  # inc/compat/ so the sources keep compiling against any NFFT version.
  target_include_directories(${TargetName} PRIVATE
    "${CMAKE_SOURCE_DIR}/inc/compat")
  # NFFT calls into FFTW so always pull FFTW in alongside.
  SetupFFTW(${TargetName})
endmacro()
