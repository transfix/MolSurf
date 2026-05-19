# SetupLapack.cmake
#
# Locate a system LAPACK + BLAS installation using CMake's stock
# FindLAPACK / FindBLAS modules. This supports OpenBLAS (Debian/Ubuntu's
# `liblapack-dev`), Apple Accelerate (macOS), and the various vcpkg
# LAPACK ports on Windows.
#
# Exposes a phony target MolSurf::lapack that bundles ${LAPACK_LIBRARIES}
# (which already includes BLAS).

if(NOT MolSurf_Lapack_Found)
  # Prefer a direct find against CMAKE_PREFIX_PATH so we pick up the
  # libcvc-deps bundle (which ships libblas.so + liblapack.so +
  # libgfortran.so.5 under <root>/lib/) without relying on CMake's
  # FindBLAS autoprobe. That autoprobe runs check_function_exists()
  # which fails when libblas.so transitively NEEDs libgfortran.so.5
  # via $ORIGIN RPATH (the linker doesn't follow $ORIGIN), even
  # though the runtime resolution works fine.
  set(_lapack_hint_dirs "")
  if(libcvc-deps_ROOT)
    list(APPEND _lapack_hint_dirs "${libcvc-deps_ROOT}/lib")
  endif()
  foreach(_p IN LISTS CMAKE_PREFIX_PATH)
    list(APPEND _lapack_hint_dirs "${_p}/lib")
  endforeach()
  message(STATUS "SetupLapack: probing for blas/lapack under: ${_lapack_hint_dirs}")
  find_library(MOLSURF_BLAS_LIBRARY   NAMES blas openblas HINTS ${_lapack_hint_dirs})
  find_library(MOLSURF_LAPACK_LIBRARY NAMES lapack        HINTS ${_lapack_hint_dirs})

  if(MOLSURF_BLAS_LIBRARY AND MOLSURF_LAPACK_LIBRARY)
    set(LAPACK_LIBRARIES "${MOLSURF_LAPACK_LIBRARY};${MOLSURF_BLAS_LIBRARY}")
    set(BLAS_LIBRARIES   "${MOLSURF_BLAS_LIBRARY}")
    # Reference BLAS / LAPACK depend on the gfortran runtime; pull
    # it in when bundled alongside.
    find_library(MOLSURF_GFORTRAN_LIBRARY
      NAMES gfortran
      HINTS ${_lapack_hint_dirs})
    if(MOLSURF_GFORTRAN_LIBRARY)
      list(APPEND LAPACK_LIBRARIES "${MOLSURF_GFORTRAN_LIBRARY}")
    elseif(libcvc-deps_ROOT AND EXISTS "${libcvc-deps_ROOT}/lib/libgfortran.so.5")
      # libgfortran in the bundle ships only as libgfortran.so.5 (no
      # unversioned .so symlink). find_library won't pick that up, so
      # reference the SONAME directly.
      list(APPEND LAPACK_LIBRARIES "${libcvc-deps_ROOT}/lib/libgfortran.so.5")
    endif()
    message(STATUS "LAPACK: using ${MOLSURF_LAPACK_LIBRARY}")
    message(STATUS "BLAS:   using ${MOLSURF_BLAS_LIBRARY}")
  else()
    message(STATUS "SetupLapack: direct probe failed (blas=${MOLSURF_BLAS_LIBRARY} lapack=${MOLSURF_LAPACK_LIBRARY}); falling back to find_package")
    find_package(LAPACK REQUIRED)
    find_package(BLAS QUIET)
  endif()

  if(NOT TARGET MolSurf::lapack)
    add_library(MolSurf::lapack INTERFACE IMPORTED GLOBAL)
    set_property(TARGET MolSurf::lapack PROPERTY
      INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}")
  endif()
  set(MolSurf_Lapack_Found TRUE CACHE INTERNAL "")

  # Legacy aliases for any straggling references.
  set(LAPACK_LIB "${LAPACK_LIBRARIES}")
  set(BLAS_LIB   "${BLAS_LIBRARIES}")
endif()

macro(SetupLapack TargetName)
  if(TARGET MolSurf::lapack)
    target_link_libraries(${TargetName} PUBLIC MolSurf::lapack)
  endif()
endmacro()
