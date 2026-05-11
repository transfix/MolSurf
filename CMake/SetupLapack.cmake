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
  find_package(LAPACK REQUIRED)
  find_package(BLAS QUIET)
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
