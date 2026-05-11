# SetupMPFR.cmake
#
# Attach the system MPFR imported target to a caller's target.

if(NOT MolSurf_MPFR_Found)
  find_package(MPFR REQUIRED)
  set(MolSurf_MPFR_Found TRUE CACHE INTERNAL "")
endif()

macro(SetupMPFR TargetName)
  if(TARGET MPFR::mpfr)
    target_link_libraries(${TargetName} PUBLIC MPFR::mpfr)
  endif()
endmacro()
