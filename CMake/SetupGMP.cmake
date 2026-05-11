# SetupGMP.cmake
#
# Attach the system GMP imported targets to a caller's target.

if(NOT MolSurf_GMP_Found)
  find_package(GMP REQUIRED)
  set(MolSurf_GMP_Found TRUE CACHE INTERNAL "")
endif()

macro(SetupGMP TargetName)
  if(TARGET GMP::gmp)
    target_link_libraries(${TargetName} PUBLIC GMP::gmp)
  endif()
  if(TARGET GMP::gmpxx)
    target_link_libraries(${TargetName} PUBLIC GMP::gmpxx)
  endif()
endmacro()
