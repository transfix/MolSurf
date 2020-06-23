#
# This macro is for setting up a sub-project to use FFTW (float or double version)
#

macro(SetupMPFR TargetName)
  find_package(MPFR)

  if(MPFR_LIB)
  target_link_libraries(${TargetName} ${MPFR_LIB})  
  endif(MPFR_LIB)
  
  IF(MPFR_INCLUDE)
    INCLUDE_DIRECTORIES( ${MPFR_INCLUDE})
  ENDIF(MPFR_INCLUDE)

#  set(MPFR_FOUND NOTFOUND)
#  find_package(MPFR)
#  MESSAGE(MPFR_LIB: ${MPFR_LIB})
#  target_link_libraries(${TargetName} ${MPFR_LIB})  
endmacro(SetupMPFR)
