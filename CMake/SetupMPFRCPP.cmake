#
# This macro is for setting up a sub-project to use FFTW (float or double version)
#

macro(SetupMPFRCPP TargetName)

  set(MPFRCPP_FOUND NOTFOUND)

  find_package(MPFRCPP)

  MESSAGE(MPFRCPP_LIB: ${MPFRCPP_LIB})

  if(MPFRCPP_LIB)
  target_link_libraries(${TargetName} ${MPFRCPP_LIB})  
  endif(MPFRCPP_LIB)
  
  IF(MPFRCPP_INCLUDE)
    INCLUDE_DIRECTORIES( ${MPFRCPP_INCLUDE})
  ENDIF(MPFRCPP_INCLUDE)

endmacro(SetupMPFRCPP)
