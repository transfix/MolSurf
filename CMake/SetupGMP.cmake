#
# This macro is for setting up a sub-project to use GMP
#

macro(SetupGMP TargetName)
  find_package(GMP)
  if(GMP_LIB)
  target_link_libraries(${TargetName} ${GMP_LIB})
  endif(GMP_LIB)
  if(GMPXX_LIB)
  target_link_libraries(${TargetName} ${GMPXX_LIB})
  endif(GMPXX_LIB)
  
  IF(GMP_INCLUDE)
    INCLUDE_DIRECTORIES( ${GMP_INCLUDE})
  ENDIF(GMP_INCLUDE)

#  set(GMP_FOUND NOTFOUND)
#  find_package(GMP)
#  MESSAGE(GMP_LIB: ${GMP_LIB})
#  target_link_libraries(${TargetName} ${GMP_LIB})
endmacro(SetupGMP)
