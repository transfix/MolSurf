#
# This macro is for setting up a sub-project to use CBlas
#

macro(SetupCBlas TargetName)

  set(CBLAS_FOUND NOTFOUND)

  find_package(CBlas)

  MESSAGE(CBLAS_LIB: ${CBLAS_LIB})

  target_link_libraries(${TargetName} ${CBLAS_LIB})
  include_directories(${CBLAS_INCLUDE})

endmacro(SetupCBlas)
