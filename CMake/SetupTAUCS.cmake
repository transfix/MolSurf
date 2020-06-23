#
# This macro is for setting up a sub-project to use GMP
#

macro(SetupTAUCS TargetName)

  set(TAUCS_FOUND NOTFOUND)

  find_package(TAUCS)

  MESSAGE(TAUCS_LIB: ${TAUCS_LIB})

  target_link_libraries(${TargetName} ${TAUCS_LIB})
  include_directories(${TAUCS_INCLUDE})

endmacro(SetupTAUCS)
