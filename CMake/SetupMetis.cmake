#
# This macro is for setting up a sub-project to use CBlas
#

macro(SetupMetis TargetName)

  set(METIS_FOUND NOTFOUND)

  find_package(Metis)

  MESSAGE(METIS_LIB: ${METIS_LIB})

  target_link_libraries(${TargetName} ${METIS_LIB})
  include_directories(${METIS_INCLUDE})

endmacro(SetupMetis)
