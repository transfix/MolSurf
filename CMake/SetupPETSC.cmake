#
# This macro is for setting up a sub-project to use PETSC
#
macro(SetupPETSC TargetName)
  set(PETSC_FOUND FOUND)
  target_link_libraries(${TargetName} petsc mpich mpl rt)

  SetupFlapack(${TargetName})
endmacro(SetupPETSC)
