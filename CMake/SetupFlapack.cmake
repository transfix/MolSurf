#
# This macro is for setting up a sub-project to use GMP
#

macro(SetupFlapack TargetName)

  # CGAL provides its own Lapack search, which will be used by default.
  # Since CGAL requires some variables to be set, find CGAL to set
  # the variables before setting up Lapack.
  #find_package(CGAL)

  find_package(FlapackLocal)
  if(FLAPACK_LIB)
  target_link_libraries(${TargetName} ${FLAPACK_LIB} gfortran)
  endif(FLAPACK_LIB)
  if(FBLAS_LIB)
  target_link_libraries(${TargetName} ${FBLAS_LIB} gfortran)
  endif(FBLAS_LIB)

  message(Target: ${TargetName})
  message(FLAPACK_LIB: ${FLAPACK_LIB})
  message(FBLAS_LIB: ${FBLAS_LIB})
endmacro(SetupFlapack)
