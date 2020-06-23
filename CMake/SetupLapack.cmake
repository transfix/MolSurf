#
# This macro is for setting up a sub-project to use GMP
#

macro(SetupLapack TargetName)

  # CGAL provides its own Lapack search, which will be used by default.
  # Since CGAL requires some variables to be set, find CGAL to set
  # the variables before setting up Lapack.
  #find_package(CGAL)

  find_package(LapackLocal)
  if(LAPACK_LIB)
  target_link_libraries(${TargetName} ${LAPACK_LIB} gfortran)
  endif(LAPACK_LIB)
  if(BLAS_LIB)
  target_link_libraries(${TargetName} ${BLAS_LIB} gfortran)
  endif(BLAS_LIB)

  message(Target: ${TargetName})
  message(LAPACK_LIB: ${LAPACK_LIB})
  message(BLAS_LIB: ${BLAS_LIB})
endmacro(SetupLapack)
