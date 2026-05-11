# FindGMP.cmake
#
# Locate a system GMP installation (libgmp-dev / Homebrew gmp / vcpkg gmp).

include(FindPackageHandleStandardArgs)

find_path(GMP_INCLUDE_DIR
  NAMES gmp.h
  PATH_SUFFIXES gmp
  DOC "GMP include directory")

find_library(GMP_LIBRARY
  NAMES gmp libgmp
  DOC "GMP library")

find_library(GMPXX_LIBRARY
  NAMES gmpxx libgmpxx
  DOC "GMP C++ wrapper library")

find_package_handle_standard_args(GMP
  REQUIRED_VARS GMP_INCLUDE_DIR GMP_LIBRARY)

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARY GMPXX_LIBRARY)

if(GMP_FOUND)
  if(NOT TARGET GMP::gmp)
    add_library(GMP::gmp UNKNOWN IMPORTED)
    set_target_properties(GMP::gmp PROPERTIES
      IMPORTED_LOCATION "${GMP_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}")
  endif()
  if(GMPXX_LIBRARY AND NOT TARGET GMP::gmpxx)
    add_library(GMP::gmpxx UNKNOWN IMPORTED)
    set_target_properties(GMP::gmpxx PROPERTIES
      IMPORTED_LOCATION "${GMPXX_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}"
      INTERFACE_LINK_LIBRARIES GMP::gmp)
  endif()

  # Legacy aliases.
  set(GMP_INCLUDE "${GMP_INCLUDE_DIR}")
  set(GMP_LIB     "${GMP_LIBRARY}")
  set(GMPXX_LIB   "${GMPXX_LIBRARY}")
endif()
