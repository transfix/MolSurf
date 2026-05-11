# FindMPFR.cmake
#
# Locate a system MPFR installation. Depends on GMP being found first.

include(FindPackageHandleStandardArgs)

find_path(MPFR_INCLUDE_DIR
  NAMES mpfr.h
  PATH_SUFFIXES mpfr
  DOC "MPFR include directory")

find_library(MPFR_LIBRARY
  NAMES mpfr libmpfr
  DOC "MPFR library")

find_package_handle_standard_args(MPFR
  REQUIRED_VARS MPFR_INCLUDE_DIR MPFR_LIBRARY)

mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARY)

if(MPFR_FOUND AND NOT TARGET MPFR::mpfr)
  add_library(MPFR::mpfr UNKNOWN IMPORTED)
  set_target_properties(MPFR::mpfr PROPERTIES
    IMPORTED_LOCATION "${MPFR_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${MPFR_INCLUDE_DIR}")
  if(TARGET GMP::gmp)
    set_property(TARGET MPFR::mpfr APPEND PROPERTY
      INTERFACE_LINK_LIBRARIES GMP::gmp)
  endif()

  # Legacy aliases.
  set(MPFR_INCLUDE "${MPFR_INCLUDE_DIR}")
  set(MPFR_LIB     "${MPFR_LIBRARY}")
endif()
