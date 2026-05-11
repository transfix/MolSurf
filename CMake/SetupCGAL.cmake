# SetupCGAL.cmake
#
# Locate a system CGAL installation (libcgal-dev / Homebrew `cgal` /
# vcpkg `cgal`) and expose it through CGAL's modern imported target
# CGAL::CGAL. GMP and MPFR are CGAL's required arithmetic dependencies
# and are pulled in alongside.
#
# Set DISABLE_CGAL=ON to turn the macro into a no-op (useful when a
# downstream component intentionally avoids CGAL on some platforms).

option(DISABLE_CGAL "Disable all CGAL-dependent components" OFF)

if(NOT DISABLE_CGAL AND NOT MolSurf_CGAL_Found)
  # CGAL configs require its arithmetic deps to be discoverable.
  include(SetupGMP)
  include(SetupMPFR)
  find_package(CGAL CONFIG REQUIRED COMPONENTS Core)
  set(MolSurf_CGAL_Found TRUE CACHE INTERNAL "")
endif()

macro(SetupCGAL TargetName)
  if(NOT DISABLE_CGAL)
    if(TARGET CGAL::CGAL)
      target_link_libraries(${TargetName} PUBLIC CGAL::CGAL)
    endif()
    if(TARGET CGAL::CGAL_Core)
      target_link_libraries(${TargetName} PUBLIC CGAL::CGAL_Core)
    endif()
    target_compile_definitions(${TargetName} PUBLIC USING_CGAL)
    # CGAL needs IEEE rounding to be honored; only meaningful with GCC.
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      target_compile_options(${TargetName} PRIVATE -frounding-math)
    endif()
    SetupMPFR(${TargetName})
    SetupGMP(${TargetName})
  endif()
endmacro()
