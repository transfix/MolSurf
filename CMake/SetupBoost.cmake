# SetupBoost.cmake
#
# Discover Boost via the system / package-manager installation (apt's
# libboost-all-dev, Homebrew's `boost`, vcpkg's `boost-*`, etc.) and attach
# the resulting imported targets to a caller's target.
#
# Usage:
#   SetupBoost(<target> [extra_components...])
#
# Default components match what MolSurf historically vendored.

set(MOLSURF_BOOST_DEFAULT_COMPONENTS
  thread date_time regex filesystem system program_options)

# Resolve Boost once per configure; cache for downstream calls.
if(NOT MolSurf_Boost_FOUND)
  find_package(Boost 1.66 REQUIRED
    COMPONENTS ${MOLSURF_BOOST_DEFAULT_COMPONENTS})
  set(MolSurf_Boost_FOUND TRUE CACHE INTERNAL "")
endif()

macro(SetupBoost TargetName)
  set(_components ${ARGN})
  if(NOT _components)
    set(_components ${MOLSURF_BOOST_DEFAULT_COMPONENTS})
  endif()
  foreach(_c IN LISTS _components)
    if(TARGET Boost::${_c})
      target_link_libraries(${TargetName} PUBLIC Boost::${_c})
    endif()
  endforeach()
  if(TARGET Boost::headers)
    target_link_libraries(${TargetName} PUBLIC Boost::headers)
  endif()
endmacro()
