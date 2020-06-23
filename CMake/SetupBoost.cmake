#
# These macros are for setting up Boost for a target
#

macro(IncludeBoost)
  set(Boost_FOUND FOUND)
  set(Boost_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/third-party/boost)
endmacro(IncludeBoost)

macro(find_boost)
  set(Boost_FOUND FOUND)
  set(Boost_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/third-party/boost)
  set(Boost_components thread date_time regex filesystem system program_options)
endmacro(find_boost)

macro(SetupBoost TargetName)
  # message("Setting up boost for target ${TargetName}")
  set(Boost_FOUND FOUND)
  set(Boost_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/third-party/boost)
# if (NOT DEFINED Boost_USE_MULTITHREADED)
  set(Boost_LIBRARIES boost_thread-mt boost_date_time-mt boost_regex-mt boost_filesystem-mt boost_system-mt boost_program_options-mt)
  target_link_libraries(${TargetName} ${Boost_LIBRARIES}) 
endmacro(SetupBoost)

##
## These macros are for setting up Boost for a target
##
#
#macro(IncludeBoost)
# find_package(Boost 1.41.0)
# if(Boost_FOUND)
#   include_directories(${Boost_INCLUDE_DIRS})
#   message("Boost includes: ${Boost_INCLUDE_DIRS}")
# else(Boost_FOUND)
#   message(SEND_ERROR "If you're having trouble finding boost, set CMake variables "
#           "BOOST_INCLUDEDIR and BOOST_LIBRARYDIR to the appropriate paths")
# endif(Boost_FOUND)
#endmacro(IncludeBoost)
#
#macro(find_boost)
# # defaults
# if (NOT DEFINED Boost_USE_STATIC)
#   set(Boost_USE_STATIC OFF)
# endif()
# if (NOT DEFINED Boost_USE_MULTITHREADED)
#   set(Boost_USE_MULTITHREADED ON)
# endif()
#
# set(Boost_components ${ARGN})
# if(NOT Boost_components)
#   # default components if none are specified
#   set(Boost_components thread date_time regex filesystem system program_options)
# endif(NOT Boost_components)
#
# message("Boost_components: ${Boost_components}")
# set(Boost_DEBUG TRUE)
#
# find_package(Boost 1.41.0 COMPONENTS ${Boost_components})
# if(Boost_FOUND)
#   include_directories(${Boost_INCLUDE_DIRS})
#   message("Boost includes: ${Boost_INCLUDE_DIRS}")
# else(Boost_FOUND)
#   message(SEND_ERROR "If you're having trouble finding boost, "
#                       "set CMake variables BOOST_INCLUDEDIR and "
#                       "BOOST_LIBRARYDIR to the appropriate paths")
# endif(Boost_FOUND)
#
#endmacro(find_boost)
#
#macro(SetupBoost TargetName)
# message("Setting up boost for target ${TargetName}")
## message("CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")
#
# find_boost(${ARGN})
# if(Boost_FOUND)
##   include_directories(${Boost_INCLUDE_DIRS})
#   set(LINK_LIBS ${LINK_LIBS} ${Boost_LIBRARIES}) 
##   message("Boost includes: ${Boost_INCLUDE_DIRS}")
#   message("Boost libraries: ${Boost_LIBRARIES}")
## else(Boost_FOUND)
##   message(SEND_ERROR "${TargetName}: If you're having trouble finding boost, "
##                       "set CMake variables BOOST_INCLUDEDIR and "
##                       "BOOST_LIBRARYDIR to the appropriate paths")
# endif(Boost_FOUND)
# 
# target_link_libraries(${TargetName} ${LINK_LIBS})
#endmacro(SetupBoost)
#