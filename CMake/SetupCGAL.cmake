#
# This macro is for setting up a sub-project to use CGAL
#

include(SetupGMP)
include(SetupMPFR)

#macro(SetupCGAL TargetName)
#  if(NOT DISABLE_CGAL)
#    message("Setting up CGAL for target ${TargetName}")
#    set(CGAL_FOUND FOUND)
#
#
#    # Note: I had to change CGAL_Core_LIBRARY in CGALConfig.cmake
#    # to have the filename of the core library.  I don't know why it wasn't
#    # included in the first place. -jme
##    find_package(CGAL COMPONENTS Core)
#
##    if(CGAL_FOUND)
#    #include(${CGAL_USE_FILE})
#    # need the following flags in case CGAL has some
#    # special compiler needs for this compiler
#    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CGAL_CXX_FLAGS_INIT}")
#
#    set(CGAL_LIBRARIES CGAL CGAL_Core)
##    set(CGAL_LIBRARIES CGAL CGAL_Core CGAL_ImageIO CGAL_Qt4)
#    message("CGAL libraries: ${CGAL_LIBRARIES}")
#    set(CGAL_3RD_PARTY_LIBRARIES gmp mpfr)
#    message("CGAL 3rd party libraries: ${CGAL_LIBRARIES}")
#
#    set(LIBS ${LIBS} ${CGAL_LIBRARIES})
#    set(LIBS ${LIBS} ${CGAL_3RD_PARTY_LIBRARIES})
#
#    add_definitions(-DUSING_CGAL)
#    if(CMAKE_COMPILER_IS_GNUCXX)
#      message("SetupCGAL: g++ detected, using -frounding-math")
#      add_definitions(-frounding-math)
#    endif(CMAKE_COMPILER_IS_GNUCXX)
#    target_link_libraries(${TargetName} ${LIBS})
#
#    #CHA: Setup GMP which used by CGAL for windows build
#    #Do we need to set this up on linux or mac as well?
#    #It seems the VolRover builds on the both withoout this.
#    if(WIN32)
#      SetupGMP(${TargetName})
#      SetupMPFR(${TargetName})
#    endif(WIN32)
#      
##    else(CGAL_FOUND)
##      message("${TargetName} is requesting CGAL but it isnt found on the system!")
##    endif(CGAL_FOUND)
#  endif(NOT DISABLE_CGAL)
#endmacro(SetupCGAL)

macro(SetupCGAL TargetName)
  if(NOT DISABLE_CGAL)
#    message("Setting up CGAL for target ${TargetName}")
    # Note: I had to change CGAL_Core_LIBRARY in CGALConfig.cmake to have
    # the filename of the core library.  I don't know why it wasn't included
    # in the first place. -jme
    find_package(CGAL COMPONENTS Core)
    #find_package(CGAL 4.0)
    if(CGAL_FOUND)
      include(${CGAL_USE_FILE})
      # need the following flags in case CGAL has some special compiler needs for this compiler
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CGAL_CXX_FLAGS_INIT}")

#      message("CGAL libraries: ${CGAL_LIBRARIES}")

      set(LIBS ${LIBS} ${CGAL_LIBRARIES})
      set(LIBS ${LIBS} ${CGAL_3RD_PARTY_LIBRARIES})

      add_definitions(-DUSING_CGAL)
      if(CMAKE_COMPILER_IS_GNUCXX)
#        message("SetupCGAL: g++ detected, using -frounding-math")
        add_definitions(-frounding-math)
      endif(CMAKE_COMPILER_IS_GNUCXX)
      target_link_libraries(${TargetName} ${LIBS})

      #CHA: Setup GMP which used by CGAL for windows build
      #Do we need to set this up on linux or mac as well?
      #It seems the VolRover builds on the both withoout this.
#      if(WIN32)
      SetupMPFR(${TargetName})
      SetupGMP(${TargetName})
#      endif(WIN32)
      
    else(CGAL_FOUND)
      message("${TargetName} is requesting CGAL but it isnt found on the system!")
    endif(CGAL_FOUND)
  endif(NOT DISABLE_CGAL)
endmacro(SetupCGAL)
