CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

#CHA: we should be able to find third-party libs from our own build
#IF (HAVE_LAPACK)
ADD_DEFINITIONS(-DHAVE_LAPACK)
#ELSE (HAVE_LAPACK)
#ADD_DEFINITIONS(-UHAVE_LAPACK)
#ENDIF (HAVE_LAPACK)

#IF (HAVE_NFFT)
ADD_DEFINITIONS(-DWITH_CG)
#ELSE (HAVE_NFFT)
#ADD_DEFINITIONS(-UWITH_CG)
#ENDIF (HAVE_NFFT)

#IF (HAVE_CGAL)
ADD_DEFINITIONS(-DWITH_CGAL)
#ELSE (DISABLE_CGAL)
#ADD_DEFINITIONS(-UWITH_CGAL)
#ENDIF (HAVE_CGAL)

IF (DISABLE_CGAL)
ADD_DEFINITIONS(-UWITH_CGAL)
ENDIF (DISABLE_CGAL)

IF(APPLE)
ADD_DEFINITIONS(-framework OpenGl)
ENDIF(APPLE)

find_package(OpenGL)


FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES ../../inc/MolSurf/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${HEADER_FILES})

INCLUDE_DIRECTORIES(
  ../../inc
)

IF(APPLE)
LINK_DIRECTORIES(/System/Library/Frameworks/OpenGL.framework/Libraries)
ENDIF(APPLE) 

SET(MOLSURF_LINK_LIBS MolSurfAPI LBIE_lib AdaptiveIsocontouring ASpline)

#IF (HAVE_NFFT)
SET(MOLSURF_LINK_LIBS ${MOLSURF_LINK_LIBS} OptimizerLib	ARM CoarseGrain Optimizer)
#ENDIF (HAVE_NFFT)

SET(MOLSURF_LINK_LIBS ${MOLSURF_LINK_LIBS} 
  Blurmaps 
  ByteOrder 
  C2C 
  arithlib
  libdjvu
  ComputeNormals
  SimpleVolumeData
  Contour
  Contouring
  Curvatures
  Decimation
  DPG
  FastMolSurf
  Geometry
  GeometryFileTypes
  GOAFileTypes
  Hausdorff
  HLevelSet
  CardinalBSpline
  LeastRMSD
  MolecularCharacteristics
#  MultiresolutionSurface
  OffsetSurface
  ColorSurfaceByFunction
  PDBParser
  Pocket
  PartitionMesh
  RepairSurfaceMesh
  SimpleCuration
  DynamicAdaptiveGrid
#  InterfaceStats
)

#IF(HAVE_LAPACK)
SET(MOLSURF_LINK_LIBS ${MOLSURF_LINK_LIBS} ${LAPACK_LIBRARIES})
#ENDIF(HAVE_LAPACK)

#IF(HAVE_CGAL)
SET(MOLSURF_LINK_LIBS ${MOLSURF_LINK_LIBS} PocketTunnel)
#ENDIF(HAVE_CGAL)

SET(MOLSURF_LINK_LIBS ${MOLSURF_LINK_LIBS} 
  SignDistanceFunction_v2
  SmoothingFunction
  SphericalPatchIntersections
  Summation
  UsefulMath
  Utility
  VolumeFileTypes)
#IF(HAVE_CGAL)
#SET(MOLSURF_LINK_LIBS ${MOLSURF_LINK_LIBS} CGAL)
#ENDIF(HAVE_CGAL)

#IF (HAVE_NFFT)
#SET(MOLSURF_LINK_LIBS ${MOLSURF_LINK_LIBS} kernelsarand FastSummation nfft3 fftw3 levmar)
find_package(FFTW)
SET(MOLSURF_LINK_LIBS ${MOLSURF_LINK_LIBS} FastSummation levmar ${FFTW_LIB})
#ENDIF (HAVE_NFFT)

SET(MOLSURF_LINK_LIBS ${MOLSURF_LINK_LIBS} 
  pthread
  #GL
  #GLU
  ${OPENGL_LIBRARIES}
)

ADD_EXECUTABLE (MolSurf MACOSX_BUNDLE 
	       ${SOURCE_FILES} ${INCLUDE_FILES})

SetupBoost(MolSurf)
SetupCGAL(MolSurf)
SetupNFFT(MolSurf)

TARGET_LINK_LIBRARIES(MolSurf ${MOLSURF_LINK_LIBS})


if(${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 2.4)
  if(APPLE)
    set_target_properties(MolSurf PROPERTIES
      OUTPUT_NAME ${CMAKE_BUNDLE_NAME})
  endif(APPLE)
  set(CMAKE_INSTALL_DESTINATION_ARGS
    BUNDLE DESTINATION "${CMAKE_BUNDLE_LOCATION}")
endif(${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 2.4)

install(TARGETS MolSurf
        RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
	LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)


# ARAND: eventually write the runMolSurf script...
if(UNIX AND NOT APPLE)
  install(PROGRAMS "${MolSurf_SOURCE_DIR}/bin/runMolSurf" DESTINATION bin)
endif(UNIX AND NOT APPLE)
if (APPLE)
  install(PROGRAMS "${MolSurf_SOURCE_DIR}/bin/runMolSurfOSX" DESTINATION ./)
endif (APPLE)

if (WIN32)
  install(PROGRAMS "${MolSurf_SOURCE_DIR}/bin/runMolSurfWin" DESTINATION ./)
endif (WIN32)

set(fixup_exe "\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/bin/MolSurf${CMAKE_EXECUTABLE_SUFFIX}")
if(APPLE)
  set(fixup_exe "\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_BUNDLE_NAME}.app")
endif(APPLE)


message("aaaa: usr/local/lib;${Boost_LIBRARY_DIRS};${BUNDLE_LIB_SEARCH_DIRS}")


install(CODE "
  include(\"${MolSurf_SOURCE_DIR}/CMake/BundleUtilities.cmake\")
  fixup_bundle(\"${fixup_exe}\" \"\" \"/usr/local/lib;/System/Library/Frameworks/Accelerate.framework;${Boost_LIBRARY_DIRS};${BUNDLE_LIB_SEARCH_DIRS}\")
")

