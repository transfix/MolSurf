CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

# CHA: we should be able to find third-party from our own build
#IF (HAVE_NFFT)
ADD_DEFINITIONS(-DWITH_CG)
#ELSE (HAVE_NFFT)
#ADD_DEFINITIONS(-UWITH_CG)
#ENDIF (HAVE_NFFT)

#IF (HAVE_CGAL)
#ADD_DEFINITIONS(-DWITH_CGAL)
#ELSE (DISABLE_CGAL)
#ADD_DEFINITIONS(-UWITH_CGAL)
#ENDIF (HAVE_CGAL)


#IF (HAVE_CGAL)
#set(LIBS "")
#find_package(CGAL)
#if(CGAL_FOUND)
#  include(${CGAL_USE_FILE})
#  # need to set the following flags in case CGAL has some special compiler needs for this compiler
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CGAL_CXX_FLAGS_INIT}")
#  set(LIBS ${LIBS} ${CGAL_LIBRARIES})
#  ADD_DEFINITIONS(-DWITH_CGAL)
#endif(CGAL_FOUND)
#endif(HAVE_CGAL)


FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES ../../inc/MolSurfAPI/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${INCLUDE_FILES})

INCLUDE_DIRECTORIES(
  ../../inc
)

ADD_LIBRARY (MolSurfAPI STATIC ${SOURCE_FILES} ${INCLUDE_FILES})

SetupBoost(MolSurfAPI)
SetupLapack(MolSurfAPI)
SetupCGAL(MolSurfAPI)
SetupNFFT(MolSurfAPI)

target_link_libraries(MolSurfAPI ${LIBS})

