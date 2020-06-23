CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES ../../inc/PocketTunnel/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${INCLUDE_FILES})

INCLUDE_DIRECTORIES(
  ../../inc
)

ADD_DEFINITIONS(-DHUGE=1000000)

ADD_LIBRARY (PocketTunnel STATIC ${SOURCE_FILES} ${INCLUDE_FILES})

set(LIBS "")

#find_package(CGAL)
#if(CGAL_FOUND)
#  include(${CGAL_USE_FILE})
  # need the following flags in case CGAL has some special compiler needs for this compiler
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CGAL_CXX_FLAGS_INIT}")
#  set(LIBS ${LIBS} ${CGAL_LIBRARIES})
#else(CGAL_FOUND)
#  message(SEND_ERROR "PocketTunnel requires CGAL!")
#endif(CGAL_FOUND)

target_link_libraries(PocketTunnel
  ${LIBS}
)

SetupCGAL(PocketTunnel)
