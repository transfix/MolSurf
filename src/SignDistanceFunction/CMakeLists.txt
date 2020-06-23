CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

SET (SRCS
    common.cpp
    compute.cpp
    head.cpp
    init.cpp
    main.cpp
    new_adjust.cpp
    octree.cpp
    propagate.cpp
    )

FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES ../../inc/SignDistanceFunction/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${HEADER_FILES})

INCLUDE_DIRECTORIES (
  ../../inc
  ../../inc/SignDistanceFunction
  ../../inc/Blurmaps
  ../../inc/Geometry
 )

ADD_LIBRARY (SignDistanceFunction STATIC ${SOURCE_FILES} ${INCLUDE_FILES})
