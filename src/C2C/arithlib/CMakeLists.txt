CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

ADD_SUBDIRECTORY (libdjvu++)

FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES ../../../inc/C2C/arithlib.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${HEADER_FILES})

#TODO: change cpp files to include with lib directory name
#e.g.: #include<UsefulMath/Ray.h> instead of #include<Ray.h>
#Then change the below to ../../inc
INCLUDE_DIRECTORIES(
  ../../../inc/
)

ADD_LIBRARY (arithlib STATIC ${SOURCE_FILES} ${INCLUDE_FILES})
