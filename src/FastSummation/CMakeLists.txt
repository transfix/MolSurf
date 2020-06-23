CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

if(CMAKE_HOST_UNIX)
  # kind of a hack... using CPP flags (i.e. include directories)
  #     for the C build also
  SET (CMAKE_C_FLAGS $ENV{CPPFLAGS})
endif(CMAKE_HOST_UNIX)

FILE(GLOB SOURCE_FILES *.cpp *.c)
FILE(GLOB INCLUDE_FILES ../../inc/FastSummation/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${HEADER_FILES})

INCLUDE_DIRECTORIES(
  ../../inc
)

ADD_LIBRARY (FastSummation STATIC ${SOURCE_FILES} ${INCLUDE_FILES})

SetupNFFT(FastSummation)
