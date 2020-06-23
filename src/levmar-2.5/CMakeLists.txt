CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

FILE(GLOB SOURCE_FILES Axb.c expfit.c lmdb.c lmblec.c lmbleic.c lm.c lmdemo.c lmlec.c misc.c)
FILE(GLOB INCLUDE_FILES ../../inc/UsefulMath/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${HEADER_FILES})

INCLUDE_DIRECTORIES(
./
)

ADD_LIBRARY (levmar STATIC ${SOURCE_FILES} ${INCLUDE_FILES})

SetupNFFT(levmar)
