CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES ../../inc/ScoringFuncs/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${INCLUDE_FILES})

INCLUDE_DIRECTORIES(
  ../../inc
)

ADD_LIBRARY (ScoringFuncs STATIC ${SOURCE_FILES} ${INCLUDE_FILES})

SET(SF_LINK_LIBS ScoringFuncs GOAFileTypes PDBParser DPG UsefulMath Utility)
ADD_EXECUTABLE (SurfaceHydrophobicityTest SurfaceHydrophobicityTest.cpp)
TARGET_LINK_LIBRARIES(SurfaceHydrophobicityTest ${SF_LINK_LIBS})
