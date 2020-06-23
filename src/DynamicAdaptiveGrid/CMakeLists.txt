CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES ../../inc/DynamicAdaptiveGrid/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${HEADER_FILES})

INCLUDE_DIRECTORIES(
  ../../inc
)

ADD_LIBRARY (DynamicAdaptiveGrid STATIC ${SOURCE_FILES} ${INCLUDE_FILES})

SET(DAG_LINK_LIBS DynamicAdaptiveGrid PDBParser DPG UsefulMath Utility)
#SET(DAG_LINK_LIBS PDBParser DPG DynamicAdaptiveGrid)

ADD_EXECUTABLE (testDynamicAdaptiveGrid TestDynamicAdaptiveGrid.cpp)
TARGET_LINK_LIBRARIES(testDynamicAdaptiveGrid ${DAG_LINK_LIBS})

