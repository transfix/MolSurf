CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES ../../inc/MolecularCharacteristics/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${HEADER_FILES})

INCLUDE_DIRECTORIES(
  ../../inc
  ../../inc/DataManager # this is a hack for TexMol.... resolve someday
)

ADD_LIBRARY (MolecularCharacteristics STATIC ${SOURCE_FILES} ${INCLUDE_FILES})
