CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES ../../../../inc/libCG/OptimizerLib/Optimizer/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${HEADER_FILES})

#TODO: change cpp files to include with lib directory name
#e.g.: #include<UsefulMath/Ray.h> instead of #include<Ray.h>
#Then change the below to ../../inc

# also, get rid of the last include and put that somewhere better...
INCLUDE_DIRECTORIES(
  ../../../../inc
#  ../../../../inc/UsefulMath
#  ../../../../inc/PDBParser
#  ../../../../inc/libCG/OptimizerLib/Optimizer
#  /h2/wzhao/Software/levmar-2.2
)

ADD_LIBRARY (Optimizer STATIC ${SOURCE_FILES} ${INCLUDE_FILES})
