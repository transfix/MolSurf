CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES ../../inc/libCG/CoarseGrain/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${HEADER_FILES})

# arand: this is not building


#TODO: change cpp files to include with lib directory name
#e.g.: #include<UsefulMath/Ray.h> instead of #include<Ray.h>
#Then change the below to ../../inc
INCLUDE_DIRECTORIES(
  ../../../inc/libCG/CoarseGrain
  ../../../inc/
  #/h2/wzhao/Software/
#  /h2/wzhao/Software/nfft-3.0.2/include/nfft
)

ADD_LIBRARY (CoarseGrain STATIC ${SOURCE_FILES} ${INCLUDE_FILES})

SetupNFFT(CoarseGrain)
