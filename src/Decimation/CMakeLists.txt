CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

FILE(GLOB SOURCE_FILES
       	*.cpp
	#Misc/*.c
	Misc/dscal.c
	Misc/lineralge.c
	Misc/testsvd.c
	Misc/sub1.c
	Misc/dsvdc.c
	Linear/*.c
	Bezier/*.c
	DecimXu/*.c
	Normal/*.c
	Decim/*.c
)
FILE(GLOB INCLUDE_FILES 
	../../inc/Decimation/*.[Hh]
	../../inc/Decimation/Misc/*.[Hh]
	../../inc/Decimation/Linear/*.[Hh]
	../../inc/Decimation/Bezier/*.[Hh]
	../../inc/Decimation/DecimXu/*.[Hh]
	../../inc/Decimation/Normal/*.[Hh]
	../../inc/Decimation/Decim/*.[Hh]
)

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${HEADER_FILES})

INCLUDE_DIRECTORIES(
  ../../inc
)

ADD_LIBRARY (Decimation STATIC ${SOURCE_FILES} ${INCLUDE_FILES})
