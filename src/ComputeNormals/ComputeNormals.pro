TEMPLATE = lib
CONFIG += qt warn_off staticlib create_prl
TARGET = Normal

SOURCES = \
		  geometryParser.cpp \  
		  normalByAverage.cpp \ 
		  normalByLoop.cpp \ 
		  produceFaceNormal.cpp \ 
		  produceNormal.cpp \
		  surface.cpp \ 
		  valenceNeighbor.cpp\ 

HEADERS = \ 
		  geometry.h \ 
		  geometryParser.h \ 
		  surface.h \ 

