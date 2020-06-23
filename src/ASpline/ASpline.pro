TEMPLATE = lib
CONFIG  += qt warn_off staticlib create_prl
TARGET   = ASpline 

SOURCES = \ 
	  aspline.cpp \ 
	  bezier.cpp \ 
	  geometryParser.cpp \ 
	  quadrature.cpp

HEADERS = \  
	  aspline.h \ 
	  geometry.h \ 
	  geometryParser.h \ 
	  quadrature.h

