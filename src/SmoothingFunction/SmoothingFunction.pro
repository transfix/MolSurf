TEMPLATE = lib
CONFIG  += qt warn_off staticlib create_prl
TARGET  += SmoothingFunction

# Input

SOURCES =  \
		CardinalBSpline.cpp \
		Gaussian.cpp \
		SmoothingFunction.cpp 


HEADERS =  \
		CardinalBSpline.h \
		Gaussian.h \
		SmoothingFunction.h 

