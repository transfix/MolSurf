TEMPLATE = lib
CONFIG   = qt warn_off staticlib create_prl
TARGET   = SignDistanceFunction

INCLUDEPATH  = ../Blurmaps
INCLUDEPATH += ../Geometry


SOURCES =  \
		common.cpp \
		compute.cpp \
		head.cpp \
		init.cpp \
		main.cpp \
		new_adjust.cpp \
		octree.cpp \
		propagate.cpp


HEADERS =  \
		common.h \
		head.h \
		sdfLib.h 


