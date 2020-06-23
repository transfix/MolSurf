TEMPLATE = lib
CONFIG  += qt warn_off staticlib create_prl
TARGET  += DecimXu

INCLUDEPATH += ../include
INCLUDEPATH += ../Decim
INCLUDEPATH += ../Linear
INCLUDEPATH += ../Misc

# Input
HEADERS = \
		fatnurb.h \
		mesh.h \
		xudecimate.h 

SOURCES = \
		cuttri.c \
		xudecimate.c
