TEMPLATE = lib
CONFIG  += qt warn_off staticlib create_prl
TARGET  += Linear

INCLUDEPATH += ../include
INCLUDEPATH += ../Control
INCLUDEPATH += ../Decim
INCLUDEPATH += ../Linear
INCLUDEPATH += ../Misc

# Input
HEADERS = \
		defs.h  \
		display.h  \
		sturm.h

SOURCES = \
		display.c  \
		root.c  \
		sturm.c
