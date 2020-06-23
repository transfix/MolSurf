TEMPLATE = lib
CONFIG  += qt warn_off staticlib create_prl
TARGET  += Normal

INCLUDEPATH += ../include
INCLUDEPATH += ../Control
INCLUDEPATH += ../Decim
INCLUDEPATH += ../Linear
INCLUDEPATH += ../Misc

# Input
HEADERS = \
		normal.h
		
SOURCES = \
                trinormal.c
