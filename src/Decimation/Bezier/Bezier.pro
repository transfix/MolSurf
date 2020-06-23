TEMPLATE = lib
CONFIG  += qt warn_off staticlib create_prl
TARGET  += Bezier

INCLUDEPATH += ../include
INCLUDEPATH += ../Control
INCLUDEPATH += ../Decim
INCLUDEPATH += ../Linear
INCLUDEPATH += ../Misc

# Input
HEADERS = \
		BezierMath.h

SOURCES = \
		BezierMath.c
