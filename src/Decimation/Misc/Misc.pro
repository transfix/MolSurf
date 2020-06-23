TEMPLATE = lib
CONFIG  += qt warn_off staticlib create_prl
TARGET  += Misc

INCLUDEPATH += ../include
INCLUDEPATH += ../Decim
INCLUDEPATH += ../Linear
INCLUDEPATH += ../Misc
INCLUDEPATH += ../Bezier

# Input
HEADERS = \
		misc.h

SOURCES = \
		dscal.c  \
		dsvdc.c  \
		lineralge.c  \
		sub1.c  \
		testsvd.c
