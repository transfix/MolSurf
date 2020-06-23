TEMPLATE = lib
CONFIG  += qt warn_off staticlib create_prl
TARGET  += Decim

INCLUDEPATH += ../include
INCLUDEPATH += ../Decim
INCLUDEPATH += ../Linear
INCLUDEPATH += ../Misc

# Input
HEADERS = \
		ct.h  \
		mesh.h  \
		sdecim.h  \
		sdecimPriv.h

SOURCES = \
		sdec.c  \
		sdecim.c  \
		sedge.c  \
		ssimple.c
