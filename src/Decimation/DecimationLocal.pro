TEMPLATE = lib
CONFIG  += qt warn_off staticlib create_prl
TARGET  += Decimation

INCLUDEPATH += ./include
INCLUDEPATH += ./DecimXu
INCLUDEPATH += ./Decim
INCLUDEPATH += ./Linear
INCLUDEPATH += ./Misc
INCLUDEPATH += ./Bezier
INCLUDEPATH += ./Normal
INCLUDEPATH += ../LBIE_lib
INCLUDEPATH += ../

# Input
HEADERS = \
		Decimation.h

SOURCES = \
		Decimation.cpp
