TEMPLATE = lib
CONFIG  += qt warn_off staticlib create_prl
TARGET  += GOAFileTypes

# Input
#INCLUDEPATH += ../

INCLUDEPATH += ../Geometry
INCLUDEPATH += ../UsefulMath
INCLUDEPATH += ../Blurmaps
INCLUDEPATH += ../PDBParser

SOURCES =  \
		GOAFileType.cpp \
		GOALoader.cpp \
		PDB1File.cpp \
		PDBFile.cpp \
		PQRFile.cpp \
		PTSFile.cpp \
		XYZFile.cpp \
		XYZRFile.cpp


HEADERS =  \
		GOAFileType.h \
		GOALoader.h \
		PDB1File.h \
		PDBFile.h \
		PQRFile.h \
		PTSFile.h \
		XYZFile.h \
		XYZRFile.h


