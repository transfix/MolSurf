TEMPLATE = lib
CONFIG  += qt warn_off staticlib create_prl
TARGET  += MolecularCharacteristics

# Input
INCLUDEPATH += ../Blurmaps
INCLUDEPATH += ../OpenGL_Viewer
INCLUDEPATH += ../PDBParser
INCLUDEPATH += ../UsefulMath
INCLUDEPATH += ../VolumeFileTypes
INCLUDEPATH += ../Geometry
INCLUDEPATH += ../C2C
INCLUDEPATH += ../

SOURCES =  \
		MolecularCharacteristics.cpp


HEADERS =  \
		MolecularCharacteristics.h 


