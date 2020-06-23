TEMPLATE = lib
CONFIG  += qt warn_off staticlib create_prl
TARGET  += LBIE_lib

INCLUDEPATH += ../
INCLUDEPATH += ../GeometryFileTypes
INCLUDEPATH += ../Geometry
INCLUDEPATH += ../PDBParser
INCLUDEPATH += ../Blurmaps
INCLUDEPATH += ../MolecularCharacteristics
INCLUDEPATH += ../VolumeLibrary
INCLUDEPATH += ../VolumeFileTypes
INCLUDEPATH += ../GOAFileTypes
INCLUDEPATH += ../SignDistanceFunction_v2
INCLUDEPATH += ../Summation
INCLUDEPATH += ../Pocket
INCLUDEPATH += ../PatchComplex
INCLUDEPATH += ../UsefulMath

# Input
HEADERS = \
		LBIE_Mesher.h \
                cellQueue.h \
                contour3d.h \
                e_face.h \
                geoframe.h \
                invertMatrix.h \
                normalspline.h \
                octree.h \
                pcio.h 

SOURCES = \
		LBIE_Mesher.cpp \
		cellQueue.cpp \
		contour3d.cpp \
		e_face.cpp \
		geoframe.cpp \
		hexa.cpp \
		interp_bit_generate.cpp \
		invertMatrix.cpp \
		normalspline.cpp \
		octree.cpp \
		pcio.cpp \
		tetra.cpp
