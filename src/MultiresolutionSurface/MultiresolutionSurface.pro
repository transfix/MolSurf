TEMPLATE = lib
CONFIG  += qt warn_off staticlib create_prl
TARGET  += MultiresolutionSurface

# Input
INCLUDEPATH += ../Geometry
INCLUDEPATH += ../OpenGL_Viewer
INCLUDEPATH += ../UsefulMath
INCLUDEPATH += ../PDBParser
INCLUDEPATH += ../SmoothingFunction
INCLUDEPATH += ../AdaptiveIsocontouring
INCLUDEPATH += ../ByteOrder

SOURCES =  \
		Gaussian.cpp \
		MultiresolutionSurface.cpp \
		nfftsummer.cpp \
		Extents.cpp \
		Blurring/blurring.cpp \
		Blurring/directSum.cpp \
		Blurring/fastfft.cpp \
		Blurring/FastInverse.cpp \
		Blurring/fftAlgo.cpp


HEADERS =  \
		MultiresolutionSurface.h \
		nfftsummer.h \
		Extents.h \
		Blurring/blurring.h \
		Blurring/fastfft.h \
		Blurring/FastInverse.h \
		Blurring/fftAlgo.h


