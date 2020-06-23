TEMPLATE = lib
CONFIG  += qt warn_off staticlib create_prl
TARGET  += Optimizer

INCLUDEPATH += ${TEXMOL}/GeometryRenderer
INCLUDEPATH += ${TEXMOL}/GeometryFileTypes
INCLUDEPATH += ${TEXMOL}/Geometry
INCLUDEPATH += ${TEXMOL}/ColorManager
INCLUDEPATH += ${TEXMOL}/Contouring
INCLUDEPATH += ${TEXMOL}/PDBParser
INCLUDEPATH += ${TEXMOL}/UsefulMath
INCLUDEPATH += ${TEXMOL}/Blurmaps
INCLUDEPATH += ${TEXMOL}/contour
INCLUDEPATH += ${TEXMOL}/DataManager
INCLUDEPATH += ${TEXMOL}/HLevelSet
INCLUDEPATH += ${TEXMOL}/ImposterRenderer
INCLUDEPATH += ${TEXMOL}/MolecularCharacteristics
INCLUDEPATH += ${TEXMOL}/ColorTable
INCLUDEPATH += ${TEXMOL}/VolumeLibrary
INCLUDEPATH += ${TEXMOL}/VolumeFileTypes
INCLUDEPATH += ${TEXMOL}/GOAFileTypes
INCLUDEPATH += ${TEXMOL}/RayTracer
INCLUDEPATH += ${TEXMOL}/Nurbs
INCLUDEPATH += ${TEXMOL}/SignDistanceFunction
INCLUDEPATH += ${TEXMOL}/SignDistanceFunction_v2
INCLUDEPATH += ${TEXMOL}/Summation
INCLUDEPATH += ${TEXMOL}/Docking
INCLUDEPATH += ${TEXMOL}/OpenGL_Viewer
INCLUDEPATH += ${TEXMOL}/Dialogs
INCLUDEPATH += ${TEXMOL}/Pocket
INCLUDEPATH += ${TEXMOL}/PatchComplex
#INCLUDEPATH += /home/albertc/BEM/nfft/include/nfft
#INCLUDEPATH += /home/albertc/BEM/fftw/include

# Input

# gmp

INCLUDEPATH += /usr/include \
               /h2/wzhao/Software/levmar-2.2

HEADERS += optimizer.h \ 
	   optimizationFunction.h \ 
	   geometry.h \ 
	   surface.h \ 

SOURCES += optimizer.cpp \ 
	   optimizationFunction.cpp \ 
	   surface.cpp \ 


