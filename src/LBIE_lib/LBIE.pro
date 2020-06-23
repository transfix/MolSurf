TEMPLATE = app
CONFIG  += qt warn_off link_prl opengl
TARGET  += LBIE_test

LIBS += \
  ${TEXMOL}/Blurmaps/libBlurmaps.a \
  ${TEXMOL}/ByteOrder/libByteOrder.a \
  ${TEXMOL}/ColorTable/libColorTable.a \
  ${TEXMOL}/Contouring/libContouring.a \
  ${TEXMOL}/Curvatures/libCurvatures.a \
  ${TEXMOL}/DataManager/libDataManager.a \
  ${TEXMOL}/C2C/libC2C.a \
  ${TEXMOL}/ColorManager/libColorManager.a \
  ${TEXMOL}/contour/libcontour.a \
  ${TEXMOL}/Docking/libDocking.a \
  ${TEXMOL}/GeometryFileTypes/libGeometryFileTypes.a \
  ${TEXMOL}/GeometryRenderer/libGeometryRenderer.a \
  ${TEXMOL}/Geometry/libGeometry.a \
  ${TEXMOL}/GOAFileTypes/libGOAFileTypes.a \
  ${TEXMOL}/HLevelSet/libHLevelSet.a \
  ${TEXMOL}/ImposterRenderer/libImposterRenderer.a \
  ${TEXMOL}/MolecularCharacteristics/libMolecularCharacteristics.a \
  ${TEXMOL}/MultiresolutionSurface/libMultiresolutionSurface.a \
  ${TEXMOL}/NMA/libNMA.a \
  ${TEXMOL}/Nurbs/libNurbs.a \
  ${TEXMOL}/ObjectRenderer/libObjectRenderer.a \
  ${TEXMOL}/OpenGL_Viewer/libOpenGL_Viewer.a \
  ${TEXMOL}/PatchComplex/libPatchComplex.a \
  ${TEXMOL}/PDBParser/libPDBParser.a \
  ${TEXMOL}/Pocket/libPocket.a \
  ${TEXMOL}/RayTracer/libRayTracer.a \
  ${TEXMOL}/SphericalPatchIntersections/libSphericalPatchIntersections.a \
  ${TEXMOL}/SignDistanceFunction/libSignDistanceFunction.a \
  ${TEXMOL}/SignDistanceFunction_v2/libSignDistanceFunction_v2.a \
  ${TEXMOL}/SmoothingFunction/libSmoothingFunction.a \
  ${TEXMOL}/SplattingRenderer/libSplattingRenderer.a \
  ${TEXMOL}/Summation/libSummation.a \
  ${TEXMOL}/UsefulMath/libUsefulMath.a \
  ${TEXMOL}/VolumeFileTypes/libVolumeFileTypes.a \
  ${TEXMOL}/VolumeLibrary/libVolumeLibrary.a \
  ${TEXMOL}/XMLParser/libXMLParser.a \
  ${TEXMOL}/ColorManager/libColorManager.a \
  -lstdc++



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
INCLUDEPATH += /home/albertc/BEM/nfft/include/nfft
INCLUDEPATH += /home/albertc/BEM/fftw/include
INCLUDEPATH += ./

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
                pcio.h \
		SimpleVolumeDataIsocontourer.h

SOURCES = \
		LBIE_lib_test.cpp \
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
		tetra.cpp \
		SimpleVolumeDataIsocontourer.cpp
