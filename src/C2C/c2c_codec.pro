TEMPLATE = lib
CONFIG  += warn_off staticlib create_prl
TARGET  += c2c_codec
INCLUDEPATH += . ./arithlib ./arithlib/libdjvu++
DEFINES += LINUX _LITTLE_ENDIAN ARITH_ENCODE ARITH_DECODE UNIQUE ZP_CODEC \
		#_DEBUG

# Input

SOURCES =  \
		c2c_codec.cpp \
		bufferedio.cpp \
		ContourGeom.cpp \
		util.cpp \
		rawslicefac.cpp \
		rawvslicefac.cpp \
		slicecache.cpp \
		slice.cpp \
		layer.cpp \
		vertex.cpp \
		arithlib/arithdecode.cpp \
		arithlib/bitbuffer.cpp \
		arithlib/coder.cpp \
		arithlib/encode.cpp \
		arithlib/stats.cpp \
		arithlib/utils.cpp \
		arithlib/libdjvu++/ByteStream.cpp \
		arithlib/libdjvu++/GContainer.cpp \
		arithlib/libdjvu++/GException.cpp \
		arithlib/libdjvu++/GSmartPointer.cpp \
		arithlib/libdjvu++/ZPCodec.cpp
		

HEADERS =  \
		c2c_codec.h \
		bufferedio.h \
		c2cbuf.h \
		ContourGeom.h \
		cubes.h \
		decode.h \
		diskio.h \
		filec2cbuf.h \
		util.h \
		vtkMarchingCubesCases.h \
		arithlib/arith_defines.h \
		arithlib/arith.h \
		arithlib/bitbuffer.h \
		arithlib/coder.h \
		arithlib/encode.h \
		arithlib/stats.h \
		arithlib/unitypes.h \
		arithlib/utils.h \
		arithlib/libdjvu++/ByteStream.h \
		arithlib/libdjvu++/GContainer.h \
		arithlib/libdjvu++/GException.h \
		arithlib/libdjvu++/GSmartPointer.h \
		arithlib/libdjvu++/ZPCodec.h

