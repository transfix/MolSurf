CMAKE_MINIMUM_REQUIRED (VERSION 2.8)

include(ExternalProject)

#------------------------------------------------------------
# gmp
#------------------------------------------------------------
ExternalProject_Add(gmp_proj
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/third-party/gmp
    CONFIGURE_COMMAND ${PROJECT_SOURCE_DIR}/third-party/gmp/configure --prefix=${PROJECT_BINARY_DIR} --disable-shared
    BUILD_COMMAND ${MAKE})

#------------------------------------------------------------
# mpfr
#------------------------------------------------------------
ExternalProject_Add(mpfr_proj
    DEPENDS gmp_proj
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/third-party/mpfr
    CONFIGURE_COMMAND ${PROJECT_SOURCE_DIR}/third-party/mpfr/configure --prefix=${PROJECT_BINARY_DIR} --disable-shared --with-gmp=${PROJECT_BINARY_DIR}
    BUILD_COMMAND ${make})

#------------------------------------------------------------
# boost
#------------------------------------------------------------
ExternalProject_Add(boost_proj
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/third-party/boost
    CMAKE_ARGS -DBOOST_INCLUDE_INSTALL_DIR=include
      -DBUILD_SOVERSIONED=OFF -DBOOST_LIB_INSTALL_DIR=lib -DENABLE_DEBUG=ON
      -DENABLE_MULTI_THREADED=ON -DENABLE_RELEASE=ON -DENABLE_SHARED=OFF
      -DENABLE_SINGLE_THREADED=OFF -DENABLE_STATIC=ON -DINSTALL_VERSIONED=OFF
      -DWITH_MPI=OFF -DWITH_PYTHON=OFF
      -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
    BUILD_COMMAND ${MAKE}
    INSTALL_DIR ${PROJECT_BINARY_DIR})

#------------------------------------------------------------
# CGAL
#------------------------------------------------------------
ExternalProject_Add(CGAL_proj
    DEPENDS boost_proj mpfr_proj gmp_proj
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/third-party/cgal
    CMAKE_ARGS -DBUILD_SHARED_LIBS=FALSE
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
    BUILD_COMMAND ${MAKE}
    INSTALL_DIR ${PROJECT_BINARY_DIR})


#------------------------------------------------------------
# lapack
#------------------------------------------------------------
ExternalProject_Add(lapack_proj
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/third-party/lapack
    CMAKE_ARGS -DBUILD_SHARED_LIBS=FALSE
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
    BUILD_COMMAND ${MAKE}
    INSTALL_DIR ${PROJECT_BINARY_DIR})
#------------------------------------------------------------
# fftw
#------------------------------------------------------------
ExternalProject_Add(fftw_proj
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/third-party/fftw
    CONFIGURE_COMMAND ${PROJECT_SOURCE_DIR}/third-party/fftw/configure
--prefix=${PROJECT_BINARY_DIR} --enable-shared=no --enable-threads --enable-static=yes
    BUILD_COMMAND ${make})

#------------------------------------------------------------
# nfft
#------------------------------------------------------------
ExternalProject_Add(nfft_proj
    DEPENDS fftw_proj
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/third-party/nfft
    CONFIGURE_COMMAND ${PROJECT_SOURCE_DIR}/third-party/nfft/configure
--prefix=${PROJECT_BINARY_DIR} --with-fftw3=${PROJECT_BINARY_DIR}
--enable-shared=no --enable-static=yes
    BUILD_COMMAND ${make})

