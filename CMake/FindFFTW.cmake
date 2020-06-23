## FFTW can be compiled and subsequently linked against
## various data types.
## There is a single set of include files, and then muttiple libraries,
## One for each type.  I.e. libfftw.a-->double, libfftwf.a-->float

## The following logic belongs in the individual package
## MARK_AS_ADVANCED(USE_FFTWD)
## OPTION(USE_FFTWD "Use double precision FFTW if found" ON)
## MARK_AS_ADVANCED(USE_FFTWF)
## OPTION(USE_FFTWF "Use single precision FFTW if found" ON)

  SET(FFTW_INC_SEARCHPATH
    ${PROJECT_BINARY_DIR}/include
    /sw/include
    /usr/include
    /usr/local/include
    /usr/include/fftw
    /usr/local/include/fftw
    $ENV{_NMI_PREREQ_fftw_ROOT}/include
    C:/MinGW/msys/1.0/local
    $ENV{TACC_FFTW3_INC}
  )

  FIND_PATH(FFTW_INCLUDE_PATH fftw3.h ${FFTW_INC_SEARCHPATH})

  IF(FFTW_INCLUDE_PATH)
    SET(FFTW_INCLUDE ${FFTW_INCLUDE_PATH})
  ENDIF (FFTW_INCLUDE_PATH)

 #message(FFTW_INCLUDE: ${FFTW_INCLUDE})

  #message(FFTW_INCLUDE_PATH: ${FFTW_INCLUDE_PATH})
  #IF(FFTW_INCLUDE)
  #  INCLUDE_DIRECTORIES( ${FFTW_INCLUDE})
  #ENDIF(FFTW_INCLUDE)

  GET_FILENAME_COMPONENT(FFTW_INSTALL_BASE_PATH ${FFTW_INCLUDE_PATH} PATH)

  SET(FFTW_LIB_SEARCHPATH
    ${FFTW_INSTALL_BASE_PATH}/lib
    ${PROJECT_BINARY_DIR}/lib
    /usr/lib/fftw
    /usr/local/lib/fftw
    C:/MinGW/msys/1.0/local
    $ENV{TACC_FFTW3_INC}
  )
 
  MARK_AS_ADVANCED(FFTW_LIB)
  FIND_LIBRARY(FFTW_LIB libfftw3.a ${FFTW_LIB_SEARCHPATH})
  if(FFTW_LIB)
     set(FFTW_FOUND FOUND)
  endif(FFTW_LIB)

  MARK_AS_ADVANCED(FFTWD_LIB FFTWD_THREADS_LIB)
  FIND_LIBRARY(FFTWD_LIB libfftw3d.a ${FFTW_LIB_SEARCHPATH}) #Double Precision Lib
  FIND_LIBRARY(FFTWD_THREADS_LIB fftw3d_threads.a ${FFTW_LIB_SEARCHPATH}) #Double Precision Lib only if compiled with threads support

  IF(FFTWD_LIB)
     SET(FFTW_LIB ${FFTWD_LIB} ${FFTW_LIB})
  ENDIF(FFTWD_LIB)
  IF(FFTWD_THREADS_LIB)
     SET(FFTW_LIB ${FFTWD_THREADS_LIB} ${FFTW_LIB})
  ENDIF(FFTWD_THREADS_LIB)

  MARK_AS_ADVANCED(FFTWF_LIB FFTWF_THREADS_LIB)
  FIND_LIBRARY(FFTWF_LIB libfftw3f.a ${FFTW_LIB_SEARCHPATH}) #Single Precision Lib
  FIND_LIBRARY(FFTWF_THREADS_LIB libfftw3f_threads.a ${FFTW_LIB_SEARCHPATH}) #Single Precision Lib only if compiled with threads support

  IF(FFTWF_LIB)
     SET(FFTW_LIB ${FFTWF_LIB} ${FFTW_LIB})
  ENDIF(FFTWF_LIB)
  IF(FFTWF_THREADS_LIB)
     SET(FFTW_LIB ${FFTWF_THREADS_LIB} ${FFTW_LIB})
  ENDIF(FFTWF_THREADS_LIB)

  #message(FFTW_LIB: ${FFTW_LIB})
  #message(FFTW_INCLUDE: ${FFTW_INCLUDE})
