#
# This macro is for setting up a sub-project to use FFTW (float or double version)
#

include(SetupFFTW)

macro(SetupNFFT TargetName)

  set(NFFT_FOUND FOUND)
  find_package(NFFT)
  #MESSAGE(NFFT_LIB: ${NFFT_LIB})

  if(NFFT_LIB)
  target_link_libraries(${TargetName} ${NFFT_LIB})
  endif(NFFT_LIB)

  IF(NFFT_INCLUDE)
    INCLUDE_DIRECTORIES( ${NFFT_INCLUDE})
  ENDIF(NFFT_INCLUDE)

  SetupFFTW(${TargetName})

endmacro(SetupNFFT)
