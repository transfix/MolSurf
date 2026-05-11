# SetupFFTW.cmake
#
# Attach the system FFTW3 to a target.

if(NOT MolSurf_FFTW_Found)
  find_package(FFTW REQUIRED)
  set(MolSurf_FFTW_Found TRUE CACHE INTERNAL "")
endif()

macro(SetupFFTW TargetName)
  if(TARGET FFTW3::fftw3)
    target_link_libraries(${TargetName} PUBLIC FFTW3::fftw3)
  endif()
  if(TARGET FFTW3::fftw3_threads)
    target_link_libraries(${TargetName} PUBLIC FFTW3::fftw3_threads)
  endif()
endmacro()
