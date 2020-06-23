#
# This macro is for setting up a sub-project to use GLUT
#

macro(SetupGLUT TargetName)

  set(GLUT_FOUND NOTFOUND)

  find_package(GLUT)
  find_package(OpenGL)

  MESSAGE(GLUT_LIB: ${GLUT_LIB})
  MESSAGE(GLUT_INCLUDE: ${GLUT_INCLUDE_DIR})
  MESSAGE(OPENGL_LIB: ${OPENGL_LIB})
  MESSAGE(OPENGL_INCLUDE: ${OPENGL_INCLUDE_DIR})

  include_directories(${GLUT_INCLUDE_DIR})
  include_directories(${OPENGL_INCLUDE_DIR})

  target_link_libraries(${TargetName} ${GLUT_LIB} ${OPENGL_LIB})

endmacro(SetupGLUT)
