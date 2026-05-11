/*
  Copyright 2011 The University of Texas at Austin

        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of MolSurf.

  MolSurf is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.


  MolSurf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with MolSurf; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#if !defined(                                                                  \
    AFX_SCENEARRAYNODE_H__1F8CDAB9_2D09_475B_9E34_A6139FF29DDC__INCLUDED_)
#define AFX_SCENEARRAYNODE_H__1F8CDAB9_2D09_475B_9E34_A6139FF29DDC__INCLUDED_

#include <Geometry/Geometry.h>
#include <UsefulMath/Matrix.h>

#ifdef _WIN32
#include <windows.h>
#endif

#if defined(__APPLE__)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

class SceneArrayNode {
public:
  SceneArrayNode(GLuint initializationCallList = 0,
                 GLuint finalizationCallList = 0);
  SceneArrayNode(Geometry *geometry, GLuint initializationCallList = 0,
                 GLuint finalizationCallList = 0);
  virtual ~SceneArrayNode();
  void translate(float tx, float ty, float tz);
  void rotation(float angle, float x, float y, float z);
  void scale(float s);
  void setGeometry(Geometry *geometry);
  Geometry *getGeometry();
  const Geometry *getGeometry() const;
  void setTransformation(const CCVOpenGLMath::Matrix &transformation);
  CCVOpenGLMath::Matrix &getTransformation();
  const CCVOpenGLMath::Matrix &getTransformation() const;
  void setInitializationCallList(GLuint initializationCallList);
  GLuint getInitializationCallList() const;
  void setFinalizationCallList(GLuint finalizationCallList);
  GLuint getFinalizationCallList() const;

protected:
  Geometry *m_Geometry;
  GLuint m_InitializationCallList;
  GLuint m_FinalizationCallList;
  CCVOpenGLMath::Matrix m_Transformation;
};

#endif
