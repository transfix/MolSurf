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
    AFX_CONTOURGEOMETRY_H__D1ABB49A_F4FB_4D74_B19A_12B037CBB7F7__INCLUDED_)
#define AFX_CONTOURGEOMETRY_H__D1ABB49A_F4FB_4D74_B19A_12B037CBB7F7__INCLUDED_

#include <Utility/utility.h>
#include <UsefulMath/Matrix.h>
#include <UsefulMath/Vector.h>
#include <Geometry/Geometry.h>

namespace CCVOpenGLMath {
class Matrix;
}
class Geometry;
struct IPolyCntl;

namespace CCVContouring {

class ContourGeometry {
public:
  ContourGeometry();
  virtual ~ContourGeometry();
  void renderContour();
  void addToGeometry(Geometry *geometry, const CCVOpenGLMath::Matrix &matrix,
                     int &nextVert, int &nextTri);
  int getNumVerts();
  int getNumTris();
  bool allocateVertexBuffers(unsigned int initialSize);
  void destroyVertexBuffers();
  bool allocateTriangleBuffers(unsigned int initialSize);
  void destroyTriangleBuffers();
  void setWireframeMode(bool state);
  void setUseColors(bool useColors);
  void setIsovalue(GLfloat isovalue);
  void setSingleColor(GLfloat R, GLfloat G, GLfloat B);
  bool useColors();
  inline bool addEdge(unsigned int &id, GLfloat v1x, GLfloat v1y, GLfloat v1z,
                      GLfloat n1x, GLfloat n1y, GLfloat n1z, GLfloat v2x,
                      GLfloat v2y, GLfloat v2z, GLfloat n2x, GLfloat n2y,
                      GLfloat n2z, GLfloat den1, GLfloat den2);
  inline bool addEdge(unsigned int &id, GLfloat v1x, GLfloat v1y, GLfloat v1z,
                      GLfloat n1x, GLfloat n1y, GLfloat n1z, GLfloat c1x,
                      GLfloat c1y, GLfloat c1z, GLfloat v2x, GLfloat v2y,
                      GLfloat v2z, GLfloat n2x, GLfloat n2y, GLfloat n2z,
                      GLfloat c2x, GLfloat c2y, GLfloat c2z, GLfloat den1,
                      GLfloat den2);
  bool addTriangle(unsigned int v1, unsigned int v2, unsigned int v3);
  inline void setVertexAndDensity(GLfloat *array, unsigned int id, GLfloat x,
                                  GLfloat y, GLfloat z, GLfloat d);
  inline void setColor(GLfloat *array, unsigned int id, GLfloat r, GLfloat g,
                       GLfloat b);
  inline void setNormal(GLfloat *array, unsigned int id, GLfloat x, GLfloat y,
                        GLfloat z);

protected:
  void setDefaults();
  // handle vertex arrays
  bool forceAllocateVertexBuffers(unsigned int size);
  bool doubleVertexBuffers();
  // handle triangle arrays
  bool forceAllocateTriangleBuffers(unsigned int size);
  bool doubleTriangleBuffers();
  bool doubleFloatArray(GLfloat *&array, unsigned int size);
  bool doubleIntArray(GLuint *&array, unsigned int size);
  void doInterpolation(GLfloat isovalue);
  inline void interpArray(GLfloat *v1, GLfloat *v2, GLfloat *n1, GLfloat *n2,
                          GLfloat d1, GLfloat d2, GLfloat isovalue);
  inline void interpArray(GLfloat *v1, GLfloat *v2, GLfloat *n1, GLfloat *n2,
                          GLfloat *c1, GLfloat *c2, GLfloat d1, GLfloat d2,
                          GLfloat isovalue);

  GLfloat *m_Normal1;
  GLfloat *m_Normal2;
  GLfloat *m_Vertex1;
  GLfloat *m_Vertex2;
  GLfloat *m_Color1;
  GLfloat *m_Color2;
  unsigned int m_NumVertices;
  unsigned int m_NumVerticesAllocated;
  GLuint *m_Triangles;
  unsigned int m_NumTriangles;
  unsigned int m_NumTrianglesAllocated;
  GLfloat m_Isovalue;
  bool m_HardwareAccelerated;
  bool m_InterpolationDone;
  bool m_UseColors;
  bool m_WireframeRender;
};
}; // namespace CCVContouring

inline bool CCVContouring::ContourGeometry::addEdge(
    unsigned int &id, GLfloat v1x, GLfloat v1y, GLfloat v1z, GLfloat n1x,
    GLfloat n1y, GLfloat n1z, GLfloat c1x, GLfloat c1y, GLfloat c1z,
    GLfloat v2x, GLfloat v2y, GLfloat v2z, GLfloat n2x, GLfloat n2y,
    GLfloat n2z, GLfloat c2x, GLfloat c2y, GLfloat c2z, GLfloat den1,
    GLfloat den2) {
  // add a single edge, doubling the array if necessary
  if (m_NumVertices == m_NumVerticesAllocated) {
    if (!doubleVertexBuffers()) {
      return false;
    }
  }
  setVertexAndDensity(m_Vertex1, m_NumVertices, v1x, v1y, v1z, den1);
  setVertexAndDensity(m_Vertex2, m_NumVertices, v2x, v2y, v2z, den2);
  setNormal(m_Normal1, m_NumVertices, n1x, n1y, n1z);
  setNormal(m_Normal2, m_NumVertices, n2x, n2y, n2z);
  setColor(m_Color1, m_NumVertices, c1x, c1y, c1z);
  setColor(m_Color2, m_NumVertices, c2x, c2y, c2z);
  id = m_NumVertices;
  m_NumVertices++;
  return true;
}

inline bool CCVContouring::ContourGeometry::addEdge(
    unsigned int &id, GLfloat v1x, GLfloat v1y, GLfloat v1z, GLfloat n1x,
    GLfloat n1y, GLfloat n1z, GLfloat v2x, GLfloat v2y, GLfloat v2z,
    GLfloat n2x, GLfloat n2y, GLfloat n2z, GLfloat den1, GLfloat den2) {
  // add a single edge, doubling the array if necessary
  if (m_NumVertices == m_NumVerticesAllocated) {
    if (!doubleVertexBuffers()) {
      return false;
    }
  }
  setVertexAndDensity(m_Vertex1, m_NumVertices, v1x, v1y, v1z, den1);
  setVertexAndDensity(m_Vertex2, m_NumVertices, v2x, v2y, v2z, den2);
  setNormal(m_Normal1, m_NumVertices, n1x, n1y, n1z);
  setNormal(m_Normal2, m_NumVertices, n2x, n2y, n2z);
  id = m_NumVertices;
  m_NumVertices++;
  return true;
}

inline void CCVContouring::ContourGeometry::setVertexAndDensity(
    GLfloat *array, unsigned int id, GLfloat x, GLfloat y, GLfloat z,
    GLfloat d) {
  // set the values in the array given x,y,z,d
  array[id * 4 + 0] = x;
  array[id * 4 + 1] = y;
  array[id * 4 + 2] = z;
  array[id * 4 + 3] = d;
}

inline void CCVContouring::ContourGeometry::setColor(GLfloat *array,
                                                     unsigned int id, GLfloat r,
                                                     GLfloat g, GLfloat b) {
  // set the values in the array given r,g,b
  array[id * 3 + 0] = r;
  array[id * 3 + 1] = g;
  array[id * 3 + 2] = b;
}

inline void CCVContouring::ContourGeometry::setNormal(GLfloat *array,
                                                      unsigned int id,
                                                      GLfloat x, GLfloat y,
                                                      GLfloat z) {
  // set the values in the array given x,y,z
  array[id * 3 + 0] = x;
  array[id * 3 + 1] = y;
  array[id * 3 + 2] = z;
}

inline void CCVContouring::ContourGeometry::interpArray(
    GLfloat *v1, GLfloat *v2, GLfloat *n1, GLfloat *n2, GLfloat d1, GLfloat d2,
    GLfloat isovalue) {
  if (d1 == d2) {
    // v1[3] = 0.0;
    //  leave v1 and n1 the way it is
  } else {
    float alpha = (isovalue - d1) / (d2 - d1);
    float oneMinusAlpha = 1.0f - alpha;
    // put the result of interpolation inside v1 and n1
    v1[0] = v1[0] * oneMinusAlpha + v2[0] * alpha;
    v1[1] = v1[1] * oneMinusAlpha + v2[1] * alpha;
    v1[2] = v1[2] * oneMinusAlpha + v2[2] * alpha;
    v1[3] = 1.0;
    n1[0] = n1[0] * oneMinusAlpha + n2[0] * alpha;
    n1[1] = n1[1] * oneMinusAlpha + n2[1] * alpha;
    n1[2] = n1[2] * oneMinusAlpha + n2[2] * alpha;
  }
}

inline void CCVContouring::ContourGeometry::interpArray(
    GLfloat *v1, GLfloat *v2, GLfloat *n1, GLfloat *n2, GLfloat *c1,
    GLfloat *c2, GLfloat d1, GLfloat d2, GLfloat isovalue) {
  if (d1 == d2) {
    // leave v1 and n1 the way it is
  } else {
    float alpha = (isovalue - d1) / (d2 - d1);
    float oneMinusAlpha = 1.0f - alpha;
    // put the result of interpolation inside v1 and n1
    v1[0] = v1[0] * oneMinusAlpha + v2[0] * alpha;
    v1[1] = v1[1] * oneMinusAlpha + v2[1] * alpha;
    v1[2] = v1[2] * oneMinusAlpha + v2[2] * alpha;
    v1[3] = 1.0;
    n1[0] = n1[0] * oneMinusAlpha + n2[0] * alpha;
    n1[1] = n1[1] * oneMinusAlpha + n2[1] * alpha;
    n1[2] = n1[2] * oneMinusAlpha + n2[2] * alpha;
    c1[0] = c1[0] * oneMinusAlpha + c2[0] * alpha;
    c1[1] = c1[1] * oneMinusAlpha + c2[1] * alpha;
    c1[2] = c1[2] * oneMinusAlpha + c2[2] * alpha;
  }
}

#endif // !defined(AFX_CONTOURGEOMETRY_H__D1ABB49A_F4FB_4D74_B19A_12B037CBB7F7__INCLUDED_)
