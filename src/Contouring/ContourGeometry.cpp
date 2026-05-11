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
#include <Contouring/ContourGeometry.h>

using CCVContouring::ContourGeometry;

ContourGeometry::ContourGeometry() { setDefaults(); }

ContourGeometry::~ContourGeometry() {
  destroyTriangleBuffers();
  destroyVertexBuffers();
}

void ContourGeometry::setWireframeMode(bool state) {
  m_WireframeRender = state;
}

void ContourGeometry::renderContour() {
  if (!m_Triangles || !m_Vertex1 || !m_Normal1) {
    return;
  }
  doInterpolation(m_Isovalue);
  glPushAttrib(GL_LIGHTING_BIT | GL_ENABLE_BIT);
  glEnable(GL_LIGHTING);
  glEnable(GL_NORMALIZE);
  if (m_WireframeRender) {
    unsigned int i;
    if (m_UseColors) // draw with color
    {
      glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
      glEnable(GL_COLOR_MATERIAL);
      for (i = 0; i < m_NumTriangles; i++) {
        glBegin(GL_LINE_LOOP);
        glColor3fv(&(m_Color1[m_Triangles[i * 3 + 0] * 3]));
        glNormal3fv(&(m_Normal1[m_Triangles[i * 3 + 0] * 3]));
        glVertex4fv(&(m_Vertex1[m_Triangles[i * 3 + 0] * 4]));
        glColor3fv(&(m_Color1[m_Triangles[i * 3 + 1] * 3]));
        glNormal3fv(&(m_Normal1[m_Triangles[i * 3 + 1] * 3]));
        glVertex4fv(&(m_Vertex1[m_Triangles[i * 3 + 1] * 4]));
        glColor3fv(&(m_Color1[m_Triangles[i * 3 + 2] * 3]));
        glNormal3fv(&(m_Normal1[m_Triangles[i * 3 + 2] * 3]));
        glVertex4fv(&(m_Vertex1[m_Triangles[i * 3 + 2] * 4]));
        glEnd();
      }
    } else // draw without color
    {
      for (i = 0; i < m_NumTriangles; i++) {
        glBegin(GL_LINE_LOOP);
        glNormal3fv(&(m_Normal1[m_Triangles[i * 3 + 0] * 3]));
        glVertex4fv(&(m_Vertex1[m_Triangles[i * 3 + 0] * 4]));
        glNormal3fv(&(m_Normal1[m_Triangles[i * 3 + 1] * 3]));
        glVertex4fv(&(m_Vertex1[m_Triangles[i * 3 + 1] * 4]));
        glNormal3fv(&(m_Normal1[m_Triangles[i * 3 + 2] * 3]));
        glVertex4fv(&(m_Vertex1[m_Triangles[i * 3 + 2] * 4]));
        glEnd();
      }
    }
  } else // triangle render
  {
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    if (m_UseColors) {
      glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
      glEnable(GL_COLOR_MATERIAL);
      glEnableClientState(GL_COLOR_ARRAY);
      glColorPointer(3, GL_FLOAT, 0, m_Color1);
    }
    glVertexPointer(4, GL_FLOAT, 0, m_Vertex1);
    glNormalPointer(GL_FLOAT, 0, m_Normal1);
    glDrawElements(GL_TRIANGLES, m_NumTriangles * 3, GL_UNSIGNED_INT,
                   m_Triangles);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
  }
  glPopAttrib();
}

void ContourGeometry::addToGeometry(Geometry *geometry,
                                    const CCVOpenGLMath::Matrix &matrix,
                                    int &nextVert, int &nextTri) {
  if (!m_Triangles || !m_Vertex1 || !m_Normal1) {
    return;
  }
  int nextVertPrivate = nextVert;
  CCVOpenGLMath::Matrix inverseTranspose = matrix.inverseTranspose();
  doInterpolation(m_Isovalue);
  unsigned int c;
  CCVOpenGLMath::Vector vector;
  CCVOpenGLMath::Vector point, normal;
  for (c = 0; c < m_NumVertices; c++) {
    point[0] = m_Vertex1[c * 4 + 0];
    point[1] = m_Vertex1[c * 4 + 1];
    point[2] = m_Vertex1[c * 4 + 2];
    point[3] = 1.0;
    normal[0] = m_Normal1[c * 3 + 0];
    normal[1] = m_Normal1[c * 3 + 1];
    normal[2] = m_Normal1[c * 3 + 2];
    normal[3] = 0.0;
    // set the data up for the right coord space
    point = matrix * point;
    normal = inverseTranspose * normal;
    geometry->m_TriVerts[nextVertPrivate * 3 + 0] = point[0];
    geometry->m_TriVerts[nextVertPrivate * 3 + 1] = point[1];
    geometry->m_TriVerts[nextVertPrivate * 3 + 2] = point[2];
    geometry->m_TriVertNormals[nextVertPrivate * 3 + 0] = normal[0];
    geometry->m_TriVertNormals[nextVertPrivate * 3 + 1] = normal[1];
    geometry->m_TriVertNormals[nextVertPrivate * 3 + 2] = normal[2];
    if (useColors()) {
      geometry->m_TriVertColorsTransparent[nextVertPrivate * 4 + 0] =
          m_Color1[c * 3 + 0];
      geometry->m_TriVertColorsTransparent[nextVertPrivate * 4 + 1] =
          m_Color1[c * 3 + 1];
      geometry->m_TriVertColorsTransparent[nextVertPrivate * 4 + 2] =
          m_Color1[c * 3 + 2];
    }
    nextVertPrivate++;
  }
  for (c = 0; c < m_NumTriangles; c++) {
    geometry->m_Tris[nextTri * 3 + 0] = m_Triangles[c * 3 + 0] + nextVert;
    geometry->m_Tris[nextTri * 3 + 1] = m_Triangles[c * 3 + 1] + nextVert;
    geometry->m_Tris[nextTri * 3 + 2] = m_Triangles[c * 3 + 2] + nextVert;
    nextTri++;
  }
  nextVert = nextVertPrivate;
}

int ContourGeometry::getNumVerts() { return m_NumVertices; }

int ContourGeometry::getNumTris() { return m_NumTriangles; }

void ContourGeometry::setDefaults() {
  // set default values for all member variables
  m_Normal1 = 0;
  m_Normal2 = 0;
  m_Vertex1 = 0;
  m_Vertex2 = 0;
  m_Color1 = 0;
  m_Color2 = 0;
  m_NumVertices = 0;
  m_NumVerticesAllocated = 0;
  m_Triangles = 0;
  m_NumTriangles = 0;
  m_NumTrianglesAllocated = 0;
  m_Isovalue = 0.0f;
  m_HardwareAccelerated = false;
  m_InterpolationDone = false;
  m_UseColors = false;
  m_WireframeRender = false;
}

void ContourGeometry::setUseColors(bool useColors) { m_UseColors = useColors; }

void ContourGeometry::setIsovalue(GLfloat isovalue) { m_Isovalue = isovalue; }

void ContourGeometry::setSingleColor(GLfloat R, GLfloat G, GLfloat B) {
  unsigned int i;
  // reassign all the vertex colors
  for (i = 0; i < m_NumVertices; i++) {
    m_Color1[i * 3 + 0] = R;
    m_Color1[i * 3 + 1] = G;
    m_Color1[i * 3 + 2] = B;
  }
}

bool ContourGeometry::useColors() { return m_UseColors; }

bool ContourGeometry::allocateVertexBuffers(unsigned int initialSize) {
  // if we already have a buffer big enough, just return true and use that
  // buffer
  if (initialSize > m_NumVerticesAllocated) {
    destroyVertexBuffers();
    m_NumVertices = 0;
    m_InterpolationDone = false;
    return forceAllocateVertexBuffers(initialSize);
  } else {
    m_NumVertices = 0;
    m_InterpolationDone = false;
    return true;
  }
}

bool ContourGeometry::forceAllocateVertexBuffers(unsigned int size) {
  // allocate the buffers without checking to see if we already have a big
  // enough buffer
  m_Vertex1 = new GLfloat[size * 4];
  m_Vertex2 = new GLfloat[size * 4];
  m_Normal1 = new GLfloat[size * 3];
  m_Normal2 = new GLfloat[size * 3];
  m_Color1 = new GLfloat[size * 3];
  m_Color2 = new GLfloat[size * 3];
  if (m_Vertex1 && m_Vertex2 && m_Normal1 && m_Normal2 && m_Color1 &&
      m_Color2) {
    m_NumVerticesAllocated = size;
    return true;
  } else {
    destroyVertexBuffers();
    return false;
  }
}

bool ContourGeometry::doubleVertexBuffers() {
  // double the size of the vertex buffer and copy the vertices to the new
  // buffer
  if (!doubleFloatArray(m_Vertex1, m_NumVerticesAllocated * 4) ||
      !doubleFloatArray(m_Vertex2, m_NumVerticesAllocated * 4) ||
      !doubleFloatArray(m_Normal1, m_NumVerticesAllocated * 3) ||
      !doubleFloatArray(m_Normal2, m_NumVerticesAllocated * 3) ||
      !doubleFloatArray(m_Color1, m_NumVerticesAllocated * 3) ||
      !doubleFloatArray(m_Color2, m_NumVerticesAllocated * 3)) {
    destroyVertexBuffers();
    return false;
  } else // success
  {
    m_NumVerticesAllocated *= 2;
    return true;
  }
}

void ContourGeometry::destroyVertexBuffers() {
  // free vertex buffer memory
  delete[] m_Vertex1;
  m_Vertex1 = 0;
  delete[] m_Vertex2;
  m_Vertex2 = 0;
  delete[] m_Normal1;
  m_Normal1 = 0;
  delete[] m_Normal2;
  m_Normal2 = 0;
  delete[] m_Color1;
  m_Color1 = 0;
  delete[] m_Color2;
  m_Color2 = 0;
  m_NumVerticesAllocated = 0;
  m_NumVertices = 0;
}

bool ContourGeometry::allocateTriangleBuffers(unsigned int initialSize) {
  // if we already have a buffer big enough, just return true and use that
  // buffer
  if (initialSize > m_NumTrianglesAllocated) {
    destroyTriangleBuffers();
    m_NumTriangles = 0;
    m_InterpolationDone = false;
    return forceAllocateTriangleBuffers(initialSize);
  } else {
    m_NumTriangles = 0;
    m_InterpolationDone = false;
    return true;
  }
}

bool ContourGeometry::forceAllocateTriangleBuffers(unsigned int size) {
  // allocate the buffers without checking to see if we already have a big
  // enough buffer
  m_Triangles = new GLuint[size * 3];
  if (m_Triangles) {
    m_NumTrianglesAllocated = size;
    return true;
  } else {
    destroyTriangleBuffers();
    return false;
  }
}

bool ContourGeometry::doubleTriangleBuffers() {
  // double the size of the triangle buffer and copy the triangles to the new
  // buffer
  if (!doubleIntArray(m_Triangles, m_NumTrianglesAllocated * 3)) {
    destroyTriangleBuffers();
    return false;
  } else // success
  {
    m_NumTrianglesAllocated *= 2;
    return true;
  }
}

void ContourGeometry::destroyTriangleBuffers() {
  // free triangle buffer memory
  delete[] m_Triangles;
  m_Triangles = 0;
  m_NumTrianglesAllocated = 0;
  m_NumTriangles = 0;
}

bool ContourGeometry::doubleFloatArray(GLfloat *&array, unsigned int size) {
  // double a single float array, copy over the values, and delete the old array
  GLfloat *newArray = new GLfloat[size * 2];
  if (!newArray) {
    return false;
  }
  unsigned int c;
  for (c = 0; c < size; c++) {
    newArray[c] = array[c];
  }
  delete[] array;
  array = newArray;
  return true;
}

bool ContourGeometry::doubleIntArray(GLuint *&array, unsigned int size) {
  // double a single uint array, copy over the values, and delete the old array
  GLuint *newArray = new GLuint[size * 2];
  if (!newArray) {
    return false;
  }
  unsigned int c;
  for (c = 0; c < size; c++) {
    newArray[c] = array[c];
  }
  delete[] array;
  array = newArray;
  return true;
}

bool ContourGeometry::addTriangle(unsigned int v1, unsigned int v2,
                                  unsigned int v3) {
  // add a single triangle, doubling the array if necessary
  if (m_NumTriangles == m_NumTrianglesAllocated) {
    if (!doubleTriangleBuffers()) {
      return false;
    }
  }
  m_Triangles[m_NumTriangles * 3 + 0] = v1;
  m_Triangles[m_NumTriangles * 3 + 1] = v2;
  m_Triangles[m_NumTriangles * 3 + 2] = v3;
  m_NumTriangles++;
  return true;
}

void ContourGeometry::doInterpolation(GLfloat isovalue) {
  if (!m_InterpolationDone) {
    unsigned int c;
    if (m_UseColors) {
      for (c = 0; c < m_NumVertices; c++) {
        interpArray(m_Vertex1 + c * 4, m_Vertex2 + c * 4, m_Normal1 + c * 3,
                    m_Normal2 + c * 3, m_Vertex1[c * 4 + 3],
                    m_Vertex2[c * 4 + 3], isovalue);
      }
    } else {
      for (c = 0; c < m_NumVertices; c++) {
        interpArray(m_Vertex1 + c * 4, m_Vertex2 + c * 4, m_Normal1 + c * 3,
                    m_Normal2 + c * 3, m_Vertex1[c * 4 + 3],
                    m_Vertex2[c * 4 + 3], isovalue);
      }
    }
    m_InterpolationDone = true;
  }
}
