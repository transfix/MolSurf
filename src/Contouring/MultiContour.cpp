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
#include <Contouring/MultiContour.h>

using CCVContouring::MultiContour;

MultiContour::MultiContour() { setDefaults(); }

MultiContour::~MultiContour() {
  removeAll();
  delete m_Matrix;
  delete m_SaveMatrix;
}

static inline double maxOfThree(double n1, double n2, double n3) {
  double max = (n1 > n2 ? n1 : n2);
  return (max > n3 ? max : n3);
}

void MultiContour::setData(unsigned char *data, unsigned int width,
                           unsigned int height, unsigned int depth,
                           double aspectX, double aspectY, double aspectZ,
                           double subMinX, double subMinY, double subMinZ,
                           double subMaxX, double subMaxY, double subMaxZ,
                           double minX, double minY, double minZ, double maxX,
                           double maxY, double maxZ) {
  m_SingleExtractor.setData(data, width, height, depth, aspectX, aspectY,
                            aspectZ, subMinX, subMinY, subMinZ, subMaxX,
                            subMaxY, subMaxZ);
  m_ContourExtractor = &m_SingleExtractor;
  resetContours();
  double maxAspect = maxOfThree(aspectX, aspectY, aspectZ);
  aspectX = aspectX / maxAspect;
  aspectY = aspectY / maxAspect;
  aspectZ = aspectZ / maxAspect;
  m_Matrix->reset();
  // scale to 0 to 1
  m_Matrix->preMultiplication(CCVOpenGLMath::Matrix::scale(
      (float)(1.0 / (double)(width - 1)), (float)(1.0 / (double)(height - 1)),
      (float)(1.0 / (double)(depth - 1))));
  m_Matrix->preMultiplication(CCVOpenGLMath::Matrix::scale(
      (float)((maxX - minX)), (float)((maxY - minY)), (float)((maxZ - minZ))));
  m_Matrix->preMultiplication(CCVOpenGLMath::Matrix::translation(
      (float)(minX), (float)(minY), (float)(minZ)));
  // set up the save matrix
  m_SaveMatrix->reset();
  // scale to 0 to 1
  m_SaveMatrix->preMultiplication(CCVOpenGLMath::Matrix::scale(
      (float)(1.0 / (double)(width - 1)), (float)(1.0 / (double)(height - 1)),
      (float)(1.0 / (double)(depth - 1))));
  // align
  m_SaveMatrix->preMultiplication(CCVOpenGLMath::Matrix::translation(
      (float)(-subMinX), (float)(-subMinY), (float)(-subMinZ)));
  // scale to min to max
  m_SaveMatrix->preMultiplication(CCVOpenGLMath::Matrix::scale(
      (float)((maxX - minX) / (subMaxX - subMinX)),
      (float)((maxY - minY) / (subMaxY - subMinY)),
      (float)((maxZ - minZ) / (subMaxZ - subMinZ))));
  // translate
  m_SaveMatrix->preMultiplication(CCVOpenGLMath::Matrix::translation(
      (float)(minX), (float)(minY), (float)(minZ)));
  m_DataLoaded = true;
}

void MultiContour::setData(unsigned char *data, unsigned char *red,
                           unsigned char *green, unsigned char *blue,
                           unsigned int width, unsigned int height,
                           unsigned int depth, double aspectX, double aspectY,
                           double aspectZ, double subMinX, double subMinY,
                           double subMinZ, double subMaxX, double subMaxY,
                           double subMaxZ, double minX, double minY,
                           double minZ, double maxX, double maxY, double maxZ) {
  m_RGBAExtractor.setData(data, red, green, blue, width, height, depth, aspectX,
                          aspectY, aspectZ, subMinX, subMinY, subMinZ, subMaxX,
                          subMaxY, subMaxZ);
  m_ContourExtractor = &m_RGBAExtractor;
  resetContours();
  double maxAspect = maxOfThree(aspectX, aspectY, aspectZ);
  aspectX = aspectX / maxAspect;
  aspectY = aspectY / maxAspect;
  aspectZ = aspectZ / maxAspect;
  m_Matrix->reset();
  m_Matrix->preMultiplication(CCVOpenGLMath::Matrix::scale(
      (float)(1.0 / (double)(width - 1)), (float)(1.0 / (double)(height - 1)),
      (float)(1.0 / (double)(depth - 1))));
  m_Matrix->preMultiplication(CCVOpenGLMath::Matrix::scale(
      (float)((maxX - minX)), (float)((maxY - minY)), (float)((maxZ - minZ))));
  m_Matrix->preMultiplication(CCVOpenGLMath::Matrix::translation(
      (float)(minX), (float)(minY), (float)(minZ)));
  // set up the save matrix
  m_SaveMatrix->reset();
  // scale to 0 to 1
  m_SaveMatrix->preMultiplication(CCVOpenGLMath::Matrix::scale(
      (float)(1.0 / (double)(width - 1)), (float)(1.0 / (double)(height - 1)),
      (float)(1.0 / (double)(depth - 1))));
  // align
  m_SaveMatrix->preMultiplication(CCVOpenGLMath::Matrix::translation(
      (float)(-subMinX), (float)(-subMinY), (float)(-subMinZ)));
  // scale to min to max
  m_SaveMatrix->preMultiplication(CCVOpenGLMath::Matrix::scale(
      (float)((maxX - minX) / (subMaxX - subMinX)),
      (float)((maxY - minY) / (subMaxY - subMinY)),
      (float)((maxZ - minZ) / (subMaxZ - subMinZ))));
  // translate
  m_SaveMatrix->preMultiplication(CCVOpenGLMath::Matrix::translation(
      (float)(minX), (float)(minY), (float)(minZ)));
  m_DataLoaded = true;
}

void MultiContour::addContour(int ID, float isovalue, float R, float G,
                              float B) {
  MultiContourNode *newNode =
      new MultiContourNode(ID, isovalue, R, G, B, m_WireframeRender, m_Head);
  m_Head = newNode;
}

void MultiContour::removeContour(int ID) {
  MultiContourNode *toDelete;
  MultiContourNode *node = m_Head;
  // empty list
  if (m_Head == 0) {
    return;
  }
  // head is a special case
  if (m_Head->m_Contour.getID() == ID) {
    toDelete = m_Head;
    m_Head = m_Head->m_Next;
    toDelete->m_Next = 0;
    delete toDelete;
    return;
  }
  // do remainder of list
  while (node->m_Next && node->m_Next->m_Contour.getID() != ID) {
    node = node->m_Next;
  }
  // remove the node
  if (node->m_Next) {
    toDelete = node->m_Next;
    node->m_Next = node->m_Next->m_Next;
    toDelete->m_Next = 0;
    delete toDelete;
    return;
  }
}

void MultiContour::removeAll() {
  delete m_Head;
  m_Head = 0;
}

void MultiContour::setIsovalue(int ID, float isovalue) {
  MultiContourNode *node = m_Head;
  while (node && node->m_Contour.getID() != ID) {
    node = node->m_Next;
  }
  if (node) {
    node->m_Contour.setIsovalue(isovalue);
  }
}

void MultiContour::setColor(int ID, float R, float G, float B) {
  MultiContourNode *node = m_Head;
  // if we're in RGBA mode, we don't want to overwrite the contour's color
  // clobberColor == false -> (m_ContourExtractor == &m_RGBAExtractor)
  bool clobberColor = (m_ContourExtractor == &m_SingleExtractor);
  while (node && node->m_Contour.getID() != ID) {
    node = node->m_Next;
  }
  if (node) {
    node->m_Contour.setSingleColor(R, G, B, clobberColor);
  }
}

void MultiContour::setWireframeMode(bool state) {
  m_WireframeRender = state;
  if (m_DataLoaded) {
    MultiContourNode *node = m_Head;
    while (node) {
      node->m_Contour.setWireframeMode(state);
      node = node->m_Next;
    }
  }
}

void MultiContour::renderContours() const {
  if (m_DataLoaded) {
    MultiContourNode *node = m_Head;
    GLint matrixMode;
    glGetIntegerv(GL_MATRIX_MODE, &matrixMode);
    glMatrixMode(GL_MODELVIEW);
    float modelview_matrix[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, modelview_matrix);
    glPushMatrix();
    /*		// get it to right pos by translating and scaling
                    // do it in right OpenGL order
                    float largestSide = getLargestSide();

                    // 3. Shift by min vals
                    glTranslatef( m_MinX, m_MinY, m_MinZ );
                    // 2. Scale by largest dimension
                    glScalef( largestSide, largestSide, largestSide );
                    // 1. Shift by 0.5 in max, something smaller elsewhere to
       make it go from 0 - <=1 glTranslatef( (m_MaxX-m_MinX)/(largestSide
       * 2.0), (m_MaxY-m_MinY)/(largestSide * 2.0), (m_MaxZ-m_MinZ)/(largestSide
       * 2.0) );
    */
    if (m_Matrix) {
      glMultMatrixf(m_Matrix->getMatrix());
    }
    while (node) {
      node->m_Contour.renderContour(*m_ContourExtractor);
      node = node->m_Next;
    }
    glPopMatrix();
    glMatrixMode(matrixMode);
  }
}

void MultiContour::forceExtraction() const {
  if (m_DataLoaded) {
    MultiContourNode *node = m_Head;
    while (node) {
      node->m_Contour.extract(*m_ContourExtractor);
      node = node->m_Next;
    }
  }
}

bool MultiContour::render(double a, double b, double c, double d) {
  double plane0[] = {a, b, c, d};
  glClipPlane(GL_CLIP_PLANE0, plane0);
  glEnable(GL_CLIP_PLANE0);
  renderContours();
  glDisable(GL_CLIP_PLANE0);
  return true;
}

Geometry *MultiContour::getGeometry() const {
  int nextVert = 0, nextTri = 0;
  Geometry *geometry = 0;
  if (m_DataLoaded && m_Head) {
    geometry = new Geometry;
    int numVerts = getNumVerts();
    int numTris = getNumTris();
    geometry->AllocateTris(numVerts, numTris);
    if (m_Head->m_Contour.useColors()) {
      geometry->AllocateTriVertColors();
    }
    MultiContourNode *node = m_Head;
    while (node) {
      node->m_Contour.addToGeometry(*m_ContourExtractor, geometry,
                                    *m_SaveMatrix, nextVert, nextTri);
      node = node->m_Next;
    }
    geometry->SetTriNormalsReady();
  }
  return geometry;
}

int MultiContour::getNumVerts() const {
  int numVerts = 0;
  if (m_DataLoaded) {
    MultiContourNode *node = m_Head;
    while (node) {
      numVerts += node->m_Contour.getNumVerts(*m_ContourExtractor);
      node = node->m_Next;
    }
  }
  return numVerts;
}

int MultiContour::getNumTris() const {
  int numVerts = 0;
  if (m_DataLoaded) {
    MultiContourNode *node = m_Head;
    while (node) {
      numVerts += node->m_Contour.getNumTris(*m_ContourExtractor);
      node = node->m_Next;
    }
  }
  return numVerts;
}

void MultiContour::setDefaults() {
  m_Head = 0;
  m_DataLoaded = false;
  m_WireframeRender = false;
  m_Matrix = new CCVOpenGLMath::Matrix;
  m_SaveMatrix = new CCVOpenGLMath::Matrix;
}

void MultiContour::resetContours() {
  MultiContourNode *node = m_Head;
  while (node) {
    node->m_Contour.resetContour();
    node = node->m_Next;
  }
}
