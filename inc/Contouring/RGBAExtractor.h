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
    AFX_RGBAEXTRACTOR_H__79C7CB25_6CC3_4FF6_870D_AA1018CEABF5__INCLUDED_)
#define AFX_RGBAEXTRACTOR_H__79C7CB25_6CC3_4FF6_870D_AA1018CEABF5__INCLUDED_

#include <Contouring/ContourExtractor.h>
#include <Contouring/ContourGeometry.h>
#include <Contouring/MarchingCubesBuffers.h>
#include <UsefulMath/Matrix.h>

namespace CCVContouring {
class RGBAExtractor : public ContourExtractor {
public:
  RGBAExtractor();
  virtual ~RGBAExtractor();
  virtual void extractContour(ContourGeometry *contourGeometry, float isovalue,
                              float R, float G, float B) const;
  void setData(unsigned char *data, unsigned char *red, unsigned char *green,
               unsigned char *blue, unsigned int width, unsigned int height,
               unsigned int depth, double aspectX, double aspectY,
               double aspectZ, double subMinX, double subMinY, double subMinZ,
               double subMaxX, double subMaxY, double subMaxZ);

protected:
  void setDefaults();
  void classifyVertices(unsigned int k, unsigned int *cacheMemory,
                        GLfloat isovalue) const;
  inline void getNormal(unsigned char *data, unsigned int i, unsigned int j,
                        unsigned int k, float &nx, float &ny, float &nz) const;
  inline void isCached(bool *cachedTable, unsigned int i, unsigned int j,
                       unsigned int k, unsigned int width,
                       unsigned int height) const;
  inline unsigned int determineCase(unsigned char *data,
                                    unsigned int *offsetTable,
                                    unsigned int index, GLfloat isovalue) const;
  inline void buildOffsetTable(unsigned int *table, unsigned int width,
                               unsigned int height) const;
  inline void buildWithinSliceOffsetTable(unsigned int *table,
                                          unsigned int width) const;
  inline void buildEdgeCacheOffsetTable(unsigned int *table,
                                        unsigned int width) const;
  inline void computeVertFromOffset(GLfloat &vx, GLfloat &vy, GLfloat &vz,
                                    unsigned int offset, unsigned int width,
                                    unsigned int height) const;
  void calcOffsets(unsigned int *offsets, unsigned int width,
                   unsigned int height);
  unsigned char *m_Data;
  unsigned char *m_Red;
  unsigned char *m_Green;
  unsigned char *m_Blue;
  // dimensions of the data
  unsigned int m_Width;
  unsigned int m_Height;
  unsigned int m_Depth;
  // buffers used to speed up marching cubes
  MarchingCubesBuffers *m_Buffers;
};
}; // namespace CCVContouring

inline void
CCVContouring::RGBAExtractor::getNormal(unsigned char *data, unsigned int i,
                                        unsigned int j, unsigned int k,
                                        float &nx, float &ny, float &nz) const {
  double dx, dy, dz, length;
  int negZOffset, negYOffset, negXOffset;
  int posZOffset, posYOffset, posXOffset;
  double scaleX = 1.0, scaleY = 1.0, scaleZ = 1.0;
  unsigned int inputindex;
  if (k == 0) // border offset
  {
    negZOffset = 0;
    posZOffset = m_Width * m_Height;
  } else if (k == m_Depth - 1) // border offset
  {
    negZOffset = -(int)(m_Width * m_Height);
    posZOffset = 0;
  } else // normal offset
  {
    negZOffset = -(int)(m_Width * m_Height);
    posZOffset = m_Width * m_Height;
    scaleZ = 0.5;
  }
  if (j == 0) // border offset
  {
    negYOffset = 0;
    posYOffset = m_Width;
  } else if (j == m_Height - 1) // border offset
  {
    negYOffset = -(int)m_Width;
    posYOffset = 0;
  } else // normal offset
  {
    negYOffset = -(int)m_Width;
    posYOffset = m_Width;
    scaleY = 0.5;
  }
  if (i == 0) // border offset
  {
    negXOffset = 0;
    posXOffset = 1;
  } else if (i == m_Width - 1) // border offset
  {
    negXOffset = -1;
    posXOffset = 0;
  } else // normal offset
  {
    negXOffset = -1;
    posXOffset = 1;
    scaleX = 0.5;
  }
  inputindex = k * m_Width * m_Height + j * m_Width + i;
  dx = (double)(data[inputindex + negXOffset] - data[inputindex + posXOffset]) *
       scaleX;
  dy = (double)(data[inputindex + negYOffset] - data[inputindex + posYOffset]) *
       scaleY;
  dz = (double)(data[inputindex + negZOffset] - data[inputindex + posZOffset]) *
       scaleZ;
  length = sqrt(dx * dx + dy * dy + dz * dz);
  if (length != 0) {
    nx = (float)(dx / length);
    ny = (float)(dy / length);
    nz = (float)(dz / length);
  }
}

inline void CCVContouring::RGBAExtractor::buildEdgeCacheOffsetTable(
    unsigned int *table, unsigned int width) const {
  // calculates where edges are cached in their respective table
  // is an offset from k*width*height + j*width + i
  // edge 0
  table[0] = 0;
  // edge 1
  table[1] = 1;
  // edge 2
  table[2] = 0;
  // edge 3
  table[3] = 0;
  // edge 4
  table[4] = width;
  // edge 5
  table[5] = width + 1;
  // edge 6
  table[6] = width;
  // edge 7
  table[7] = width;
  // edge 8
  table[8] = 0;
  // edge 9
  table[9] = 1;
  // edge 10
  table[10] = 0;
  // edge 11
  table[11] = 1;
}

inline void CCVContouring::RGBAExtractor::isCached(
    bool *cachedTable, unsigned int i, unsigned int j, unsigned int k,
    unsigned int width, unsigned int height) const {
  // bool ret = false;
  // ret = j!=0 || k!=0;
  // edge 0
  cachedTable[0] = j != 0 || k != 0;
  // edge 1
  cachedTable[1] = j != 0;
  // edge 2
  cachedTable[2] = j != 0;
  // edge 3
  cachedTable[3] = i != 0 || j != 0;
  // edge 4
  cachedTable[4] = k != 0;
  // edge 5
  cachedTable[5] = false;
  // edge 6
  cachedTable[6] = false;
  // edge 7
  cachedTable[7] = i != 0;
  // edge 8
  cachedTable[8] = k != 0 || i != 0;
  // edge 9
  cachedTable[9] = k != 0;
  // edge 10
  cachedTable[10] = i != 0;
  // edge 11
  cachedTable[11] = false;
}

inline void CCVContouring::RGBAExtractor::buildOffsetTable(
    unsigned int *table, unsigned int width, unsigned int height) const {
  // determine the offset of each vertex of the cube give the width and height
  // of the data
  unsigned int c;
  for (c = 0; c < 8; c++) {
    table[c] = verts[c][2] + width * verts[c][1] + width * height * verts[c][0];
  }
}

inline void CCVContouring::RGBAExtractor::buildWithinSliceOffsetTable(
    unsigned int *table, unsigned int width) const {
  // determine the offset of each vertex of the cube give the width and height
  // of the data
  unsigned int c;
  for (c = 0; c < 8; c++) {
    table[c] = verts[c][2] + width * verts[c][1];
  }
}

inline void CCVContouring::RGBAExtractor::computeVertFromOffset(
    GLfloat &vx, GLfloat &vy, GLfloat &vz, unsigned int offset,
    unsigned int width, unsigned int height) const {
  vx = (GLfloat)(offset % (width));
  vy = (GLfloat)((offset / width) % height);
  vz = (GLfloat)(offset / (width * height));
}

#endif
