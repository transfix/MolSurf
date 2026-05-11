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
#include <Blurmaps/SkinRegion2.h>
#include <PDBParser/FlattenGOA.h>

using namespace CVCUtility;

SkinRegion2::SkinRegion2() {
  m_SasInteriorVertices = 0;
  m_SasBoundaryVertices = 0;
  m_VDWInteriorVertices = 0;
  m_VDWBoundaryVertices = 0;
  m_SasBoundaryCells = 0;
  m_IndexOfBoundaryCell = 0;
  m_SphereIndices = 0;
  m_SubdividedNodeSasInteriorVertex = 0;
  m_SubdividedNodeSasBoundaryCell = 0;
  m_Regions = 0;
  m_SES = 0;
  m_Red = 0;
  m_Green = 0;
  m_Blue = 0;
  enlargedRadius = 0.0;
}

SkinRegion2::~SkinRegion2() {
  if (m_SasInteriorVertices) {
    delete[] m_SasInteriorVertices;
    m_SasInteriorVertices = 0;
  }
  if (m_SasBoundaryVertices) {
    delete[] m_SasBoundaryVertices;
    m_SasBoundaryVertices = 0;
  }
  if (m_VDWInteriorVertices) {
    delete[] m_VDWInteriorVertices;
    m_VDWInteriorVertices = 0;
  }
  if (m_VDWBoundaryVertices) {
    delete[] m_VDWBoundaryVertices;
    m_VDWBoundaryVertices = 0;
  }
  if (m_SasBoundaryCells) {
    delete[] m_SasBoundaryCells;
    m_SasBoundaryCells = 0;
  }
  if (m_IndexOfBoundaryCell) {
    delete[] m_IndexOfBoundaryCell;
    m_IndexOfBoundaryCell = 0;
  }
  if (m_SphereIndices) {
    delete[] m_SphereIndices;
    m_SphereIndices = 0;
  }
  if (m_SubdividedNodeSasInteriorVertex) {
    delete[] m_SubdividedNodeSasInteriorVertex;
    m_SubdividedNodeSasInteriorVertex = 0;
  }
  if (m_SubdividedNodeSasBoundaryCell) {
    delete[] m_SubdividedNodeSasBoundaryCell;
    m_SubdividedNodeSasBoundaryCell = 0;
  }
  if (m_Regions) {
    delete[] m_Regions;
    m_Regions = 0;
  }
  {
    int n = m_AtomList.size();
    int i;
    for (i = 0; i < n; i++) {
      delete m_AtomList[i];
    }
    m_AtomList.clear();
  }
}

/*bool SkinRegion2::intersects( double* center, double radius, int iIndex, int
jIndex, int kIndex )
{
        double Bmin[3], Bmax[3];
        Bmin[0] = m_Orig[0] + m_Span[0] * iIndex;
        Bmin[1] = m_Orig[1] + m_Span[1] * jIndex;
        Bmin[2] = m_Orig[2] + m_Span[2] * kIndex;

        Bmax[0] = Bmin[0] + m_Span[0];
        Bmax[1] = Bmin[1] + m_Span[1];
        Bmax[2] = Bmin[2] + m_Span[2];

        double dmin = 0;
        double dmax = 0;
        bool face = false;
        int i;
        for( i = 0; i < 3; i++ )
        {
                double a = sqrt( center[i] - Bmin[i] );
                double b = sqrt( center[i] - Bmax[i] );
                if( a>b ) dmax += a;
                else dmax += b;

                if( center[i] < Bmin[i] ) {
                        face = true;
                        dmin += a;
                }
                else if( center[i] > Bmax[i] ) {
                        face = true;
                        dmin += b;
                }
                else
                {
                        double minval = a;
                        if( b<a ) minval = b;
                        if( minval <= radius*radius ) face = true;
                }
        }
        if(face && ( dmin <= radius*radius ) && ( radius*radius <= dmax)) return
true; return false;
}*/

void SkinRegion2::setEnlargedRadius(double r) { enlargedRadius = r; }

bool SkinRegion2::intersects(double *center, double radius, int iIndex,
                             int jIndex, int kIndex) {
  double Bmin[3], Bmax[3];
  Bmin[0] = m_Orig[0] + m_Span[0] * iIndex;
  Bmin[1] = m_Orig[1] + m_Span[1] * jIndex;
  Bmin[2] = m_Orig[2] + m_Span[2] * kIndex;
  Bmax[0] = Bmin[0] + m_Span[0];
  Bmax[1] = Bmin[1] + m_Span[1];
  Bmax[2] = Bmin[2] + m_Span[2];
  double dmin = 0;
  int i;
  for (i = 0; i < 3; i++) {
    if (center[i] < Bmin[i]) {
      dmin += (center[i] - Bmin[i]) * (center[i] - Bmin[i]);
    } else if (center[i] > Bmax[i]) {
      dmin += (center[i] - Bmax[i]) * (center[i] - Bmax[i]);
    }
  }
  if (dmin <= radius * radius) {
    return true;
  }
  return false;
}

/*!
Cell is defined by the corner vertex.
Algorithm:
1. All vertices inside, return false
2. One vertex intersects, other doesnt, return true
3. Center of circle within edge, dist from center to edge < radius, return true
4. Center of circle within face, dist from center to face < radius, return true
5. return false
*/
bool SkinRegion2::intersects(double *center, double radius, int cell) {
  {
    int cellIndices[8];
    TrilinearGrid::getCellVertices(cell, cellIndices, m_Dim);
    int v;
    bool oneIn = false;
    bool oneOut = false;
    for (v = 0; v < 8; v++) {
      int i, j, k;
      TrilinearGrid::vtx2xyz(cellIndices[v], i, j, k, m_Dim);
      double xpos = m_Orig[0] + m_Span[0] * i;
      double ypos = m_Orig[1] + m_Span[1] * j;
      double zpos = m_Orig[2] + m_Span[2] * k;
      double dist = (xpos - center[0]) * (xpos - center[0]) +
                    (ypos - center[1]) * (ypos - center[1]) +
                    (zpos - center[2]) * (zpos - center[2]);
      if (dist < radius * radius) {
        oneIn = true;
      } else {
        oneOut = true;
      }
    }
    // case 1
    if (!oneOut) {
      return false;
    }
    // case 2
    if (oneOut && oneIn) {
      return true;
    }
  }
  // case 3
  return false;
}

bool SkinRegion2::isBoundaryCell(bool *sasCellInterior) {
  bool interiorFlag, exteriorFlag;
  interiorFlag = exteriorFlag = false;
  int i;
  for (i = 0; i < 8; i++) {
    if (sasCellInterior[i] == true) {
      interiorFlag = true;
    } else {
      exteriorFlag = true;
    }
  }
  return (interiorFlag & exteriorFlag);
}

/*!
Once vertices are classified as interior or not, we find and mark the boundary
cells and mark all vertices belonging to a boundary cell as a boundary vertex

  Also first classify m_Regions as VDW boundary vertices and then as SAS
boundary vertices. This priority was chosen arbitrarily
*/
void SkinRegion2::getBoundaryVerticesCells() {
  {
    int i;
    for (i = 0; i < m_NumVerts; i++) {
      m_SasBoundaryVertices[i] = false;
    }
    for (i = 0; i < m_NumCells; i++) {
      m_SasBoundaryCells[i] = false;
    }
  }
  int cellVertexArray[8];
  {
    bool vdwCellInterior[8];
    int i;
    for (i = 0; i < m_NumCells; i++) {
      TrilinearGrid::getCellVertices(i, cellVertexArray, m_Dim);
      int j;
      for (j = 0; j < 8; j++) {
        vdwCellInterior[j] = m_VDWInteriorVertices[cellVertexArray[j]];
      }
      if (isBoundaryCell(vdwCellInterior)) {
        for (j = 0; j < 8; j++) {
          m_Regions[cellVertexArray[j]] = VDW_SURFACE;
          m_VDWBoundaryVertices[cellVertexArray[j]] = true;
        }
      }
    }
  }
  {
    int i;
    bool sasCellInterior[8];
    for (i = 0; i < m_NumCells; i++) {
      TrilinearGrid::getCellVertices(i, cellVertexArray, m_Dim);
      int j;
      for (j = 0; j < 8; j++) {
        sasCellInterior[j] = m_SasInteriorVertices[cellVertexArray[j]];
      }
      if (isBoundaryCell(sasCellInterior)) {
        for (j = 0; j < 8; j++) {
          if (m_Regions[cellVertexArray[j]] == SAS_VOLUME) {
            m_SasBoundaryVertices[cellVertexArray[j]] = true;
            m_Regions[cellVertexArray[j]] = SAS_SURFACE;
          } else {
            // eiher in vdw, or outside! if outside, mark it as outside
            // if( m_Regions[cellVertexArray[j]] != VDW_VOLUME &&
            // m_Regions[cellVertexArray[j]] != VDW_SURFACE )
            //        m_Regions[cellVertexArray[j]] = SAS_OUTER_SURFACE;
          }
        }
        m_SasBoundaryCells[i] = true;
      }
    }
  }
}

// Given two grid points, we return the distance between the two of them in real
// space Craig: this is called intensively so I'm inlining it
double inline SkinRegion2::getDist(double x1, double y1, double z1, double x2,
                                   double y2, double z2) {
  double x = fabs(x1 - x2) * m_Span[0];
  double y = fabs(y1 - y2) * m_Span[1];
  double z = fabs(z1 - z2) * m_Span[2];
  return sqrt(x * x + y * y + z * z);
}

/*!
Given a set of cells, we give an index to each that is a boundary cell.
We also make the indices of other cells as -1.
We return the number of boundary cells.
*/
int SkinRegion2::fillBoundaryIndices() {
  int count = 0;
  int i;
  for (i = 0; i < m_NumCells; i++) {
    if (m_SasBoundaryCells[i]) {
      m_IndexOfBoundaryCell[i] = count;
      count++;
    } else {
      m_IndexOfBoundaryCell[i] = -1;
    }
  }
  return count;
}

/*!
Given a grid point, we return a bounding box to search for the SES from it. We
only need to look probe radius away
*/
void SkinRegion2::getBoundingBox(int x, int y, int z, int *imin, int *jmin,
                                 int *kmin, int *imax, int *jmax, int *kmax,
                                 int iProbeRadius, int jProbeRadius,
                                 int kProbeRadius) {
  *imin = x - iProbeRadius - 1;
  if (*imin < 0) {
    *imin = 0;
  }
  *jmin = y - jProbeRadius - 1;
  if (*jmin < 0) {
    *jmin = 0;
  }
  *kmin = z - kProbeRadius - 1;
  if (*kmin < 0) {
    *kmin = 0;
  }
  *imax = x + iProbeRadius + 1;
  if (*imax > m_Dim[0] - 1) {
    *imax = m_Dim[0] - 1;
  }
  *jmax = y + jProbeRadius + 1;
  if (*jmax > m_Dim[1] - 1) {
    *jmax = m_Dim[1] - 1;
  }
  *kmax = z + kProbeRadius + 1;
  if (*kmax > m_Dim[2] - 1) {
    *kmax = m_Dim[2] - 1;
  }
}

/*!
For each boundary cell, we have the sphere that pass through it.
For each sphere in a cell, we find the subvertices of the cell that are interior
to it.
*/
void SkinRegion2::getSubdividedNodesSasInteriorVertices() {
  int n = pow(2.0f, (m_Depth - 1)) + 1; // number of vertices in each side
  {
    int c = m_NumBoundaryCells * m_NumSmallVerts;
    for (int i = 0; i < c; i++) {
      m_SubdividedNodeSasInteriorVertex[i] = false;
    }
  }
  unsigned int small_dim[3];
  small_dim[0] = small_dim[1] = small_dim[2] = n;
  for (int c = 0; c < m_NumCells; c++) {
    if (!m_SasBoundaryCells[c]) {
      continue;
    }
    int cellIndex = m_IndexOfBoundaryCell[c];
    int celli, cellj, cellk;
    TrilinearGrid::cell2xyz(c, celli, cellj, cellk, m_Dim);
    for (int s = 0; s < m_SphereIndices[cellIndex].size(); s++) {
      PDBParser::Atom *at = m_AtomList[m_SphereIndices[cellIndex].at(s)];
      float radius = at->getRadius(m_RadiusType);
      radius += enlargedRadius;
      radius += m_ProbeRadius;
      double sx = at->m_Position[0];
      double sy = at->m_Position[1];
      double sz = at->m_Position[2];
      for (int sv = 0; sv < m_NumSmallVerts; sv++) {
        int i, j, k;
        TrilinearGrid::vtx2xyz(sv, i, j, k, small_dim);
        double vx = m_Orig[0] + (celli + i / ((double)n)) * m_Span[0];
        double vy = m_Orig[1] + (cellj + j / ((double)n)) * m_Span[1];
        double vz = m_Orig[2] + (cellk + k / ((double)n)) * m_Span[2];
        double dist = (sx - vx) * (sx - vx) + (sy - vy) * (sy - vy) +
                      (sz - vz) * (sz - vz);
        if (dist < radius * radius) {
          m_SubdividedNodeSasInteriorVertex[cellIndex * m_NumSmallVerts + sv] =
              true;
        }
      }
    }
  }
}

/*!
We already know which small vertices are in and out. We now classify the cells
as boundary or not.
*/
void SkinRegion2::getSubdividedNodesSasBoundaryCells() {
  int n = pow(2.0f, (m_Depth - 1)) + 1; // number of small verts in each side
  unsigned int small_dim[3];
  small_dim[0] = small_dim[1] = small_dim[2] = n;
  {
    int c = m_NumBoundaryCells * m_NumSmallCells;
    int i;
    for (i = 0; i < c; i++) {
      m_SubdividedNodeSasBoundaryCell[i] = false;
    }
  }
  int c;
  for (c = 0; c < m_NumCells; c++) {
    if (!m_SasBoundaryCells[c]) {
      continue;
    }
    int cellIndex = m_IndexOfBoundaryCell[c];
    int celli, cellj, cellk;
    TrilinearGrid::cell2xyz(c, celli, cellj, cellk, m_Dim);
    int sc;
    for (sc = 0; sc < m_NumSmallCells; sc++) {
      int cellVertexArray[8];
      TrilinearGrid::getCellVertices(sc, cellVertexArray, small_dim);
      bool sasCellInterior[8];
      int j;
      for (j = 0; j < 8; j++) {
        sasCellInterior[j] =
            m_SubdividedNodeSasInteriorVertex[cellIndex * m_NumSmallVerts +
                                              cellVertexArray[j]];
      }
      if (isBoundaryCell(sasCellInterior)) {
        m_SubdividedNodeSasBoundaryCell[cellIndex * m_NumSmallCells + sc] =
            true;
      }
    }
  }
}

void SkinRegion2::getSphereRadiusCenter(PDBParser::Atom *at, double *radius,
                                        double *center) {
  if (!at || !radius || !center) {
    return;
  }
  *radius = at->getRadius(m_RadiusType) + enlargedRadius;
  center[0] = at->m_Position[0];
  center[1] = at->m_Position[1];
  center[2] = at->m_Position[2];
}

/*!
For each sphere, find all the cells it intersects.
If it intersects a boundary cell, add its index to the m_SphereIndices array at
the right location
*/
void SkinRegion2::getBoundarySphereIndices() {
  int numAdded = 0;
  int numBCells = 0;
  int s;
  for (s = 0; s < m_AtomList.size(); s++) {
    PDBParser::Atom *at = m_AtomList[s];
    double radius, center[3];
    getSphereRadiusCenter(at, &radius, center);
    int bbmin[3], bbmax[3];
    {
      int i;
      for (i = 0; i < 3; i++) {
        bbmin[i] =
            floor((at->m_Position[i] - (radius + m_ProbeRadius) - m_Orig[i]) /
                  m_Span[i]);
        bbmin[i]--;
        if (bbmin[i] < 0) {
          bbmin[i] = 0;
        }
        bbmax[i] =
            floor((at->m_Position[i] + (radius + m_ProbeRadius) - m_Orig[i]) /
                  m_Span[i]); // we take floor because we only need the lower
                              // vertex for the cell
        bbmax[i]++;
        if (bbmax[i] > m_Dim[i] - 2) {
          bbmax[i] = m_Dim[i] - 2;
        }
      }
    }
    int i, j, k;
    for (i = bbmin[0]; i <= bbmax[0]; i++) {
      for (j = bbmin[1]; j <= bbmax[1]; j++) {
        for (k = bbmin[2]; k <= bbmax[2]; k++) {
          int c = TrilinearGrid::xyz2cell(i, j, k, m_Dim);
          if (!m_SasBoundaryCells[c]) {
            continue;
          }
          numBCells++;
          // if( intersects( center, radius+m_ProbeRadius, c ) )
          if (intersects(center, radius + m_ProbeRadius, i, j, k)) {
            m_SphereIndices[m_IndexOfBoundaryCell[c]].push_back(s);
            numAdded++;
          }
        }
      }
    }
  }
  // printf("Number added = %d\n", numAdded );
}

void SkinRegion2::getSasInterior() {
  {
    int c;
    for (c = 0; c < m_NumVerts; c++) {
      m_SasInteriorVertices[c] = false;
      m_VDWInteriorVertices[c] = false;
      m_VDWBoundaryVertices[c] = false;
      m_Regions[c] = EXTERIOR_VOLUME;
      // m_SES[c] = m_ProbeRadius*4.0;
      m_SES[c] = 0;
    }
  }
  unsigned int zoff = 0;
  int a;
  for (a = 0; a < m_AtomList.size(); a++) {
    PDBParser::Atom *at = m_AtomList[a];
    double radius, c[3];
    getSphereRadiusCenter(at, &radius, c);
    // calculate the maximum radius for this atom
    double maxRad = radius + m_ProbeRadius + 1;
    // compute the dataset coordinates of the atom's center
    unsigned int amax[3], amin[3];
    {
      int j;
      double apos[3];
      for (j = 0; j < 3; j++) {
        apos[j] = (at->m_Position[j] - m_Orig[j]) / m_Span[j];
        apos[j] = ((apos[j] - floor(apos[j])) >= 0.5) ? ceil(apos[j])
                                                      : floor(apos[j]);
      }
      // then compute the bounding box of the atom (maxRad^3)
      for (j = 0; j < 3; j++) {
        int tempAMin = (int)(apos[j] - (maxRad / m_Span[j]) - 1);
        tempAMin = (tempAMin < 0) ? 0 : tempAMin;
        amin[j] = tempAMin;
        int tempAMax = (int)(apos[j] + (maxRad / m_Span[j]) + 1);
        tempAMax = (tempAMax > m_Dim[j]) ? m_Dim[j] : tempAMax;
        amax[j] = tempAMax;
      }
    }
    // begin blurring kernel
    int i, j, k;
    for (i = amin[0]; i < amax[0]; i++) {
      for (j = amin[1]; j < amax[1]; j++) {
        for (k = amin[2]; k < amax[2]; k++) {
          int v = TrilinearGrid::xyz2vtx(i, j, k, m_Dim);
          float pnt[3];
          pnt[0] = m_Orig[0] + i * m_Span[0];
          pnt[1] = m_Orig[1] + j * m_Span[1];
          pnt[2] = m_Orig[2] + k * m_Span[2];
          double r = (c[0] - pnt[0]) * (c[0] - pnt[0]) +
                     (c[1] - pnt[1]) * (c[1] - pnt[1]) +
                     (c[2] - pnt[2]) * (c[2] - pnt[2]);
          if (r < (radius + m_ProbeRadius) * (radius + m_ProbeRadius)) {
            m_SasInteriorVertices[v] = true;
            if (m_Regions[v] != VDW_VOLUME) {
              m_Regions[v] = SAS_VOLUME;
            }
          }
          if (r < (radius) * (radius)) {
            m_VDWInteriorVertices[v] = true;
            m_Regions[v] = VDW_VOLUME;
          }
        }
      }
    }
    // end blurring kernel
  }
}

double SkinRegion2::getDistToSas2(int vertx, int verty, int vertz,
                                  int neighborx, int neighbory, int neighborz) {
  return getDist(vertx, verty, vertz, neighborx + 0.5, neighbory + 0.5,
                 neighborz + 0.5);
}

void SkinRegion2::getColors(PDBParser::Atom *atom, unsigned char *red,
                            unsigned char *green, unsigned char *blue) {
  if (!atom || !red || !green || !blue) {
    return;
  }
  // aAlpha = ((float)(atom->m_UserDefinedColor & 0xff)) / 255.0;
  *red = (unsigned char)((atom->m_UserDefinedColor & 0xff000000) >> 24);
  *green = (unsigned char)((atom->m_UserDefinedColor & 0xff0000) >> 16);
  *blue = (unsigned char)((atom->m_UserDefinedColor & 0xff00) >> 8);
}

/*!
Find the distance of the grid point vertex to the SAS. Search in the small cells
of the given cell Take the distance to be the distance to the center of the
small cell.
*/
double SkinRegion2::getDistToSas(int vertx, int verty, int vertz, int neighborx,
                                 int neighbory, int neighborz,
                                 unsigned char *red, unsigned char *green,
                                 unsigned char *blue, bool color) {
  double curDist = 10000000;
  int cell = TrilinearGrid::xyz2cell(neighborx, neighbory, neighborz, m_Dim);
  int cellIndex = m_IndexOfBoundaryCell[cell];
  {
    if (m_SphereIndices[cellIndex].size() == 0) {
      return curDist;
    }
    if (m_SphereIndices[cellIndex].size() == 1) {
      PDBParser::Atom *at = m_AtomList[m_SphereIndices[cellIndex].at(0)];
      float radius = at->getRadius() + enlargedRadius + m_ProbeRadius;
      double center[3];
      center[0] = at->m_Position[0];
      center[1] = at->m_Position[1];
      center[2] = at->m_Position[2];
      double point[3];
      point[0] = m_Orig[0] + vertx * m_Span[0];
      point[1] = m_Orig[1] + verty * m_Span[1];
      point[2] = m_Orig[2] + vertz * m_Span[2];
      double bmin[3];
      bmin[0] = m_Orig[0] + neighborx * m_Span[0];
      bmin[1] = m_Orig[1] + neighbory * m_Span[1];
      bmin[2] = m_Orig[2] + neighborz * m_Span[2];
      double bmax[3];
      bmax[0] = bmin[0] + m_Span[0];
      bmax[1] = bmin[1] + m_Span[1];
      bmax[2] = bmin[2] + m_Span[2];
      SphericalPatch *sphericalPatch = new SphericalPatch();
      double d =
          sphericalPatch->getClosestDist(radius, center, point, bmin, bmax);
      delete sphericalPatch;
      if (color) {
        getColors(at, red, green, blue);
      }
      if (d >= 0) {
        return d;
      }
      return curDist;
    }
  }
  if (color) {
    PDBParser::Atom *at = m_AtomList[m_SphereIndices[cellIndex].at(0)];
    getColors(at, red, green, blue);
  }
  unsigned int dim[3];
  dim[0] = dim[1] = dim[2] = pow(2.0f, m_Depth - 1) + 1;
  for (int i = 0; i < m_NumSmallCells; i++) {
    int smallCellIndex = cellIndex * m_NumSmallCells + i;
    if (!m_SubdividedNodeSasBoundaryCell[smallCellIndex]) {
      continue;
    }
    int cellx, celly, cellz;
    TrilinearGrid::cell2xyz(i, cellx, celly, cellz, dim);
    if ((cellx + 0.5) / ((double)dim[0]) > 1.01) {
      error("Error getting distance to SAS!");
    }
    if ((celly + 0.5) / ((double)dim[1]) > 1.01) {
      error("Error getting distance to SAS!");
    }
    if ((cellz + 0.5) / ((double)dim[2]) > 1.01) {
      error("Error getting distance to SAS!");
    }
    double d = getDist(vertx, verty, vertz,
                       neighborx + (cellx + 0.5) / ((double)dim[0]),
                       neighbory + (celly + 0.5) / ((double)dim[1]),
                       neighborz + (cellz + 0.5) / ((double)dim[2]));
    if (curDist > d) {
      curDist = d;
    }
  }
  return curDist;
}

/*!
For each vertex belonging to the S_SAS, V_SAS, S_VDW, S_SES, S_VDW:
Search in a region around it.
If you found a boundary cell
Search within the boundary cell
Store the shortest distance from the vertex to the center of the small cell
containing a S_SAS
*/
void SkinRegion2::getSesSurfaceAndVolume() {
  int imin = 0, jmin = 0, kmin = 0, imax = 0, jmax = 0, kmax = 0;
  int iProbeRadius = ceil(m_ProbeRadius / m_Span[0]);
  int jProbeRadius = ceil(m_ProbeRadius / m_Span[1]);
  int kProbeRadius = ceil(m_ProbeRadius / m_Span[2]);
  float interiorDistVal = 10 * m_ProbeRadius;
  float maxExtDistVal = -m_ProbeRadius;
  for (int v = 0; v < m_NumVerts; v++) {
    if (m_Regions[v] == EXTERIOR_VOLUME) {
      m_SES[v] = maxExtDistVal;
      continue;
    }
    if (m_Regions[v] == VDW_VOLUME) {
      m_SES[v] = interiorDistVal;
      continue;
    }
    double minDistToSas = 100000000;
    unsigned char red = 0, green = 0, blue = 0;
    int x, y, z;
    TrilinearGrid::vtx2xyz(v, x, y, z, m_Dim);
    getBoundingBox(x, y, z, &imin, &jmin, &kmin, &imax, &jmax, &kmax,
                   iProbeRadius, jProbeRadius, kProbeRadius);
    for (int i = imin; i < imax; i++) {
      for (int j = jmin; j < jmax; j++) {
        for (int k = kmin; k < kmax; k++) {
          int neighborVert = TrilinearGrid::xyz2cell(i, j, k, m_Dim);
          if (!m_SasBoundaryCells[neighborVert]) {
            continue;
          }
          unsigned char curRed, curGreen, curBlue;
          double curMinDistToSas =
              getDistToSas(x, y, z, i, j, k, &curRed, &curGreen, &curBlue,
                           m_Red && m_Green && m_Blue);
          if (minDistToSas > curMinDistToSas) {
            minDistToSas = curMinDistToSas;
            red = curRed;
            green = curGreen, blue = curBlue;
          }
        }
      }
    }
    if (m_Red && m_Green && m_Blue) {
      m_Red[v] = red;
      m_Green[v] = green;
      m_Blue[v] = blue;
    }
    {
      if (minDistToSas > m_ProbeRadius) {
        m_Regions[v] = SES_VOLUME;
        m_SES[v] = minDistToSas;
        if (minDistToSas > interiorDistVal) {
          m_SES[v] = interiorDistVal;
        }
      } else {
        m_SES[v] = minDistToSas;
        double scaleFactor = (VDW_SURFACE - SAS_SURFACE) / (interiorDistVal);
        double dist = minDistToSas * scaleFactor;
        if (dist > VDW_SURFACE) {
          dist = VDW_SURFACE;
        }
        m_Regions[v] = SAS_SURFACE + dist;
      }
    }
    // printf("Finished vertex %d out of %d\n", v, m_NumVerts );
    /* // Craig: this is prohibitively expensive
    if(((v+1) % 500) == 0 || (v+1) == m_NumVerts)
    {
            printf("%2.2f%% done \r", 100.0*(v+1)/(float)m_NumVerts);
            fflush(stdout);
            }*/
  }
}

bool SkinRegion2::allocateSpace(bool withColor, SimpleVolumeData *vol) {
  vol->setDimensions(m_Dim);
  int numVerts = m_Dim[0] * m_Dim[1] * m_Dim[2];
  if (withColor) {
    m_Red = new unsigned char[numVerts], m_Green = new unsigned char[numVerts],
    m_Blue = new unsigned char[numVerts];
    m_SES = new float[numVerts];
    vol->setNumberOfVariables(4);
    vol->setData(0, m_Red);
    vol->setType(0, SimpleVolumeData::UCHAR);
    vol->setName(0, "red");
    vol->setData(1, m_Green);
    vol->setType(1, SimpleVolumeData::UCHAR);
    vol->setName(1, "green");
    vol->setData(2, m_Blue);
    vol->setType(2, SimpleVolumeData::UCHAR);
    vol->setName(2, "blue");
    vol->setData(3, m_SES);
    vol->setType(3, SimpleVolumeData::FLOAT);
    vol->setName(3, "density");
  } else {
    m_SES = new float[numVerts];
    vol->setNumberOfVariables(1);
    vol->setData(0, m_SES);
    vol->setType(0, SimpleVolumeData::FLOAT);
    vol->setName(0, "density");
  }
  vol->setMinExtent(m_Min);
  vol->setMaxExtent(m_Max);
  vol->setVariablesToZero();
  return true;
}

bool SkinRegion2::getColoredAtomList(
    PDBParser::GroupOfAtoms *molecule, const char *cmapFile,
    PDBParser::GroupOfAtoms::GOA_TYPE colorLevel) {
  if (!molecule) {
    return false;
  }
  {
    /*
                    PDBParser::GOAColor* goaColor = new PDBParser::GOAColor();
                    if( cmapFile && strlen(cmapFile))
                    {
                            goaColor->parseColormapFile((char *)cmapFile);
                            goaColor->ApplyColormap(molecule);
                    }
                    else
                            goaColor->ApplyColorsByGOALevel(molecule,
       colorLevel); delete goaColor;
    */
  }
  bool colorBySubunits = false;
  if (colorLevel == PDBParser::GroupOfAtoms::PROTEIN ||
      colorLevel == PDBParser::GroupOfAtoms::COLLECTION) {
    colorBySubunits = true;
  }
  // "flatten" the GOA
  PDBParser::CollectionData *collectionData = 0;
  if (molecule->type == PDBParser::COLLECTION_TYPE) {
    collectionData = molecule->m_CollectionData;
  }
  FlattenGOA(molecule, m_AtomList, collectionData, 0, 0, 0, m_RadiusType,
             PDBParser::ATOM_TYPE, colorBySubunits);
  /*{
          FILE* fp = fopen( "checkColor.txt","w" );

          int i;
          for( i=0; i<m_AtomList.size(); i++ )
          {
                  unsigned char r, g, b;
                  getColors( m_AtomList[i], &r, &g, &b );
                  fprintf( fp, "%s %s %d %d %d\n", m_AtomList[i]->getName(),
                          m_AtomList[i]->getResidueName(),
                          r, g, b );
          }
          fclose( fp );
  }*/
  return true;
}

bool SkinRegion2::getSkinRegion(PDBParser::GroupOfAtoms *molecule, int dim1,
                                int dim2, int dim3, double probeRadius,
                                PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType,
                                SimpleVolumeData *vol, int depth,
                                const char *cmapFile, bool withColor,
                                PDBParser::GroupOfAtoms::GOA_TYPE colorLevel) {
  if (!molecule || !vol) {
    return false;
  }
  m_ProbeRadius = probeRadius;
  m_RadiusType = radiusType;
  m_Depth = depth;
  if (!getColoredAtomList(molecule, cmapFile, colorLevel)) {
    return false;
  }
  // compute the bounding box
  {
    m_Min[0] = m_Min[1] = m_Min[2] = 0.;
    m_Max[0] = m_Max[1] = m_Max[2] = 0.;
    BlurMapsDataManager::getBoundingBox(m_AtomList, m_Min, m_Max, m_RadiusType,
                                        m_ProbeRadius * 2 + 1, 0);
    m_Dim[0] = dim1;
    m_Dim[1] = dim2;
    m_Dim[2] = dim3;
    // origin,
    m_Orig[0] = m_Min[0];
    m_Orig[1] = m_Min[1];
    m_Orig[2] = m_Min[2];
    // and span
    m_Span[0] = (m_Max[0] - m_Min[0]) / (float)(dim1 - 1);
    m_Span[1] = (m_Max[1] - m_Min[1]) / (float)(dim2 - 1);
    m_Span[2] = (m_Max[2] - m_Min[2]) / (float)(dim3 - 1);
  }
  // allocate space
  if (!allocateSpace(withColor, vol)) {
    return false;
  }
  double start_time = getTime();
  m_NumCells = (dim1 - 1) * (dim2 - 1) * (dim3 - 1);
  m_NumVerts = dim1 * dim2 * dim3;
  m_NumSmallCells = pow(pow(2.0f, m_Depth - 1), 3);
  m_NumSmallVerts = pow(pow(2.0f, m_Depth - 1) + 1, 3);

  m_SasInteriorVertices = new bool[m_NumVerts];
  m_SasBoundaryVertices = new bool[m_NumVerts];
  m_VDWInteriorVertices = new bool[m_NumVerts];
  m_VDWBoundaryVertices = new bool[m_NumVerts];
  m_SasBoundaryCells = new bool[m_NumCells];
  m_IndexOfBoundaryCell = new int[m_NumCells];
  m_Regions = new unsigned char[m_NumVerts];
  getSasInterior();
  getBoundaryVerticesCells();
  m_NumBoundaryCells = fillBoundaryIndices();
  m_SphereIndices = new intvec[m_NumBoundaryCells];
  m_SubdividedNodeSasInteriorVertex =
      new bool[m_NumBoundaryCells * m_NumSmallVerts];
  m_SubdividedNodeSasBoundaryCell =
      new bool[m_NumBoundaryCells * m_NumSmallCells];
  getBoundarySphereIndices();
  getSubdividedNodesSasInteriorVertices();
  getSubdividedNodesSasBoundaryCells();
  /*{
          // count how many boundary cells did not have small boundary cells in
  them. int count = 0; int countNonSpheresCells = 0;

          int numSpheres[100];
          {
                  int i;
                  for( i=0; i<100; i++ ) numSpheres[i] = 0;
          }
          {
                  int c;
                  for( c=0; c<m_NumCells; c++ )
                  {
                          if( m_SasBoundaryCells[c] )
                          {
                                  int cellIndex = m_IndexOfBoundaryCell[c];
                                  numSpheres[m_SphereIndices[cellIndex].size()]++;
                          }
                  }
          }
          {
                  printf("Number of cells with given number of spheres:\n");
                  int i;
                  for( i=0; i<100; i++ )
                  {
                          printf("%d spheres: %d cells.\n", i, numSpheres[i] );
                  }
          }
          int c;
          for( c=0; c<m_NumCells; c++ )
          {
                  if( m_SasBoundaryCells[c] )
                  {
                          int cellIndex = m_IndexOfBoundaryCell[c];
                          int numSpheres = m_SphereIndices[cellIndex].size();
                          if( numSpheres == 0 )
                                  countNonSpheresCells++;

                          int sc;
                          bool found = false;
                          for( sc=0; sc<m_NumSmallCells; sc++ )
                          {
                                  if(
  m_SubdividedNodeSasBoundaryCell[cellIndex*m_NumSmallCells + sc] )
                                  {
                                          found = true;
                                          break;
                                  }
                          }
                          if( found ) count++;
                  }
          }

          printf("%d out of %d were true boundary cells\n", count,
  m_NumBoundaryCells ); printf("%d out of %d had no spheres in them\n",
  countNonSpheresCells, m_NumBoundaryCells );
  }*/
  double t1 = getTime() - start_time;
  printf("First part took %lf\n", getTime() - start_time);
  start_time = getTime();
  getSesSurfaceAndVolume();
  double t2 = getTime() - start_time;
  printf("Second part took %lf\n", getTime() - start_time);
  delete[] m_SasInteriorVertices;
  m_SasInteriorVertices = 0;
  delete[] m_SasBoundaryVertices;
  m_SasBoundaryVertices = 0;
  delete[] m_VDWInteriorVertices;
  m_VDWInteriorVertices = 0;
  delete[] m_VDWBoundaryVertices;
  m_VDWBoundaryVertices = 0;
  delete[] m_SasBoundaryCells;
  m_SasBoundaryCells = 0;
  delete[] m_IndexOfBoundaryCell;
  m_IndexOfBoundaryCell = 0;
  delete[] m_SphereIndices;
  m_SphereIndices = 0;
  delete[] m_SubdividedNodeSasInteriorVertex;
  m_SubdividedNodeSasInteriorVertex = 0;
  delete[] m_SubdividedNodeSasBoundaryCell;
  m_SubdividedNodeSasBoundaryCell = 0;
  delete[] m_Regions;
  m_Regions = 0;
  return true;
}
