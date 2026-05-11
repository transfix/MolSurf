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
#include <Blurmaps/SurfaceAtomExtractor.h>
#include <PDBParser/FlattenGOA.h>

using CCVOpenGLMath::TrilinearGrid;

extern double getTime();

SurfaceAtomExtractor::SurfaceAtomExtractor() {
  m_SasInteriorVertices = 0;
  m_SasBoundaryVertices = 0;
  m_SasBoundaryCells = 0;
  m_IndexOfBoundaryCell = 0;
}

SurfaceAtomExtractor::~SurfaceAtomExtractor() {
  if (m_SasInteriorVertices) {
    delete[] m_SasInteriorVertices;
    m_SasInteriorVertices = 0;
  }
  if (m_SasBoundaryVertices) {
    delete[] m_SasBoundaryVertices;
    m_SasBoundaryVertices = 0;
  }
  if (m_SasBoundaryCells) {
    delete[] m_SasBoundaryCells;
    m_SasBoundaryCells = 0;
  }
  if (m_IndexOfBoundaryCell) {
    delete[] m_IndexOfBoundaryCell;
    m_IndexOfBoundaryCell = 0;
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

bool SurfaceAtomExtractor::intersects(double *center, double radius, int iIndex,
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
bool SurfaceAtomExtractor::intersects(double *center, double radius, int cell) {
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

bool SurfaceAtomExtractor::isBoundaryCell(bool *sasCellInterior) {
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

void SurfaceAtomExtractor::getBoundaryVerticesCells() {
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
          m_SasBoundaryVertices[cellVertexArray[j]] = true;
        }
        m_SasBoundaryCells[i] = true;
      }
    }
  }
}

/*!
Given a set of cells, we give an index to each that is a boundary cell.
We also make the indices of other cells as -1.
We return the number of boundary cells.
*/
int SurfaceAtomExtractor::fillBoundaryIndices() {
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
void SurfaceAtomExtractor::getBoundingBox(int x, int y, int z, int *imin,
                                          int *jmin, int *kmin, int *imax,
                                          int *jmax, int *kmax,
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

void SurfaceAtomExtractor::getSphereRadiusCenter(PDBParser::Atom *at,
                                                 double *radius,
                                                 double *center) {
  if (!at || !radius || !center) {
    return;
  }
  *radius = at->getRadius(m_RadiusType);
  center[0] = at->m_Position[0];
  center[1] = at->m_Position[1];
  center[2] = at->m_Position[2];
}

void SurfaceAtomExtractor::getBoundarySphereIndices(bool **boundaryAtom) {
  int numAdded = 0;
  int numBCells = 0;
  int s;
  for (s = 0; s < m_AtomList.size(); s++) {
    (*boundaryAtom)[s] = false;
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
            (*boundaryAtom)[s] = true;
            numAdded++;
          }
        }
      }
    }
  }
  printf("Number added = %d\n", numAdded);
}

void SurfaceAtomExtractor::getSasInterior() {
  {
    int c;
    for (c = 0; c < m_NumVerts; c++) {
      m_SasInteriorVertices[c] = false;
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
          }
        }
      }
    }
    // end blurring kernel
  }
}

bool SurfaceAtomExtractor::getBoundaryAtoms(
    PDBParser::GroupOfAtoms *molecule, int dim1, int dim2, int dim3,
    double probeRadius, PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType,
    bool **boundaryAtom) {
  m_ProbeRadius = probeRadius;
  m_RadiusType = radiusType;
  // "flatten" the GOA
  PDBParser::CollectionData *collectionData = 0;
  if (molecule->type == PDBParser::COLLECTION_TYPE) {
    collectionData = molecule->m_CollectionData;
  }
  FlattenGOA(molecule, m_AtomList, collectionData, 0, 0, 0, radiusType,
             PDBParser::ATOM_TYPE, false);
  (*boundaryAtom) = new bool[m_AtomList.size()];
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
  m_NumCells = (dim1 - 1) * (dim2 - 1) * (dim3 - 1);
  m_NumVerts = dim1 * dim2 * dim3;
  m_SasInteriorVertices = new bool[m_NumVerts];
  m_SasBoundaryVertices = new bool[m_NumVerts];
  m_SasBoundaryCells = new bool[m_NumCells];
  m_IndexOfBoundaryCell = new int[m_NumCells];
  getSasInterior();
  getBoundaryVerticesCells();
  m_NumBoundaryCells = fillBoundaryIndices();
  getBoundarySphereIndices(boundaryAtom);
  delete[] m_SasInteriorVertices;
  m_SasInteriorVertices = 0;
  delete[] m_SasBoundaryVertices;
  m_SasBoundaryVertices = 0;
  delete[] m_SasBoundaryCells;
  m_SasBoundaryCells = 0;
  delete[] m_IndexOfBoundaryCell;
  m_IndexOfBoundaryCell = 0;
  return true;
}
