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
#include <Blurmaps/SignDistanceMolSurface.h>
#include <PDBParser/FlattenGOA.h>

SignDistanceMolSurface::SignDistanceMolSurface() {
  m_SignDistanceGrid = 0;
  m_SphereIndices = 0;
  m_PrimitivesClosestToVertex = 0;
  m_SasInteriorVertex = 0;
  m_NumberOfAssociatedSpheres = 0;
  m_NumberOfAssociatedPrimitives = 0;
  m_ProbeRadius = 1.4;
  m_RadiusType = PDBParser::GroupOfAtoms::VDW_RADIUS;
}

SignDistanceMolSurface::~SignDistanceMolSurface() {
  if (m_SphereIndices) {
    delete[] m_SphereIndices;
    m_SphereIndices = 0;
  }
  if (m_PrimitivesClosestToVertex) {
    delete[] m_PrimitivesClosestToVertex;
    m_PrimitivesClosestToVertex = 0;
  }
  if (m_SasInteriorVertex) {
    delete[] m_SasInteriorVertex;
    m_SasInteriorVertex = 0;
  }
  if (m_NumberOfAssociatedSpheres) {
    delete[] m_NumberOfAssociatedSpheres;
    m_NumberOfAssociatedSpheres = 0;
  }
  if (m_NumberOfAssociatedPrimitives) {
    delete[] m_NumberOfAssociatedPrimitives;
    m_NumberOfAssociatedPrimitives = 0;
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

void SignDistanceMolSurface::assignSign() {
  int c;
  for (c = 0; c < m_Dim[0] * m_Dim[1] * m_Dim[2]; c++) {
    if (m_SasInteriorVertex[c]) {
      if (m_SignDistanceGrid[c] > 0) {
        m_SignDistanceGrid[c] *= -1;
      }
    }
  }
}

bool SignDistanceMolSurface::getsignDistanceGrid() {
  {
    int c;
    for (c = 0; c < m_Dim[0] * m_Dim[1] * m_Dim[2]; c++) {
      m_SignDistanceGrid[c] = _INFINITE;
      m_NumberOfAssociatedSpheres[c] = 0;
      m_NumberOfAssociatedPrimitives[c] = 0;
      {
        int i;
        for (i = 0; i < NUMBER_OF_SPHERES_PER_VERTEX; i++) {
          m_SphereIndices[c * NUMBER_OF_SPHERES_PER_VERTEX * 2 + i + 0] = -1;
          m_SphereIndices[c * NUMBER_OF_SPHERES_PER_VERTEX * 2 + i + 1] = -1;
          m_PrimitivesClosestToVertex[c * NUMBER_OF_SPHERES_PER_VERTEX * 2 + i +
                                      0] = -1;
          m_PrimitivesClosestToVertex[c * NUMBER_OF_SPHERES_PER_VERTEX * 2 + i +
                                      1] = -1;
        }
      }
      m_SasInteriorVertex[c] = false;
    }
  }
  unsigned int zoff = 0;
  int i;
  for (i = 0; i < m_AtomList.size(); i++) {
    PDBParser::Atom *at = m_AtomList[i];
    float radius = at->getRadius(m_RadiusType);
    radius += m_ProbeRadius;
    double c[3], maxRad;
    const float eps = 1.0e-3f;
    unsigned int amax[3], amin[3], j, k, l;
    // calculate the maximum radius for this atom
    maxRad = radius + 1;
    // compute the dataset coordinates of the atom's center
    for (j = 0; j < 3; j++) {
      c[j] = (at->m_Position[j] - m_Orig[j]) / m_Span[j];
      c[j] = ((c[j] - floor(c[j])) >= 0.5) ? ceil(c[j]) : floor(c[j]);
    }
    // then compute the bounding box of the atom (maxRad^3)
    for (j = 0; j < 3; j++) {
      int tempAMin = (int)(c[j] - (maxRad / m_Span[j]) - 1);
      tempAMin = (tempAMin < 0) ? 0 : tempAMin;
      amin[j] = tempAMin;
      int tempAMax = (int)(c[j] + (maxRad / m_Span[j]) + 1);
      tempAMax = (tempAMax > m_Dim[j]) ? m_Dim[j] : tempAMax;
      amax[j] = tempAMax;
    }
    c[0] = at->m_Position[0]; // m_Orig[0] + c[0]*m_Span[0];
    c[1] = at->m_Position[1]; // m_Orig[1] + c[1]*m_Span[1];
    c[2] = at->m_Position[2]; // m_Orig[2] + c[2]*m_Span[2];
    // begin blurring kernel
    for (l = amin[2]; l < amax[2]; l++) {
      for (k = amin[1]; k < amax[1]; k++) {
        for (j = amin[0]; j < amax[0]; j++) {
          int n = j + k * m_Dim[0] + (l - zoff) * m_Dim[0] * m_Dim[1];
          if (n < 0 || n >= (m_Dim[0] * m_Dim[1] * m_Dim[2])) {
            printf("Out of range index");
          }
          float pnt[3];
          pnt[0] = m_Orig[0] + j * m_Span[0];
          pnt[1] = m_Orig[1] + k * m_Span[1];
          pnt[2] = m_Orig[2] + l * m_Span[2];
          double r = (c[0] - pnt[0]) * (c[0] - pnt[0]) +
                     (c[1] - pnt[1]) * (c[1] - pnt[1]) +
                     (c[2] - pnt[2]) * (c[2] - pnt[2]);
          if (r < radius * radius) {
            m_SasInteriorVertex[n] = true;
            if (m_NumberOfAssociatedSpheres[n] < NUMBER_OF_SPHERES_PER_VERTEX) {
              // associate <sphere,vertex> with vertex.
              int vertexIndex = TrilinearGrid::xyz2vtx(j, k, l, m_Dim);
              m_SphereIndices[n * NUMBER_OF_SPHERES_PER_VERTEX * 2 +
                              m_NumberOfAssociatedSpheres[n] + 0] = i;
              m_SphereIndices[n * NUMBER_OF_SPHERES_PER_VERTEX * 2 +
                              m_NumberOfAssociatedSpheres[n] + 1] = vertexIndex;
              m_NumberOfAssociatedSpheres[n]++;
            }
          }
        }
      }
    }
    // end blurring kernel
  }
  bool *BoundaryVerticesArray = new bool[m_Dim[0] * m_Dim[1] * m_Dim[2]];
  {
    int i;
    for (i = 0; i < m_Dim[0] * m_Dim[1] * m_Dim[2]; i++) {
      BoundaryVerticesArray[i] = false;
    }
  }
  // get boundary vertices and also associate the spheres to its primitives list
  getBoundaryVertices(BoundaryVerticesArray);
  {
    // count total number of primitives associated with all vertices.
    int i;
    int *total = new int[NUMBER_OF_SPHERES_PER_VERTEX + 1];
    for (i = 0; i < NUMBER_OF_SPHERES_PER_VERTEX + 1; i++) {
      total[i] = 0;
    }
    for (i = 0; i < m_Dim[0] * m_Dim[1] * m_Dim[2]; i++) {
      total[m_NumberOfAssociatedPrimitives[i]]++;
    }
    for (i = 0; i < NUMBER_OF_SPHERES_PER_VERTEX + 1; i++) {
      printf("%d primitives in %d cases\n", i, total[i]);
    }
    delete[] total;
  }
  for (i = 0; i < m_Dim[0] * m_Dim[1] * m_Dim[2]; i++) {
    if (BoundaryVerticesArray[i] == true) {
      m_SignDistanceGrid[i] = computeSignedDistance(i, i);
    }
  }
  distancePropagation(BoundaryVerticesArray);
  assignSign();
  delete[] BoundaryVerticesArray;
  BoundaryVerticesArray = 0;
  return true;
}

double SignDistanceMolSurface::signed_dist_pt2sphere(double x_pos, double y_pos,
                                                     double z_pos,
                                                     double radius, double c_x,
                                                     double c_y, double c_z) {
  return sqrt((x_pos - c_x) * (x_pos - c_x) + (y_pos - c_y) * (y_pos - c_y) +
              (z_pos - c_z) * (z_pos - c_z)) -
         radius;
}

bool SignDistanceMolSurface::isBoundaryCell(bool *sasCellInterior) {
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

void SignDistanceMolSurface::addPrimitive(int srcVertex, int destVertex) {
  int i, j;
  bool skipFlag;
  /// the two below define a patch
  int currentSphere;
  int currentVertex;
  for (j = 0; j < m_NumberOfAssociatedSpheres[srcVertex]; j++) {
    skipFlag = false;
    currentSphere =
        m_SphereIndices[srcVertex * NUMBER_OF_SPHERES_PER_VERTEX * 2 + j + 0];
    currentVertex =
        m_SphereIndices[srcVertex * NUMBER_OF_SPHERES_PER_VERTEX * 2 + j + 1];
    for (i = 0; i < m_NumberOfAssociatedPrimitives[destVertex]; i++) {
      if ((currentSphere ==
           m_PrimitivesClosestToVertex
               [destVertex * NUMBER_OF_SPHERES_PER_VERTEX * 2 + i + 0]) &&
          (currentVertex ==
           m_PrimitivesClosestToVertex
               [destVertex * NUMBER_OF_SPHERES_PER_VERTEX * 2 + i + 1])) {
        skipFlag = true;
        break;
      }
    }
    if (skipFlag == true) {
      continue;
    } else {
      if (NUMBER_OF_SPHERES_PER_VERTEX ==
          m_NumberOfAssociatedPrimitives[destVertex]) {
        continue;
      }
      m_PrimitivesClosestToVertex[destVertex * NUMBER_OF_SPHERES_PER_VERTEX *
                                      2 +
                                  i + 0] = currentSphere;
      m_PrimitivesClosestToVertex[destVertex * NUMBER_OF_SPHERES_PER_VERTEX *
                                      2 +
                                  i + 1] = currentVertex;
      m_NumberOfAssociatedPrimitives[destVertex]++;
    }
  }
}

double SignDistanceMolSurface::computeSignedDistance(int source_vert,
                                                     int dest_vert) {
  int primitiveIndex = dest_vert * NUMBER_OF_SPHERES_PER_VERTEX;
  int i;
  int x_idx, y_idx, z_idx;
  double x_pos, y_pos, z_pos;
  TrilinearGrid::vtx2xyz(source_vert, x_idx, y_idx, z_idx, m_Dim);
  gridToPos(x_idx, y_idx, z_idx, &x_pos, &y_pos, &z_pos);
  double min_distance, current_distance;
  min_distance = _INFINITE;
  for (i = 0; i < m_NumberOfAssociatedPrimitives[dest_vert]; i++) {
    int atomIndex = m_PrimitivesClosestToVertex[primitiveIndex * 2 + i + 0];
    int vertexIndex = m_PrimitivesClosestToVertex[primitiveIndex * 2 + i + 1];
    double radius = m_AtomList[atomIndex]->getRadius() + m_ProbeRadius;
    double x = m_AtomList[atomIndex]->m_Position[0];
    double y = m_AtomList[atomIndex]->m_Position[1];
    double z = m_AtomList[atomIndex]->m_Position[2];
    current_distance =
        signed_dist_pt2sphere(x_pos, y_pos, z_pos, radius, x, y, z);
    if (i == 0 || current_distance < min_distance) {
      min_distance = current_distance;
    }
  }
  //     assert(i==0);
  return min_distance;
}

void SignDistanceMolSurface::assignPrimitives(int vertexIndex,
                                              int neighborVertexIndex) {
  int i;
  int no = m_NumberOfAssociatedPrimitives[vertexIndex] =
      m_NumberOfAssociatedPrimitives[neighborVertexIndex];
  for (i = 0; i < no; i++) {
    m_PrimitivesClosestToVertex[vertexIndex * NUMBER_OF_SPHERES_PER_VERTEX * 2 +
                                i + 0] =
        m_PrimitivesClosestToVertex[neighborVertexIndex *
                                        NUMBER_OF_SPHERES_PER_VERTEX * 2 +
                                    i + 0];
    m_PrimitivesClosestToVertex[vertexIndex * NUMBER_OF_SPHERES_PER_VERTEX * 2 +
                                i + 1] =
        m_PrimitivesClosestToVertex[neighborVertexIndex *
                                        NUMBER_OF_SPHERES_PER_VERTEX * 2 +
                                    i + 1];
  }
}

void SignDistanceMolSurface::propagateDistanceInOneSetOfDirections(
    int i_start, int i_end, int i_inc, int j_start, int j_end, int j_inc,
    int k_start, int k_end, int k_inc, bool *BoundaryVerticesArray) {
  int i, j, k;
  for (k = k_start; k != k_end; k += k_inc) {
    for (j = j_start; j != j_end; j += j_inc) {
      for (i = i_start; i != i_end; i += i_inc) {
        int vertexIndex = TrilinearGrid::xyz2vtx(i, j, k, m_Dim);
        if (BoundaryVerticesArray[vertexIndex] == true) {
          continue;
        }
        double x_pos, y_pos, z_pos;
        gridToPos(i, j, k, &x_pos, &y_pos, &z_pos);
        int neighborArray[6];
        int numberNeighbor =
            TrilinearGrid::getNeighbor(i, j, k, neighborArray, m_Dim);
        double temp_distance[6];
        for (int neighborIndex = 0; neighborIndex < numberNeighbor;
             neighborIndex++) {
          temp_distance[neighborIndex] =
              computeSignedDistance(vertexIndex, neighborArray[neighborIndex]);
        }
        double min_distance = temp_distance[0];
        int lambda = 0;
        for (int l = 1; l < numberNeighbor; l++) {
          if (temp_distance[l] < min_distance) {
            lambda = l;
            min_distance = temp_distance[l];
          }
        }
        // consider sign
        if (m_SignDistanceGrid[vertexIndex] > min_distance) {
          m_SignDistanceGrid[vertexIndex] = min_distance;
          assignPrimitives(vertexIndex, neighborArray[lambda]);
        }
      }
    }
  }
}

void SignDistanceMolSurface::distancePropagation(bool *BoundaryVerticesArray) {
  int i_start, i_end, i_inc;
  int j_start, j_end, j_inc;
  int k_start, k_end, k_inc;
  // (x+,y+,z+)
  i_start = 0;
  i_end = m_Dim[0];
  i_inc = 1;
  j_start = 0;
  j_end = m_Dim[1];
  j_inc = 1;
  k_start = 0;
  k_end = m_Dim[2];
  k_inc = 1;
  propagateDistanceInOneSetOfDirections(i_start, i_end, i_inc, j_start, j_end,
                                        j_inc, k_start, k_end, k_inc,
                                        BoundaryVerticesArray);
  // (x-,y+,z+)
  i_start = m_Dim[0] - 1;
  i_end = -1;
  i_inc = -1;
  j_start = 0;
  j_end = m_Dim[1];
  j_inc = 1;
  k_start = 0;
  k_end = m_Dim[2];
  k_inc = 1;
  propagateDistanceInOneSetOfDirections(i_start, i_end, i_inc, j_start, j_end,
                                        j_inc, k_start, k_end, k_inc,
                                        BoundaryVerticesArray);
  // (x+,y-,z+)
  i_start = 0;
  i_end = m_Dim[0];
  i_inc = 1;
  j_start = m_Dim[1] - 1;
  j_end = -1;
  j_inc = -1;
  k_start = 0;
  k_end = m_Dim[2];
  k_inc = 1;
  propagateDistanceInOneSetOfDirections(i_start, i_end, i_inc, j_start, j_end,
                                        j_inc, k_start, k_end, k_inc,
                                        BoundaryVerticesArray);
  // (x-,y-,z+)
  i_start = m_Dim[0] - 1;
  i_end = -1;
  i_inc = -1;
  j_start = m_Dim[1] - 1;
  j_end = -1;
  j_inc = -1;
  k_start = 0;
  k_end = m_Dim[2];
  k_inc = 1;
  propagateDistanceInOneSetOfDirections(i_start, i_end, i_inc, j_start, j_end,
                                        j_inc, k_start, k_end, k_inc,
                                        BoundaryVerticesArray);
  // (x+,y+,z-)
  i_start = 0;
  i_end = m_Dim[0];
  i_inc = 1;
  j_start = 0;
  j_end = m_Dim[1];
  j_inc = 1;
  k_start = m_Dim[2] - 1;
  k_end = -1;
  k_inc = -1;
  propagateDistanceInOneSetOfDirections(i_start, i_end, i_inc, j_start, j_end,
                                        j_inc, k_start, k_end, k_inc,
                                        BoundaryVerticesArray);
  // (x-,y+,z-)
  i_start = m_Dim[0] - 1;
  i_end = -1;
  i_inc = -1;
  j_start = 0;
  j_end = m_Dim[1];
  j_inc = 1;
  k_start = m_Dim[2] - 1;
  k_end = -1;
  k_inc = -1;
  propagateDistanceInOneSetOfDirections(i_start, i_end, i_inc, j_start, j_end,
                                        j_inc, k_start, k_end, k_inc,
                                        BoundaryVerticesArray);
  // (x+,y-,z-)
  i_start = 0;
  i_end = m_Dim[0];
  i_inc = 1;
  j_start = m_Dim[1] - 1;
  j_end = -1;
  j_inc = -1;
  k_start = m_Dim[2] - 1;
  k_end = -1;
  k_inc = -1;
  propagateDistanceInOneSetOfDirections(i_start, i_end, i_inc, j_start, j_end,
                                        j_inc, k_start, k_end, k_inc,
                                        BoundaryVerticesArray);
  // (x-,y-,z-)
  i_start = m_Dim[0] - 1;
  i_end = -1;
  i_inc = -1;
  j_start = m_Dim[1] - 1;
  j_end = -1;
  j_inc = -1;
  k_start = m_Dim[2] - 1;
  k_end = -1;
  k_inc = -1;
  propagateDistanceInOneSetOfDirections(i_start, i_end, i_inc, j_start, j_end,
                                        j_inc, k_start, k_end, k_inc,
                                        BoundaryVerticesArray);
}

void SignDistanceMolSurface::getBoundaryVertices(bool *BoundaryVerticesArray) {
  int numberCell, i, j, k;
  numberCell = (m_Dim[0] - 1) * (m_Dim[1] - 1) * (m_Dim[2] - 1);
  int cellVertexArray[8];
  bool sasCellInterior[8];
  for (i = 0; i < numberCell; i++) {
    TrilinearGrid::getCellVertices(i, cellVertexArray, m_Dim);
    for (j = 0; j < 8; j++) {
      sasCellInterior[j] = m_SasInteriorVertex[cellVertexArray[j]];
    }
    if (isBoundaryCell(sasCellInterior)) {
      for (j = 0; j < 8; j++) {
        int srcVertex = cellVertexArray[j];
        BoundaryVerticesArray[srcVertex] = true;
        for (k = 0; k < 8; k++) {
          int destVertex = cellVertexArray[k];
          addPrimitive(srcVertex, destVertex);
        }
      }
    }
  }
}

bool SignDistanceMolSurface::getMolecularSurface(
    PDBParser::GroupOfAtoms *molecule, int dim1, int dim2, int dim3,
    double probeRadius, PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType,
    SimpleVolumeData *vol) {
  if (!molecule || (dim1 < 2) || (dim1 < 2) || (dim1 < 2) ||
      (probeRadius <= 0) || !vol) {
    return false;
  }
  m_ProbeRadius = probeRadius;
  m_RadiusType = radiusType;
  // "flatten" the GOA
  PDBParser::CollectionData *collectionData = 0;
  if (molecule->type == PDBParser::COLLECTION_TYPE) {
    collectionData = molecule->m_CollectionData;
  }
  FlattenGOA(molecule, m_AtomList, collectionData, 0, 0, 0, radiusType,
             PDBParser::ATOM_TYPE, false);
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
  m_SignDistanceGrid = new float[dim1 * dim2 * dim3];
  m_SphereIndices =
      new int[dim1 * dim2 * dim3 * NUMBER_OF_SPHERES_PER_VERTEX * 2];
  m_PrimitivesClosestToVertex =
      new int[dim1 * dim2 * dim3 * NUMBER_OF_SPHERES_PER_VERTEX * 2];
  m_SasInteriorVertex = new bool[dim1 * dim2 * dim3];
  m_NumberOfAssociatedSpheres = new unsigned char[dim1 * dim2 * dim3];
  m_NumberOfAssociatedPrimitives = new unsigned char[dim1 * dim2 * dim3];
  if (!getsignDistanceGrid()) {
    delete[] m_NumberOfAssociatedPrimitives;
    m_NumberOfAssociatedPrimitives = 0;
    delete[] m_SasInteriorVertex;
    m_SasInteriorVertex = 0;
    delete[] m_NumberOfAssociatedSpheres;
    m_NumberOfAssociatedSpheres = 0;
    delete[] m_SphereIndices;
    m_SphereIndices = 0;
    delete[] m_PrimitivesClosestToVertex;
    m_PrimitivesClosestToVertex = 0;
    delete[] m_SignDistanceGrid;
    m_SignDistanceGrid = 0;
    return false;
  }
  delete[] m_NumberOfAssociatedPrimitives;
  m_NumberOfAssociatedPrimitives = 0;
  delete[] m_SasInteriorVertex;
  m_SasInteriorVertex = 0;
  delete[] m_NumberOfAssociatedSpheres;
  m_NumberOfAssociatedSpheres = 0;
  delete[] m_SphereIndices;
  m_SphereIndices = 0;
  delete[] m_PrimitivesClosestToVertex;
  m_PrimitivesClosestToVertex = 0;
  vol->setNumberOfVariables(1);
  vol->setData(0, m_SignDistanceGrid);
  vol->setType(0, SimpleVolumeData::FLOAT);
  vol->setName(0, "SignDistanceMolSurface");
  vol->setMinExtent(m_Min);
  vol->setMaxExtent(m_Max);
  return true;
}
