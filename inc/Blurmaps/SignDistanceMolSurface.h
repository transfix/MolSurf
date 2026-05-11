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
#ifndef CCV_SIGN_DISTANCE_MOL_SURFACE_H
#define CCV_SIGN_DISTANCE_MOL_SURFACE_H

#include <Utility/utility.h>
#include <Blurmaps/BlurMapsDataManager.h>
#include <Blurmaps/GOABlur.h>
#include <PDBParser/Atom.h>
#include <PDBParser/CollectionData.h>
#include <PDBParser/ElementInformation.h>
#include <PDBParser/GroupOfAtoms.h>
#include <PDBParser/ResidueInformation.h>
#include <SimpleVolumeData/SimpleVolumeData.h>
#include <UsefulMath/TrilinearGrid.h>

using CCVOpenGLMath::TrilinearGrid;

namespace PDBParser {
class GroupOfAtoms;
class Atom;
}; // namespace PDBParser

class SimpleVolumeData;

const int NUMBER_OF_SPHERES_PER_VERTEX = 5;
const double _INFINITE = 100000000.;

class SignDistanceMolSurface {
public:
  SignDistanceMolSurface();
  virtual ~SignDistanceMolSurface();

  bool getMolecularSurface(PDBParser::GroupOfAtoms *molecule, int dim1,
                           int dim2, int dim3, double probeRadius,
                           PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType,
                           SimpleVolumeData *vol);

protected:
  bool getsignDistanceGrid();

  void distancePropagation(bool *BoundaryVerticesArray);
  void propagateDistanceInOneSetOfDirections(int i_start, int i_end, int i_inc,
                                             int j_start, int j_end, int j_inc,
                                             int k_start, int k_end, int k_inc,
                                             bool *BoundaryVerticesArray);

  void getBoundaryVertices(bool *BoundaryVerticesArray);
  double computeSignedDistance(int source_vert, int dest_vert);
  void addPrimitive(int srcVertex, int destVertex);
  void assignPrimitives(int vertexIndex, int neighborVertexIndex);
  void assignSign();

  inline bool isBoundary(int i, int j, int k, bool *interiorArray);
  inline void gridToPos(int i, int j, int k, double *x, double *y, double *z);
  inline void gridToPos(int index, double *pos, int dimension);
  inline bool isInterior(int i, int j, int k, bool *interiorArray);

  double signed_dist_pt2sphere(double x_pos, double y_pos, double z_pos,
                               double radius, double c_x, double c_y,
                               double c_z);
  bool isBoundaryCell(bool *sasCellInterior);

  float *m_SignDistanceGrid;

  ////// the sphere index and the vertex index is stored to represent a patch.
  //////////
  int *m_SphereIndices;
  int *m_PrimitivesClosestToVertex;
  ////////////////////////////////////////////////////////////////////////////////////

  bool *m_SasInteriorVertex;
  unsigned char *m_NumberOfAssociatedSpheres;
  unsigned char *m_NumberOfAssociatedPrimitives;

  vector<PDBParser::Atom *> m_AtomList;
  vector<PDBParser::Atom *> m_SkinAtomList;
  float m_Min[3], m_Max[3];
  float m_Orig[3], m_Span[3];
  unsigned int m_Dim[3];
  double m_ProbeRadius;
  PDBParser::GroupOfAtoms::RADIUS_TYPE m_RadiusType;
};

inline bool SignDistanceMolSurface::isInterior(int i, int j, int k,
                                               bool *interiorArray) {
  if (!interiorArray) {
    return false;
  }
  if (i < 0 || i >= m_Dim[0]) {
    return false;
  }
  if (j < 1 || j >= m_Dim[1]) {
    return false;
  }
  if (k < 2 || k >= m_Dim[2]) {
    return false;
  }
  int v = i * m_Dim[0] * m_Dim[0] + j * m_Dim[1] + k;
  if (interiorArray[v]) {
    return true;
  }
  return false;
}

// if current vert not in interior and has a neighbor which is, return true.
inline bool SignDistanceMolSurface::isBoundary(int i, int j, int k,
                                               bool *interiorArray) {
  if (!interiorArray) {
    return false;
  }
  int v = i * m_Dim[0] * m_Dim[0] + j * m_Dim[1] + k;
  if (interiorArray[v]) {
    return false;
  }
  // check all 27 neighbors.
  if (isInterior(i - 1, j - 1, k - 1, interiorArray)) {
    return true;
  }
  if (isInterior(i - 1, j, k - 1, interiorArray)) {
    return true;
  }
  if (isInterior(i - 1, j + 1, k - 1, interiorArray)) {
    return true;
  }
  if (isInterior(i - 1, j - 1, k, interiorArray)) {
    return true;
  }
  if (isInterior(i - 1, j, k, interiorArray)) {
    return true;
  }
  if (isInterior(i - 1, j + 1, k, interiorArray)) {
    return true;
  }
  if (isInterior(i - 1, j - 1, k + 1, interiorArray)) {
    return true;
  }
  if (isInterior(i - 1, j, k + 1, interiorArray)) {
    return true;
  }
  if (isInterior(i - 1, j + 1, k + 1, interiorArray)) {
    return true;
  }
  if (isInterior(i, j - 1, k - 1, interiorArray)) {
    return true;
  }
  if (isInterior(i, j, k - 1, interiorArray)) {
    return true;
  }
  if (isInterior(i, j + 1, k - 1, interiorArray)) {
    return true;
  }
  if (isInterior(i, j - 1, k, interiorArray)) {
    return true;
  }
  if (isInterior(i, j, k, interiorArray)) {
    return true;
  }
  if (isInterior(i, j + 1, k, interiorArray)) {
    return true;
  }
  if (isInterior(i, j - 1, k + 1, interiorArray)) {
    return true;
  }
  if (isInterior(i, j, k + 1, interiorArray)) {
    return true;
  }
  if (isInterior(i, j + 1, k + 1, interiorArray)) {
    return true;
  }
  if (isInterior(i + 1, j - 1, k - 1, interiorArray)) {
    return true;
  }
  if (isInterior(i + 1, j, k - 1, interiorArray)) {
    return true;
  }
  if (isInterior(i + 1, j + 1, k - 1, interiorArray)) {
    return true;
  }
  if (isInterior(i + 1, j - 1, k, interiorArray)) {
    return true;
  }
  if (isInterior(i + 1, j, k, interiorArray)) {
    return true;
  }
  if (isInterior(i + 1, j + 1, k, interiorArray)) {
    return true;
  }
  if (isInterior(i + 1, j - 1, k + 1, interiorArray)) {
    return true;
  }
  if (isInterior(i + 1, j, k + 1, interiorArray)) {
    return true;
  }
  if (isInterior(i + 1, j + 1, k + 1, interiorArray)) {
    return true;
  }
  return false;
}

inline void SignDistanceMolSurface::gridToPos(int i, int j, int k, double *x,
                                              double *y, double *z) {
  *x = m_Orig[0] + m_Span[0] * i;
  *y = m_Orig[1] + m_Span[1] * j;
  *z = m_Orig[2] + m_Span[2] * k;
}

inline void SignDistanceMolSurface::gridToPos(int index, double *pos,
                                              int dimension) {
  *pos = m_Orig[dimension] + m_Span[dimension] * index;
}

#endif
