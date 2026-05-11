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
#ifndef CCV_ATOM_POPULATOR_H
#define CCV_ATOM_POPULATOR_H

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

class AtomPopulator {
public:
  AtomPopulator(PDBParser::GroupOfAtoms *molecule, int dim1, int dim2, int dim3,
                double probeRadius,
                PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType);
  AtomPopulator(PDBParser::GroupOfAtoms *molecule, int dim1, int dim2, int dim3,
                double probeRadius, double floatingBandWidth,
                PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType);
  virtual ~AtomPopulator();
  bool populate(PDBParser::GroupOfAtoms *molecule);
  bool populate(const char *optFileName, bool printMiscFiles = false);

protected:
  bool populateAtomList();
  bool addAtom(int cellx, int celly, int cellz, bool *atomPresent,
               int numSpheres);
  bool printSkin(const char *filename, bool printMiscFiles = false);
  void getBoundingBox(vector<PDBParser::Atom *> &grp, float min[3],
                      float max[3], int radiusType, int probeRadius);
  bool getSasInterior();
  bool getSasBoundary(unsigned char *sasBoundary);
  bool getSasBoundaryCells();
  bool getBoundaryCells();
  inline bool isInterior(int i, int j, int k, unsigned char *interiorArray);
  inline bool isBoundary(int i, int j, int k, unsigned char *interiorArray);
  bool isBoundaryCell(unsigned char *sasCellInterior);
  unsigned char getMaxInArray(unsigned char *sasCellInterior, int n);
  unsigned char *m_SasBoundaryCells;
  unsigned char *m_SasInteriorVertex;
  vector<PDBParser::Atom *> m_AtomList;
  vector<PDBParser::Atom *> m_SkinAtomList;
  float m_Min[3], m_Max[3];
  float m_Orig[3], m_Span[3];
  unsigned int m_Dim[3];
  double m_ProbeRadius;
  double m_floatingBandWidth;
  PDBParser::GroupOfAtoms::RADIUS_TYPE m_RadiusType;
};

inline bool AtomPopulator::isInterior(int i, int j, int k,
                                      unsigned char *interiorArray) {
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
inline bool AtomPopulator::isBoundary(int i, int j, int k,
                                      unsigned char *interiorArray) {
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

#endif
