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
#include <PDBParser/GroupOfAtoms.h>
// For circular definition:
#include <PDBParser/Atom.h>
#include <PDBParser/GOAMolGeometry.h>

using namespace PDBParser;
using namespace CCVOpenGLMath;

GroupOfAtoms::GroupOfAtoms() {
  type = PROTEIN_TYPE;
  resID = -9999;
  resPos = -9999;
  m_ResSeq = 0;
  resHeader = 0;
  startRes = 0;
  endRes = 0;
  multi_r = 0.0;
  multi_g = 0.0;
  multi_b = 0.0;
  // ADD to copy constructer also !
  m_ResidueData = 0;
  m_SecondaryStructureData = 0;
  m_ProteinData = 0;
  m_CollectionData = 0;
  m_Min[0] = 1;
  m_Min[1] = 1;
  m_Min[2] = 1;
  m_Max[0] = 0;
  m_Max[1] = 0;
  m_Max[2] = 0;
}

GroupOfAtoms::GroupOfAtoms(const GroupOfAtoms &goa) {
  int i, n;
  n = goa.m_SubGroups.size();
  for (i = 0; i < n; i++) {
    GroupOfAtoms *subCopy =
        new GroupOfAtoms(*(goa.m_SubGroups.at(i))); // Recursive copy
    m_SubGroups.push_back(subCopy);
  }
  n = goa.m_Atoms.size();
  for (i = 0; i < n; i++) {
    Atom *atomCopy = new Atom(*(goa.m_Atoms.at(i)));
    m_Atoms.push_back(atomCopy);
  }
  n = goa.m_MultiresAtoms.size();
  for (i = 0; i < n; i++) {
    Atom *atomCopy = new Atom(*(goa.m_MultiresAtoms[i]));
    m_MultiresAtoms.push_back(atomCopy);
  }
  int numAtoms = goa.m_Atoms.size();
  type = goa.type;
  resID = goa.resID;
  m_ResSeq = goa.m_ResSeq;
  resPos = goa.resPos;
  if (goa.resHeader) {
    resHeader = new ResidueHeader(*(goa.resHeader));
  } else {
    resHeader = 0;
  }
  if (goa.startRes) {
    startRes = new ResidueHeader(*(goa.startRes));
  } else {
    startRes = 0;
  }
  if (goa.endRes) {
    endRes = new ResidueHeader(*(goa.endRes));
  } else {
    endRes = 0;
  }
  multi_r = goa.multi_r;
  multi_g = goa.multi_g;
  multi_b = goa.multi_b;
  if (goa.m_ResidueData) {
    m_ResidueData = new ResidueData(*(goa.m_ResidueData));
  } else {
    m_ResidueData = 0;
  }
  if (goa.m_SecondaryStructureData) {
    m_SecondaryStructureData =
        new SecondaryStructureData(*(goa.m_SecondaryStructureData));
  } else {
    m_SecondaryStructureData = 0;
  }
  if (goa.m_ProteinData) {
    m_ProteinData = new ProteinData(*(goa.m_ProteinData));
  } else {
    m_ProteinData = 0;
  }
  if (goa.m_CollectionData) {
    m_CollectionData = new CollectionData(*(goa.m_CollectionData));
  } else {
    m_CollectionData = 0;
  }
  m_Min[0] = goa.m_Min[0];
  m_Min[1] = goa.m_Min[1];
  m_Min[2] = goa.m_Min[2];
  m_Max[0] = goa.m_Max[0];
  m_Max[1] = goa.m_Max[1];
  m_Max[2] = goa.m_Max[2];
}

GroupOfAtoms::~GroupOfAtoms() {
  // delete the sub groups
  {
    std::vector<GroupOfAtoms *>::iterator iter = m_SubGroups.begin(),
                                          end = m_SubGroups.end();
    for (; iter != end; ++iter) {
      delete *iter;
    }
    m_SubGroups.clear();
  }
  // delete the atoms
  {
    std::vector<Atom *>::iterator iter = m_Atoms.begin(), end = m_Atoms.end();
    for (; iter != end; ++iter) {
      delete *iter;
    }
    m_Atoms.clear();
  }
  delete resHeader;
  resHeader = 0;
  delete startRes;
  startRes = 0;
  delete endRes;
  endRes = 0;
  delete m_ResidueData;
  m_ResidueData = 0;
  delete m_SecondaryStructureData;
  m_SecondaryStructureData = 0;
  delete m_ProteinData;
  m_ProteinData = 0;
  delete m_CollectionData;
  m_CollectionData = 0;
}

void GroupOfAtoms::addAtom(int serial, int elementIndex, int resNum, char *name,
                           unsigned char resId, char chain, double x, double y,
                           double z, double radius, float charge,
                           bool boundaryAtom) {
  Atom *atom = new Atom();
  atom->setAtomNum(serial);
  atom->setElementIndex(elementIndex);
  atom->setPosition(x, y, z);
  atom->setResidueNum(resNum);
  atom->setChain(chain);
  atom->setBoundaryFlag(boundaryAtom);
  atom->setGOAptr(this);
  if (elementIndex < 0 || elementIndex >= PDBParser::MAX_BIOCHEM_ELEMENTS) {
    atom->setRadius(radius);
    atom->setName(name);
    atom->setResidueId(resId);
    atom->setCharge(charge);
  } else {
    // the info is picked up from a standard table.
  }
  m_Atoms.push_back(atom);
  m_ResSeq = resNum;
}

// Color should be on [0.0, 1.0].  Values outside the range are capped.
void GroupOfAtoms::setGOAColor(double r, double g, double b) {
  if (r < 0.0) {
    r = 0.0;
  }
  if (r > 1.0) {
    r = 1.0;
  }
  if (g < 0.0) {
    g = 0.0;
  }
  if (g > 1.0) {
    g = 1.0;
  }
  if (b < 0.0) {
    b = 0.0;
  }
  if (b > 1.0) {
    b = 1.0;
  }
  // Set the color used for this group of atoms when it is drawn in
  // multiresolution
  multi_r = r;
  multi_g = g;
  multi_b = b;
}

// Based on the residue ID, apply the appropriate color
void GroupOfAtoms::setResidueColor(int arg_resID) {
  double r, g, b;
  if (arg_resID < 0 || arg_resID >= TOTAL_RES_TYPE) {
    // Use default residue for unidentified ID's. It is the 0th residue
    r = resColorR[0];
    g = resColorG[0];
    b = resColorB[0];
  } else {
    // Lookup the color
    r = resColorR[arg_resID];
    g = resColorG[arg_resID];
    b = resColorB[arg_resID];
  }
  setGOAColor(r, g, b);
}

// Based on the secondary structure ID, apply the appropriate color
void GroupOfAtoms::setSSColor(int ssID) {
  double r, g, b;
  double randFactor = rand() / (RAND_MAX + 1.0) / 2.0;
  if (ssID < 0 || ssID >= TOTAL_SS_TYPE) {
    // Use black*(0 - 0.25) as a default color for unidentified ID's
    r = randFactor;
    g = randFactor;
    b = randFactor;
  } else {
    // Lookup the color , multiply by (0.75-1.0) to show separate SSs in the
    // molecule.
    r = ssColorR[ssID] * (randFactor + 0.5);
    g = ssColorG[ssID] * (randFactor + 0.5);
    b = ssColorB[ssID] * (randFactor + 0.5);
  }
  setGOAColor(r, g, b);
}

// Apply the default molecule-level color
void GroupOfAtoms::setMoleculeColor(void) {
  setGOAColor(moleculeColor[0], moleculeColor[1], moleculeColor[2]);
}

// Cycles through all component atoms and subgroups of the current GOA and
// creates a single bounding sphere that is guanranteed to contain all of the
// components.  Recursively creates bounding spheres for all subgroups.
//
// A bounding sphere is first created for all of the subgroups, and then a
// bounding sphere is created for the current group of atoms.
void GroupOfAtoms::createSingleBoundingSphere(double xmid, double ymid,
                                              double zmid) {
  int i;
  int n; // Size of the m_Atoms array
  int m; // Size of the m_SubGroups array
  // Recursively create a bounding sphere for all of the subgroups
  n = (int)m_Atoms.size();
  m = (int)m_SubGroups.size();
  for (i = 0; i < m; i++) {
    GroupOfAtoms *currSubGroup = m_SubGroups.at(i);
    currSubGroup->createSingleBoundingSphere(xmid, ymid, zmid);
  }
  if (m_MultiresAtoms.size() > 1) {
    printf("This is wierd! Residue in more than one SS\n");
  } else {
    // The average of the atom's centers
    float x_avg = 0;
    float y_avg = 0;
    float z_avg = 0;
    float r_max = 0;
    // Tally the centers and radii from the subgroups
    for (i = 0; i < m; i++) {
      // There is only one multires. atom for each subgroup, that is, the single
      // bounding sphere. In other words, for each subgroup, get it bounding
      // sphere, which is stored as the first item in its multires array.
      Atom *atom = m_SubGroups.at(i)->m_MultiresAtoms[0];
      x_avg += atom->m_Position[0];
      y_avg += atom->m_Position[1];
      z_avg += atom->m_Position[2];
      if (atom->getRadius() > r_max) {
        r_max = atom->getRadius();
      }
    }
    // Tally the centers and radii from the atoms
    for (i = 0; i < n; i++) {
      Atom *atom = m_Atoms.at(i);
      x_avg += atom->m_Position[0];
      y_avg += atom->m_Position[1];
      z_avg += atom->m_Position[2];
      if (atom->getRadius() > r_max) {
        r_max = atom->getRadius();
      }
    }
    // Get the averages
    x_avg /= (n + m);
    y_avg /= (n + m);
    z_avg /= (n + m);
    // Cycle through the arrays again, this time looking for the atom
    // furthest away from (x_avg, y_avg, z_avg).  Record the max distance.
    float x, y, z;
    float curr_d; // Actually, contains the square of the current distance, for
                  // speed purposes
    float d_max = 0;
    // Subgroup array
    for (i = 0; i < m; i++) {
      Atom *atom = m_SubGroups.at(i)->m_MultiresAtoms[0];
      // Get the offsets
      x = atom->m_Position[0] - x_avg;
      y = atom->m_Position[1] - y_avg;
      z = atom->m_Position[2] - z_avg;
      // Get the square of the distance to (x_avg, y_avg, z_avg)
      curr_d = (x * x) + (y * y) + (z * z);
      if (curr_d > d_max) {
        d_max = curr_d;
      }
    }
    // Atom array
    for (i = 0; i < n; i++) {
      Atom *atom = m_Atoms.at(i);
      // Get the offsets
      x = atom->m_Position[0] - x_avg;
      y = atom->m_Position[1] - y_avg;
      z = atom->m_Position[2] - z_avg;
      // Get the square of the distance to (x_avg, y_avg, z_avg)
      curr_d = (x * x) + (y * y) + (z * z);
      if (curr_d > d_max) {
        d_max = curr_d;
      }
    }
    // The value in d_max is currently the square of the appropriate distance,
    // so take the square root
    d_max = sqrt(d_max);
    // The center of the bounding sphere will be (x_avg, y_avg, z_avg).  The
    // radius of the bounding sphere will be the distance to the furtherst
    // atomic center plus the largest radius, thus guaranteeing that all
    // component atoms are within the bounding sphere.
    Atom *currAtom;
    currAtom = new Atom();
    currAtom->m_Position[0] = x_avg;
    currAtom->m_Position[1] = y_avg;
    currAtom->m_Position[2] = z_avg;
    currAtom->setRadius(d_max + r_max);
    currAtom->m_ElementIndex = -1; // no longer a normal atom
    m_MultiresAtoms.push_back(currAtom);
    // Currently, there is only one atom in the multiresolution representation
    // for each residue
    int spot = 0;
    int numAtoms = 1;
    float r;
    // Load the information from the new atom that was just created
    x = currAtom->m_Position[0];
    y = currAtom->m_Position[1];
    z = currAtom->m_Position[2];
    r = currAtom->getRadius();
    // No need to set the color.  Color is already set inside the current GOA.
  }
}

// After all of the component atoms have been added to a molecule, call this
// method to initialize the multiresolution representation of subgroups.
void GroupOfAtoms::createMultires(double x_mid, double y_mid, double z_mid,
                                  int multires_type) {
  if (multires_type == MULTIRES_BOUNDING_SPHERE) {
    // Create a hierarchical multiresolution tree using bounding spheres
    createSingleBoundingSphere(x_mid, y_mid, z_mid);
  }
}

bool GroupOfAtoms::getAngles(unsigned int n, double *psi, double *phi,
                             unsigned char *resType) {
  // SKVINAY rewrite this to return appropriate angles from ResidueData
  return false;
}

// This routine currently deals with amino acid residues
// It takes in a molecule.
// For each residue, it rearranges the atoms according to a template
// (residues CA's are moved to the origin and the CA-CB bond is lined up with
// the +z axis)
void GroupOfAtoms::normalizeResidues() {}

bool GroupOfAtoms::getNumberOfAtomsRecursive(int *count) {
  if (!count) {
    return false;
  }
  *count += (int)m_Atoms.size();
  int i, m;
  m = (int)m_SubGroups.size();
  for (i = 0; i < m; i++) {
    GroupOfAtoms *currSubGroup = m_SubGroups.at(i);
    currSubGroup->getNumberOfAtomsRecursive(count);
  }
  return true;
}

// if residue and is within the required strand only!
bool GroupOfAtoms::getNumberOfAtoms(int *count, int strandNumber) {
  // this should be called only for Sheets !
  if (type != SS_TYPE) {
    return false;
  }
  if (resID != SS_SHEET) {
    return false;
  }
  if (m_SecondaryStructureData == 0) {
    return false;
  }
  // for given strand, fill up atoms
  int startResId = m_SecondaryStructureData
                       ->m_ResSequenceNumberPairList[strandNumber * 2 + 0];
  int endResId = m_SecondaryStructureData
                     ->m_ResSequenceNumberPairList[strandNumber * 2 + 1];
  int i;
  int m;
  m = (int)m_SubGroups.size();
  for (i = 0; i < m; i++) {
    GroupOfAtoms *currSubGroup = m_SubGroups.at(i);
    if ((currSubGroup->m_ResSeq >= startResId) &&
        (currSubGroup->m_ResSeq <= endResId)) {
      *count += (int)currSubGroup->m_Atoms.size();
    }
  }
  return true;
}

void GroupOfAtoms::fillUpAtoms(double *x, double *y, double *z, int *index) {
  int i;
  int n = (int)m_Atoms.size();
  for (i = 0; i < n; i++) {
    Atom *currAtom = m_Atoms.at(i);
    x[*index] = currAtom->m_Position[0];
    y[*index] = currAtom->m_Position[1];
    z[*index] = currAtom->m_Position[2];
    *index += 1;
  }
  int m;
  m = (int)m_SubGroups.size();
  for (i = 0; i < m; i++) {
    GroupOfAtoms *currSubGroup = m_SubGroups.at(i);
    currSubGroup->fillUpAtoms(x, y, z, index);
  }
}

bool GroupOfAtoms::fillUpAtoms(double *x, double *y, double *z, int *index,
                               int strandNumber) {
  // this should be called only for Sheets !
  if (type != SS_TYPE) {
    return false;
  }
  if (resID != SS_SHEET) {
    return false;
  }
  if (m_SecondaryStructureData == 0) {
    return false;
  }
  // for given strand, fill up atoms
  int startResId = m_SecondaryStructureData
                       ->m_ResSequenceNumberPairList[strandNumber * 2 + 0];
  int endResId = m_SecondaryStructureData
                     ->m_ResSequenceNumberPairList[strandNumber * 2 + 1];
  int m;
  m = (int)m_SubGroups.size();
  int i;
  for (i = 0; i < m; i++) {
    GroupOfAtoms *currSubGroup = m_SubGroups.at(i);
    if ((currSubGroup->m_ResSeq >= startResId) &&
        (currSubGroup->m_ResSeq <= endResId)) {
      int n = (int)currSubGroup->m_Atoms.size();
      int j;
      for (j = 0; j < n; j++) {
        Atom *currAtom = currSubGroup->m_Atoms.at(j);
        x[*index] = currAtom->m_Position[0];
        y[*index] = currAtom->m_Position[1];
        z[*index] = currAtom->m_Position[2];
        *index += 1;
      }
    }
  }
  return true;
}

void GroupOfAtoms::transformAtoms(double *transformation) {
  int i;
  int n = (int)m_Atoms.size();
  for (i = 0; i < n; i++) {
    Atom *currAtom = m_Atoms.at(i);
    double x = currAtom->m_Position[0];
    double y = currAtom->m_Position[1];
    double z = currAtom->m_Position[2];
    double x1 =
        x * transformation[0] + y * transformation[3] + z * transformation[6];
    double y1 =
        x * transformation[1] + y * transformation[4] + z * transformation[7];
    double z1 =
        x * transformation[2] + y * transformation[5] + z * transformation[8];
    currAtom->m_Position[0] = x1;
    currAtom->m_Position[1] = y1;
    currAtom->m_Position[2] = z1;
  }
  n = (int)m_MultiresAtoms.size();
  for (i = 0; i < n; i++) {
    Atom *currAtom = m_MultiresAtoms[i];
    double x = currAtom->m_Position[0];
    double y = currAtom->m_Position[1];
    double z = currAtom->m_Position[2];
    double x1 =
        x * transformation[0] + y * transformation[3] + z * transformation[6];
    double y1 =
        x * transformation[1] + y * transformation[4] + z * transformation[7];
    double z1 =
        x * transformation[2] + y * transformation[5] + z * transformation[8];
    currAtom->m_Position[0] = x1;
    currAtom->m_Position[1] = y1;
    currAtom->m_Position[2] = z1;
  }
  int m;
  m = (int)m_SubGroups.size();
  for (i = 0; i < m; i++) {
    GroupOfAtoms *currSubGroup = m_SubGroups.at(i);
    currSubGroup->transformAtoms(transformation);
  }
  // now the positions array in the GOA, which need to be removed someday, are
  // wrong! SKVINAY
}

bool GroupOfAtoms::addChain(GroupOfAtoms *chain) {
  if (type == PROTEIN_TYPE) {
    m_SubGroups.push_back(chain);
    return true;
  }
  int m;
  m = (int)m_SubGroups.size();
  if (m > 0) {
    GroupOfAtoms *currSubGroup = m_SubGroups.at(0);
    return currSubGroup->addChain(chain);
  } else {
    return false;
  }
}

//---------------------------------------------------------------------------------------
// CHA: Code for build hierarchical information, from very old TexMol
// (~wzhao/Software/TexMol) If there is any problem with this part, first take a
// look the old version.
//---------------------------------------------------------------------------------------

// return the first N atom and the last residues C atom positions
bool GroupOfAtoms::getEndPoints(CCVOpenGLMath::Vector *v1,
                                CCVOpenGLMath::Vector *v2) {
  int m;
  m = (int)m_SubGroups.size();

  {
    // get the first point
    GroupOfAtoms *currSubGroup = m_SubGroups.at(0);
    if (currSubGroup->type != RESIDUE_TYPE)
      return false;

    int n = (int)currSubGroup->m_Atoms.size();
    if (n < 3)
      return false;

    double x, y, z;
    Atom *currAtom = currSubGroup->m_Atoms.at(0);
    x = currAtom->m_Position[0];
    y = currAtom->m_Position[1];
    z = currAtom->m_Position[2];
    v1->set(x, y, z, 1.0);
  }

  {
    // get the last point
    GroupOfAtoms *currSubGroup = m_SubGroups.at(m - 1);
    if (currSubGroup->type != RESIDUE_TYPE)
      return false;

    int n = (int)currSubGroup->m_Atoms.size();
    if (n < 3)
      return false;

    double x, y, z;
    Atom *currAtom = currSubGroup->m_Atoms.at(2);
    x = currAtom->m_Position[0];
    y = currAtom->m_Position[1];
    z = currAtom->m_Position[2];
    v2->set(x, y, z, 1.0);
  }

  return true;
}

// return the first N atom and the last residues C atom positions in the residue
// sequence given
bool GroupOfAtoms::getEndPoints(CCVOpenGLMath::Vector *v1,
                                CCVOpenGLMath::Vector *v2, int startingResidue,
                                int endingResidue) {
  int m;
  m = (int)m_SubGroups.size();

  if (startingResidue >= m)
    return false;
  if (endingResidue >= m)
    return false;
  if (startingResidue > endingResidue)
    return false;

  {
    // get the first point
    GroupOfAtoms *currSubGroup = m_SubGroups.at(startingResidue);
    if (currSubGroup->type != RESIDUE_TYPE)
      return false;

    int n = (int)currSubGroup->m_Atoms.size();
    if (n < 3)
      return false;

    double x, y, z;
    Atom *currAtom = currSubGroup->m_Atoms.at(0);
    x = currAtom->m_Position[0];
    y = currAtom->m_Position[1];
    z = currAtom->m_Position[2];
    v1->set(x, y, z, 1.0);
  }

  {
    // get the last point
    GroupOfAtoms *currSubGroup = m_SubGroups.at(endingResidue);
    if (currSubGroup->type != RESIDUE_TYPE)
      return false;

    int n = (int)currSubGroup->m_Atoms.size();
    if (n < 3)
      return false;

    double x, y, z;
    Atom *currAtom = currSubGroup->m_Atoms.at(2);
    x = currAtom->m_Position[0];
    y = currAtom->m_Position[1];
    z = currAtom->m_Position[2];
    v2->set(x, y, z, 1.0);
  }

  return true;
}

bool GroupOfAtoms::getEndPoints(CCVOpenGLMath::Vector *v1,
                                CCVOpenGLMath::Vector *v2, int strandNumber) {
  if (type != SS_TYPE)
    return false;
  if (resID != SS_SHEET)
    return false;

  int m;
  m = (int)m_SubGroups.size();

  if (m_SecondaryStructureData == 0)
    return false;

  int startResId = m_SecondaryStructureData
                       ->m_ResSequenceNumberPairList[strandNumber * 2 + 0];
  int endResId = m_SecondaryStructureData
                     ->m_ResSequenceNumberPairList[strandNumber * 2 + 1];

  // get the first point
  bool found = false;
  int i;
  for (i = 0; i < m; i++) {
    GroupOfAtoms *currSubGroup = m_SubGroups.at(i);
    if (currSubGroup->m_ResSeq == startResId) {
      int n = (int)currSubGroup->m_Atoms.size();
      if (n < 3)
        return false;

      Atom *currAtom = currSubGroup->m_Atoms.at(0);
      double x, y, z;
      x = currAtom->m_Position[0];
      y = currAtom->m_Position[1];
      z = currAtom->m_Position[2];
      v1->set(x, y, z, 1.0);
      found = true;
      break;
    }
  }
  if (!found)
    return false;

  // get the last point
  found = false;
  for (i = 0; i < m; i++) {
    GroupOfAtoms *currSubGroup = m_SubGroups.at(i);
    if (currSubGroup->m_ResSeq == endResId) {
      int n = (int)currSubGroup->m_Atoms.size();
      if (n < 3)
        return false;

      Atom *currAtom = currSubGroup->m_Atoms.at(2);
      double x, y, z;
      x = currAtom->m_Position[0];
      y = currAtom->m_Position[1];
      z = currAtom->m_Position[2];
      v2->set(x, y, z, 1.0);
      found = true;
      break;
    }
  }
  if (!found)
    return false;

  return true;
}

bool GroupOfAtoms::createHelicesForResidues(
    std::vector<double> *m_HelixGeometry, int startingResidue, int length) {
  if (type != SS_TYPE)
    return false;

  //     In one pass, count atoms
  int endingResidue = startingResidue + length;
  if (endingResidue >= m_SubGroups.size())
    endingResidue = m_SubGroups.size() - 1;

  int count = 0;
  int i;
  for (i = startingResidue; i < m_SubGroups.size(); i++) {
    if (!(m_SubGroups.at(i))->getNumberOfAtomsRecursive(&count))
      return false;
  }
  if (count < 5)
    return false;

  //     create arrays x, y, z
  double *x = new double[count];
  double *y = new double[count];
  double *z = new double[count];

  //     Next pass, fill up x, y, z
  int index = 0;
  for (i = startingResidue; i < m_SubGroups.size(); i++) {
    (m_SubGroups.at(i))->fillUpAtoms(x, y, z, &index);
  }

  LinearAlgebra linearAlgebra;
  CCVOpenGLMath::Vector endPoint1;
  CCVOpenGLMath::Vector endPoint2;
  double radius;

  linearAlgebra.getCylinderFit(count, x, y, z, &endPoint1, &endPoint1, &radius);

  getEndPoints(&endPoint1, &endPoint2, startingResidue, endingResidue);

  m_HelixGeometry->push_back(endPoint1[0]);
  m_HelixGeometry->push_back(endPoint1[1]);
  m_HelixGeometry->push_back(endPoint1[2]);
  m_HelixGeometry->push_back(endPoint2[0]);
  m_HelixGeometry->push_back(endPoint2[1]);
  m_HelixGeometry->push_back(endPoint2[2]);
  m_HelixGeometry->push_back(radius);

  return true;
}

void GroupOfAtoms::createSS() {
  if (type != SS_TYPE)
    return;

  if (this->resID == SS_HELIX) {
    if (m_SecondaryStructureData == 0)
      return; // this is an error in the parsing!

    int i;
    int length = 3;

    for (i = 0; i < m_SubGroups.size(); i += length) {
      createHelicesForResidues(&(m_SecondaryStructureData->m_HelixGeometry), i,
                               length);
    }
  } else {
    if (m_SecondaryStructureData == 0)
      return; // this is an error in the parsing!
    //    linearAlgebra.getCylinderFit( 4, x, y, z, &m_HelixEndPoint1,
    //    &m_HelixEndPoint2, &m_HelixRadius );
    // get number of strands
    // for each strand
    // get number of atoms
    // fill up atoms
    // create cylinder
    //    In one pass, count atoms
    int numStrands =
        m_SecondaryStructureData->m_ResSequenceNumberPairList.size() / 2;
    int i;
    for (i = 0; i < numStrands; i++) {
      int numAtoms = 0;
      if (!getNumberOfAtoms(
              &numAtoms,
              i)) // overloaded to show this is for strands, not for helix
        continue;

      if (numAtoms < 5)
        continue; // why not 5 !

      //   create arrays x, y, z
      double *x = new double[numAtoms];
      double *y = new double[numAtoms];
      double *z = new double[numAtoms];

      //   Next pass, fill up x, y, z
      int index = 0;
      if (!fillUpAtoms(x, y, z, &index,
                       i)) // overloaded, passing in strand index
        continue;

      Vector p1, p2;
      double radius;
      LinearAlgebra linearAlgebra;
      linearAlgebra.getCylinderFit(numAtoms, x, y, z, &p1, &p2, &radius);

      getEndPoints(&p1, &p2, i); // should make this as return false !
      // add values to array of strands geometry
      m_SecondaryStructureData->m_StrandGeometry.push_back(p1[0]);
      m_SecondaryStructureData->m_StrandGeometry.push_back(p1[1]);
      m_SecondaryStructureData->m_StrandGeometry.push_back(p1[2]);
      m_SecondaryStructureData->m_StrandGeometry.push_back(p2[0]);
      m_SecondaryStructureData->m_StrandGeometry.push_back(p2[1]);
      m_SecondaryStructureData->m_StrandGeometry.push_back(p2[2]);
      m_SecondaryStructureData->m_StrandGeometry.push_back(radius);
    }
  }
}

void GroupOfAtoms::createSSs() {
  if (this->type == SS_TYPE) {
    // the SS data should be created by now
    if (m_SecondaryStructureData == 0)
      return;

    if ((this->resID == SS_HELIX) || (this->resID == SS_SHEET))
      createSS();
  }

  int i, m;
  m = (int)m_SubGroups.size();
  for (i = 0; i < m; i++) {
    GroupOfAtoms *currSubGroup = m_SubGroups.at(i);
    currSubGroup->createSSs();
  }
}

//---------------------------------------------------------------------------------------
void GroupOfAtoms::buildHierarchicalInformation(double xmid, double ymid,
                                                double zmid) {
  createMultires(xmid, ymid, zmid, MULTIRES_BOUNDING_SPHERE);
  createSSs();
}

void GroupOfAtoms::getMinMaxRecursive(double *minX, double *minY, double *minZ,
                                      double *maxX, double *maxY,
                                      double *maxZ) {
  if (!minX || !minY || !minZ || !maxX || !maxY || !maxZ) {
    return;
  }
  int m, i;
  m = (int)m_SubGroups.size();
  for (i = 0; i < m; i++) {
    GroupOfAtoms *currSubGroup = m_SubGroups.at(i);
    if (currSubGroup) {
      double subgroup_minX = 1e10, subgroup_minY = 1e10, subgroup_minZ = 1e10,
             subgroup_maxX = -1e10, subgroup_maxY = -1e10,
             subgroup_maxZ = -1e10;
      currSubGroup->getMinMaxRecursive(&subgroup_minX, &subgroup_minY,
                                       &subgroup_minZ, &subgroup_maxX,
                                       &subgroup_maxY, &subgroup_maxZ);
      if (subgroup_minX < *minX) {
        *minX = subgroup_minX;
      }
      if (subgroup_minY < *minY) {
        *minY = subgroup_minY;
      }
      if (subgroup_minZ < *minZ) {
        *minZ = subgroup_minZ;
      }
      if (subgroup_maxX > *maxX) {
        *maxX = subgroup_maxX;
      }
      if (subgroup_maxY > *maxY) {
        *maxY = subgroup_maxY;
      }
      if (subgroup_maxZ > *maxZ) {
        *maxZ = subgroup_maxZ;
      }
    }
  }
  int n = (int)m_Atoms.size();
  for (i = 0; i < n; i++) {
    Atom *currAtom = m_Atoms.at(i);
    if (currAtom) {
      double subgroup_minX = 1e10, subgroup_minY = 1e10, subgroup_minZ = 1e10,
             subgroup_maxX = -1e10, subgroup_maxY = -1e10,
             subgroup_maxZ = -1e10;
      currAtom->getMinMax(&subgroup_minX, &subgroup_minY, &subgroup_minZ,
                          &subgroup_maxX, &subgroup_maxY, &subgroup_maxZ);
      if (subgroup_minX < *minX) {
        *minX = subgroup_minX;
      }
      if (subgroup_minY < *minY) {
        *minY = subgroup_minY;
      }
      if (subgroup_minZ < *minZ) {
        *minZ = subgroup_minZ;
      }
      if (subgroup_maxX > *maxX) {
        *maxX = subgroup_maxX;
      }
      if (subgroup_maxY > *maxY) {
        *maxY = subgroup_maxY;
      }
      if (subgroup_maxZ > *maxZ) {
        *maxZ = subgroup_maxZ;
      }
    }
  }
}

void GroupOfAtoms::getMinMax(double *minX, double *minY, double *minZ,
                             double *maxX, double *maxY, double *maxZ) {
  if (!minX || !minY || !minZ || !maxX || !maxY || !maxZ) {
    return;
  }
  *minX = 1e10, *minY = 1e10, *minZ = 1e10, *maxX = -1e10, *maxY = -1e10,
  *maxZ = -1e10;
  if (m_Min[0] <= m_Max[0] && m_Min[1] <= m_Max[1] && m_Min[2] <= m_Max[2]) {
    *minX = m_Min[0];
    *minY = m_Min[1];
    *minZ = m_Min[2];
    *maxX = m_Max[0];
    *maxY = m_Max[1];
    *maxZ = m_Max[2];
    return;
  }
  getMinMaxRecursive(minX, minY, minZ, maxX, maxY, maxZ);
  m_Min[0] = *minX;
  m_Min[1] = *minY;
  m_Min[2] = *minZ;
  m_Max[0] = *maxX;
  m_Max[1] = *maxY;
  m_Max[2] = *maxZ;
}

void GroupOfAtoms::getMaximumDistanceFromPoint(double xOrigin, double yOrigin,
                                               double zOrigin,
                                               double *distance) {
  if (!distance) {
    return;
  }
  int m, i;
  m = (int)m_SubGroups.size();
  for (i = 0; i < m; i++) {
    GroupOfAtoms *currSubGroup = m_SubGroups.at(i);
    if (currSubGroup) {
      currSubGroup->getMaximumDistanceFromPoint(xOrigin, yOrigin, zOrigin,
                                                distance);
    }
  }
  int n = (int)m_Atoms.size();
  for (i = 0; i < n; i++) {
    Atom *currAtom = m_Atoms.at(i);
    if (currAtom) {
      currAtom->getMaximumDistanceFromPoint(xOrigin, yOrigin, zOrigin,
                                            distance);
    }
  }
}

bool GroupOfAtoms::intersects(CCVOpenGLMath::Ray targetVector) {
  // if there was no multires computed, we have to check the children //
  if (m_MultiresAtoms.size() < 1) {
    return true;
  }
  // intersect the multires enclosing sphere with the ray //
  Vector center(m_MultiresAtoms[0]->m_Position[0],
                m_MultiresAtoms[0]->m_Position[1],
                m_MultiresAtoms[0]->m_Position[2], 1);
  float radius = m_MultiresAtoms[0]->getRadius();
  Vector intersectionPoint1;
  Vector intersectionPoint2;
  float distance1;
  float distance2;
  if (targetVector.intersectSphere(center, radius, &intersectionPoint1,
                                   &intersectionPoint2, &distance1,
                                   &distance2)) {
    if (distance1 > 0 || distance2 > 0) {
      return true;
    }
  }
  return false;
}

void GroupOfAtoms::recursiveGetClosestAtom(Atom **atom, double *distance,
                                           CCVOpenGLMath::Ray targetVector) {
  int m, i;
  m = (int)m_SubGroups.size();
  for (i = 0; i < m; i++) {
    GroupOfAtoms *currSubGroup = m_SubGroups.at(i);
    if (currSubGroup) {
      if (currSubGroup->intersects(targetVector)) {
        currSubGroup->recursiveGetClosestAtom(atom, distance, targetVector);
      }
    }
  }
  int n = (int)m_Atoms.size();
  for (i = 0; i < n; i++) {
    Atom *cur_atom = m_Atoms.at(i);
    if (cur_atom) {
      double cur_distance = cur_atom->getDistanceFromRay(
          targetVector); // returns -1 when something's wrong.
      if (cur_distance < 0) {
        continue;
      }
      // either its the first time, or we are updating
      if ((cur_distance < *distance) || (*distance < 0)) {
        (*atom) = cur_atom;
        // printf("Got %s\n", (*atom)->getName());
        *distance = cur_distance;
      }
    }
  }
}

Atom *GroupOfAtoms::getSelectedAtom(CCVOpenGLMath::Ray targetVector) {
  // recursively search to see which atom is closest to the ray ( in the forward
  // direction ) At each level, we use the multires information.
  Atom *atom = 0;
  double distance = -1;
  recursiveGetClosestAtom(&atom, &distance, targetVector);
  return atom;
}

bool GroupOfAtoms::applyColorMap(const char *fileName) {
  if (!fileName) {
    return false;
  }
  return true;
}

GroupOfAtoms *GroupOfAtoms::getEmptyHierarchy() {
  PDBParser::GroupOfAtoms *temp = new PDBParser::GroupOfAtoms();
  temp->type = PDBParser::COLLECTION_TYPE;
  temp->m_CollectionData = new PDBParser::CollectionData();
  PDBParser::GroupOfAtoms *temp_molecule = new PDBParser::GroupOfAtoms();
  temp_molecule->type = PDBParser::PROTEIN_TYPE;
  temp->m_SubGroups.push_back(temp_molecule);
  PDBParser::GroupOfAtoms *temp_chain = new PDBParser::GroupOfAtoms();
  temp_chain->type = PDBParser::CHAIN_TYPE;
  temp_molecule->m_SubGroups.push_back(temp_chain);
  PDBParser::GroupOfAtoms *temp_ss = new PDBParser::GroupOfAtoms();
  temp_ss->type = PDBParser::SS_TYPE;
  temp_chain->m_SubGroups.push_back(temp_ss);
  PDBParser::GroupOfAtoms *temp_residue = new PDBParser::GroupOfAtoms();
  temp_residue->type = PDBParser::RESIDUE_TYPE;
  temp_ss->m_SubGroups.push_back(temp_residue);
  return temp;
}

PDBParser::GroupOfAtoms::TYPE_OF_RESIDUE GroupOfAtoms::getResidueType() const {
  if (resID == 0) {
    return PDBParser::GroupOfAtoms::UNKNOWN_RESIDUE;
  }
  if (resID >= PDBParser::RES_TYPE_FIRST_NUCLEOTIDE) {
    return PDBParser::GroupOfAtoms::RNA_RESIDUE;
  }
  return PDBParser::GroupOfAtoms::AMINO_ACID_RESIDUE;
}

// if this is a residue, set the required atom whose atom name is given, if more
// than one match, set first one.
Atom *GroupOfAtoms::getAtom(const char *atomName) {
  if (type != PDBParser::RESIDUE_TYPE) {
    return 0;
  }
  int i;
  for (i = 0; i < m_Atoms.size(); i++) {
    PDBParser::Atom *curAtom = m_Atoms.at(i);
    if (strcmp(curAtom->getName(), atomName) == 0) {
      return curAtom;
    }
  }
  return 0;
}

bool GroupOfAtoms::add(std::vector<double> *list, double x, double y, double z,
                       double r) {
  if (!list) {
    return false;
  }
  list->push_back(x);
  list->push_back(y);
  list->push_back(z);
  list->push_back(r);
  return true;
}

bool GroupOfAtoms::intersects(double *min, double *max) {
  if (!min || !max) {
    return false;
  }
  if (m_Min[0] > max[0] || m_Min[1] > max[1] || m_Min[2] > max[2]) {
    return false;
  }
  if (min[0] > m_Max[0] || min[1] > m_Max[1] || min[2] > m_Max[2]) {
    return false;
  }
  return true;
}

// checks if current GOA is entirely within min max
bool GroupOfAtoms::within(double *min, double *max) {
  if (!min || !max) {
    return false;
  }
  if (m_Min[0] < min[0] || m_Min[1] < min[1] || m_Min[2] < min[2]) {
    return false;
  }
  if (m_Max[0] > max[0] || m_Max[1] > max[1] || m_Max[2] > max[2]) {
    return false;
  }
  return true;
}

// return atoms inside / outside / inside and outside this bounding box.
void GroupOfAtoms::getAtoms(double *min, double *max,
                            std::vector<double> *interiorAtomList,
                            std::vector<double> *exteriorAtomList,
                            POSITION position) {
  if (min[0] > max[0]) // my way of seeing if its initialized SKVINAY
  {
    double minX, minY, minZ, maxX, maxY, maxZ;
    getMinMax(&minX, &minY, &minZ, &maxX, &maxY, &maxZ);
  }
  if (!min || !max || (position == INSIDE && !interiorAtomList) ||
      (position == OUTSIDE && !exteriorAtomList) ||
      (position == BOTH && (!interiorAtomList || !exteriorAtomList))) {
    return;
  } else {
    getAtomsRecursive(min, max, interiorAtomList, exteriorAtomList, position);
  }
}

void GroupOfAtoms::getAtomsRecursive(double *min, double *max,
                                     std::vector<double> *interiorAtomList,
                                     std::vector<double> *exteriorAtomList,
                                     POSITION position) {
  // check if this GOA is going to contribute or not //
  {
    bool valid = false;
    if (position == BOTH) {
      valid = true;
    } else if ((position == INSIDE) && intersects(min, max)) {
      valid = true;
    } else if ((position == OUTSIDE) && !within(min, max)) {
      valid = true;
    }
    if (!valid) {
      return;
    }
  }
  int m, i;
  m = (int)m_SubGroups.size();
  for (i = 0; i < m; i++) {
    GroupOfAtoms *currSubGroup = m_SubGroups.at(i);
    if (currSubGroup) {
      currSubGroup->getAtomsRecursive(min, max, interiorAtomList,
                                      exteriorAtomList, position);
    }
  }
  int n = (int)m_Atoms.size();
  for (i = 0; i < n; i++) {
    Atom *cur_atom = m_Atoms.at(i);
    double x = cur_atom->m_Position[0];
    double y = cur_atom->m_Position[1];
    double z = cur_atom->m_Position[2];
    double r = cur_atom->getRadius();
    if (cur_atom) {
      if (position == BOTH) {
        if (x >= min[0] && x <= max[0] && y >= min[1] && y <= max[1] &&
            z >= min[2] && z <= max[2]) {
          add(interiorAtomList, x, y, z, r);
        } else {
          add(exteriorAtomList, x, y, z, r);
        }
      } else if ((position == INSIDE)) {
        if (x >= min[0] && x <= max[0] && y >= min[1] && y <= max[1] &&
            z >= min[2] && z <= max[2]) {
          add(interiorAtomList, x, y, z, r);
        }
      } else if ((position == OUTSIDE)) {
        if (!(x >= min[0] && x <= max[0] && y >= min[1] && y <= max[1] &&
              z >= min[2] && z <= max[2])) {
          add(exteriorAtomList, x, y, z, r);
        }
      }
    }
  }
  return;
}

void GroupOfAtoms::getAtomsRecursive(double *posRadius, int maxLength,
                                     int *curPosition) {
  if ((*curPosition) >= maxLength) {
    return;
  }
  int m, i;
  m = (int)m_SubGroups.size();
  for (i = 0; i < m; i++) {
    GroupOfAtoms *currSubGroup = m_SubGroups.at(i);
    if (currSubGroup) {
      currSubGroup->getAtomsRecursive(posRadius, maxLength, curPosition);
    }
  }
  int n = (int)m_Atoms.size();
  for (i = 0; i < n; i++) {
    Atom *cur_atom = m_Atoms.at(i);
    posRadius[(*curPosition)] = cur_atom->m_Position[0];
    (*curPosition)++;
    posRadius[(*curPosition)] = cur_atom->m_Position[1];
    (*curPosition)++;
    posRadius[(*curPosition)] = cur_atom->m_Position[2];
    (*curPosition)++;
    posRadius[(*curPosition)] = cur_atom->getRadius();
    (*curPosition)++;
  }
}

bool GroupOfAtoms::intToFunctionType(FUNCTIONS *functionType,
                                     int intFunctionType) {
  if (intFunctionType >= NUMBER_OF_FUNCTIONS) {
    *functionType = ELECTRON_DENSITY;
    return false;
  }
  switch (intFunctionType) {
  case 0:
    *functionType = ELECTRON_DENSITY;
    return true;
  case 1:
    *functionType = PER_ATOM_HYDROPHOBICITY;
    return true;
  case 2:
    *functionType = PER_RESIDUE_HYDROPHOBICITY;
    return true;
  case 3:
    *functionType = ELECTROSTATIC_POTENTIAL;
  default:
    *functionType = ELECTRON_DENSITY;
    return false;
  }
}

bool GroupOfAtoms::intToRadiusType(RADIUS_TYPE *radiusType, int intRadiusType) {
  if (intRadiusType >= NUMBER_OF_RADII_TYPES) {
    *radiusType = VDW_RADIUS;
    return false;
  }
  switch (intRadiusType) {
  case 0:
    *radiusType = VDW_RADIUS;
    return true;
  case 1:
    *radiusType = ROUX_RADIUS;
    return true;
  default:
    *radiusType = VDW_RADIUS;
    return false;
  }
}

bool GroupOfAtoms::intToGOAType(GOA_TYPE *goaType, int intGoaType) {
  if (intGoaType >= NUMBER_OF_GOA_TYPES) {
    *goaType = ATOM;
    return false;
  }
  switch (intGoaType) {
  case 0:
    *goaType = ATOM;
    return true;
  case 1:
    *goaType = RESIDUE;
    return true;
  case 2:
    *goaType = SECONDARY_STRUCTURE;
    return true;
  case 3:
    *goaType = CHAIN;
    return true;
  case 4:
    *goaType = PROTEIN;
    return true;
  case 5:
    *goaType = COLLECTION;
    return true;
  default:
    *goaType = ATOM;
    return false;
  }
}

double GroupOfAtoms::getFunction(PDBParser::GroupOfAtoms::FUNCTIONS function) {
  if (function == ELECTRON_DENSITY) {
    return 1;
  }
  if (type != RESIDUE_TYPE) {
    return 0;
  }
  // we have a residue
  if (function == PER_ATOM_HYDROPHOBICITY ||
      function == PER_RESIDUE_HYDROPHOBICITY) {
    if (resID < 0 || resID >= TOTAL_RES_TYPE) {
      return 0;
    }
    return residueTable[resID].hydrophobicity;
  }
  return 0;
}

void GroupOfAtoms::getAttributesRecursive(
    double *points, double *radii, double *min, double *max, double *weights,
    double *maxRadius, PDBParser::GroupOfAtoms::FUNCTIONS FUNCTION,
    int GOA_level, PDBParser::GroupOfAtoms::RADIUS_TYPE radius_type,
    int *current_atom_position, int maxLength) {
  if ((*current_atom_position) >= maxLength) {
    printf("This is odd\n");
    return;
  }
  int m, i;
  m = (int)m_SubGroups.size();
  for (i = 0; i < m; i++) {
    GroupOfAtoms *currSubGroup = m_SubGroups.at(i);
    if (currSubGroup) {
      currSubGroup->getAttributesRecursive(
          points, radii, min, max, weights, maxRadius, FUNCTION, GOA_level,
          radius_type, current_atom_position, maxLength);
    }
  }
  if (type == GOA_level) {
    int i;
    for (i = 0; i < m_MultiresAtoms.size(); i++) {
      double x = m_MultiresAtoms[i]->m_Position[0];
      double y = m_MultiresAtoms[i]->m_Position[1];
      double z = m_MultiresAtoms[i]->m_Position[2];
      if (points) {
        points[(*current_atom_position) * 3 + 0] = x;
        points[(*current_atom_position) * 3 + 1] = y;
        points[(*current_atom_position) * 3 + 2] = z;
      }
      double r = m_MultiresAtoms[i]->getRadius();
      if (radii) {
        radii[(*current_atom_position)] = r;
      }
      if (maxRadius) {
        if ((*maxRadius) < r) {
          (*maxRadius) = r;
        }
      }
      if (min && max) {
        if (min[0] < x) {
          min[0] = x;
        }
        if (min[1] < y) {
          min[1] = y;
        }
        if (min[2] < z) {
          min[2] = z;
        }
        if (max[0] > x) {
          max[0] = x;
        }
        if (max[1] > y) {
          max[1] = y;
        }
        if (max[2] > z) {
          max[2] = z;
        }
      }
      if (weights) {
        weights[(*current_atom_position)] = getFunction(FUNCTION);
      }
      (*current_atom_position)++;
    }
  }
  if ((GOA_level == ATOM_TYPE) && (type == RESIDUE_TYPE)) {
    int n = (int)m_Atoms.size();
    for (i = 0; i < n; i++) {
      Atom *cur_atom = m_Atoms.at(i);
      double x = cur_atom->m_Position[0];
      double y = cur_atom->m_Position[1];
      double z = cur_atom->m_Position[2];
      if (points) {
        points[(*current_atom_position) * 3 + 0] = x;
        points[(*current_atom_position) * 3 + 1] = y;
        points[(*current_atom_position) * 3 + 2] = z;
      }
      double r = cur_atom->getRadius();
      if (radii) {
        radii[(*current_atom_position)] = r;
      }
      if (maxRadius) {
        if ((*maxRadius) < r) {
          (*maxRadius) = r;
        }
      }
      if (min && max) {
        if (min[0] > x) {
          min[0] = x;
        }
        if (min[1] > y) {
          min[1] = y;
        }
        if (min[2] > z) {
          min[2] = z;
        }
        if (max[0] < x) {
          max[0] = x;
        }
        if (max[1] < y) {
          max[1] = y;
        }
        if (max[2] < z) {
          max[2] = z;
        }
      }
      if (weights) {
        weights[(*current_atom_position)] = cur_atom->getFunction(FUNCTION);
      }
      (*current_atom_position)++;
    }
  }
}

// Any of the variables can be NULL, stating that we do not need to compute them
bool GroupOfAtoms::getAttributes(
    double *points, double *radii, double *min, double *max, double *weights,
    double *maxRadius, PDBParser::GroupOfAtoms::FUNCTIONS function,
    int GOA_level, PDBParser::GroupOfAtoms::RADIUS_TYPE radius_type,
    int maxLength) {
  if (maxLength < 1) {
    return false;
  }
  int current_atom_position = 0;
  getAttributesRecursive(points, radii, min, max, weights, maxRadius, function,
                         GOA_level, radius_type, &current_atom_position,
                         maxLength);
  return true;
}

//    we assume that the transformation does not contain any scale or shear
bool GroupOfAtoms::transform(CCVOpenGLMath::Matrix transformationMatrix) {
  int i, n;
  n = m_SubGroups.size();
  for (i = 0; i < n; i++) {
    GroupOfAtoms *subGroup = m_SubGroups.at(i);
    if (subGroup) {
      subGroup->transform(transformationMatrix);
    }
  }
  n = m_Atoms.size();
  for (i = 0; i < n; i++) {
    Atom *atom = m_Atoms.at(i);
    if (atom) {
      atom->transform(transformationMatrix);
    }
  }
  n = m_MultiresAtoms.size();
  for (i = 0; i < n; i++) {
    Atom *atom = m_MultiresAtoms[i];
    if (atom) {
      atom->transform(transformationMatrix);
    }
  }
  // update min max
  {
    double minX, minY, minZ, maxX, maxY, maxZ;
    getMinMax(&minX, &minY, &minZ, &maxX, &maxY, &maxZ);
  }
  return true;
}

bool GroupOfAtoms::getNumberOfGOAsRecursive(
    int *count, PDBParser::GroupOfAtoms::GOA_TYPE goaType) {
  if (!count) {
    return false;
  }
  GOA_TYPE currentType;
  if (!intToGOAType(&currentType, type - 1)) {
    return false; // type-1 because there is some useless type 0, should be
                  // removed. SKVINAY
  }
  if (currentType == goaType) {
    (*count)++;
    return true;
  }
  if (goaType == PDBParser::GroupOfAtoms::ATOM) {
    *count += (int)m_Atoms.size();
  }
  int i, m;
  m = (int)m_SubGroups.size();
  for (i = 0; i < m; i++) {
    GroupOfAtoms *currSubGroup = m_SubGroups.at(i);
    currSubGroup->getNumberOfGOAsRecursive(count, goaType);
  }
  return true;
}

//  This does not use the transformations! SKVINAY
int GroupOfAtoms::getNumberOfGOAs(PDBParser::GroupOfAtoms::GOA_TYPE goaType) {
  int count = 0;
  if (!getNumberOfGOAsRecursive(&count, goaType)) {
    return 0;
  }
  return count;
}

void GroupOfAtoms::printColors(FILE *fp) {
  if (!fp) {
    return;
  }
  {
    int n = m_Atoms.size();
    int i;
    for (i = 0; i < n; i++) {
      float r, g, b, a;
      m_Atoms[i]->getColor(&r, &g, &b, &a);
      fprintf(fp, "[%c -> %2.1f %2.1f %2.1f %2.1f]\n", m_Atoms[i]->getChain(),
              r, g, b, a);
    }
  }
  int i, m;
  m = (int)m_SubGroups.size();
  for (i = 0; i < m; i++) {
    GroupOfAtoms *currSubGroup = m_SubGroups.at(i);
    currSubGroup->printColors(fp);
  }
}

bool GroupOfAtoms::isValidGOA(GroupOfAtoms *molecule, double overlapAllowed) {
  if (!molecule) {
    return false;
  }
  return true;
}

int GroupOfAtoms::getMassCenter(double *x, double *y, double *z) {
  // x,y,z will contain the sum of all coordinates
  // returns the number of atoms
  // divide x,y,z by the return value to get the average center
  int count = 0;
  double tempx, tempy, tempz;
  if (m_SubGroups.size() > 0) {
    for (int i = 0; i < m_SubGroups.size(); i++) {
      tempx = 0;
      tempy = 0;
      tempz = 0;
      count += m_SubGroups[i]->getMassCenter(&tempx, &tempy, &tempz);
      *x += tempx;
      *y += tempy;
      *z += tempz;
    }
    return count;
  }
  if (m_Atoms.size() > 0) {
    for (int i = 0; i < m_Atoms.size(); i++) {
      *x += m_Atoms[i]->m_Position[0];
      *y += m_Atoms[i]->m_Position[1];
      *z += m_Atoms[i]->m_Position[2];
    }
    count = m_Atoms.size();
    return count;
  }
}
