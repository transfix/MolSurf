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
#include <PDBParser/Residues.h>

using namespace PDBParser;
using namespace CCVOpenGLMath;

Residues::Residues() {}

Residues::~Residues() {}

int Residues::mapResidue(char *resName) {
  int i;
  for (i = 1; i < TOTAL_RES_TYPE;
       i++) // 0th residue is a default, dont compare with that
  {
    if (strcmpCaseInsensitive(resName, residueTable[i].residueName)) {
      return i;
    }
  }
  // Print a warning with the name of the unknown residue.
  // This block of code is a pretty round-about way to be able to print a
  // warning with a constant string followed by a variable string.
  // fprintf( stderr, "mapResidue: unknown residue - %s\n", resName );

  // arand commented since verbosity control has been removed.
  // printf("mapResidue: unknown residue - %s\n", resName);
  return RES_TYPE_DEFAULT;
}

char Residues::lookupShortResID(char *resName) {
  int i;
  for (i = 0; i < TOTAL_RES_TYPE; i++) {
    if (strcmpCaseInsensitive(resName, residueTable[i].residueName)) {
      return residueTable[i].shortResidueName;
    }
  }
  // Print a warning with the name of the unknown residue.
  // This block of code is a pretty round-about way to be able to print a
  // warning with a constant string followed by a variable string.
  // fprintf( stderr, "mapResidue: unknown residue - %s\n", resName );

  // arand commented since verbosity control has been removed.
  // printf("mapResidue: unknown residue - %s\n", resName);
  return RES_TYPE_DEFAULT_SHORT_ID;
}

const char *Residues::lookupLongResID(char resShortID) {
  int i;
  for (i = 0; i < TOTAL_RES_TYPE; i++) {
    if (resShortID == residueTable[i].shortResidueName) {
      return residueTable[i].residueName;
    }
  }
  // Print a warning with the name of the unknown residue.
  // This block of code is a pretty round-about way to be able to print a
  // warning with a constant string followed by a variable string.
  // fprintf( stderr, "mapResidue: unknown residue - %s\n", resName );

  // arand commented since verbosity control has been removed.
  // printf("mapResidue: unknown residue - %c\n", resShortID);
  return RES_TYPE_DEFAULT_ID;
}

Atom *Residues::getAtom(GroupOfAtoms *residue, const char *atomName) {
  int j, m;
  m = (int)residue->m_Atoms.size();
  for (j = 0; j < m; j++) {
    Atom *currAtom = residue->m_Atoms.at(j);
    if (strcmpCaseInsensitive(currAtom->getName(), atomName)) {
      return currAtom;
    }
  }
  return 0;
}

bool Residues::getAngle(int angleNumber, GroupOfAtoms *residue, double *angle) {
  *angle = 0;
  if (residue->type != RESIDUE_TYPE) {
    return false;
  }
  const char *atomName1;
  const char *atomName2;
  const char *atomName3;
  const char *atomName4;
  if (angleNumber == AA_CHI1) {
    atomName1 = AminoAcidChi1Table[residue->resID][1];
    atomName2 = AminoAcidChi1Table[residue->resID][2];
    atomName3 = AminoAcidChi1Table[residue->resID][3];
    atomName4 = AminoAcidChi1Table[residue->resID][4];
  } else if (angleNumber == AA_CHI2) {
    atomName1 = AminoAcidChi2Table[residue->resID][1];
    atomName2 = AminoAcidChi2Table[residue->resID][2];
    atomName3 = AminoAcidChi2Table[residue->resID][3];
    atomName4 = AminoAcidChi2Table[residue->resID][4];
  } else if (angleNumber == AA_CHI3) {
    atomName1 = AminoAcidChi3Table[residue->resID][1];
    atomName2 = AminoAcidChi3Table[residue->resID][2];
    atomName3 = AminoAcidChi3Table[residue->resID][3];
    atomName4 = AminoAcidChi3Table[residue->resID][4];
  } else if (angleNumber == AA_CHI4) {
    atomName1 = AminoAcidChi4Table[residue->resID][1];
    atomName2 = AminoAcidChi4Table[residue->resID][2];
    atomName3 = AminoAcidChi4Table[residue->resID][3];
    atomName4 = AminoAcidChi4Table[residue->resID][4];
  } else if (angleNumber == AA_CHI5) {
    atomName1 = AminoAcidChi5Table[residue->resID][1];
    atomName2 = AminoAcidChi5Table[residue->resID][2];
    atomName3 = AminoAcidChi5Table[residue->resID][3];
    atomName4 = AminoAcidChi5Table[residue->resID][4];
  } else {
    return false;
  }
  if (!strcmp(atomName1, "    ")) {
    return false;
  }
  if (!strcmp(atomName2, "    ")) {
    return false;
  }
  if (!strcmp(atomName3, "    ")) {
    return false;
  }
  if (!strcmp(atomName4, "    ")) {
    return false;
  }
  Atom *nx = getAtom(residue, atomName1);
  Atom *nr = getAtom(residue, atomName2);
  Atom *fr = getAtom(residue, atomName3);
  Atom *fx = getAtom(residue, atomName4);
  if ((nx == 0) || (nr == 0) || (fr == 0) || (fx == 0)) {
    return false;
  }
  *angle = findTorsion(nx, nr, fr, fx);
  return true;
}

bool Residues::transformFrom(GroupOfAtoms *residue, int angleType, Matrix *mat,
                             const char *pivotAtomName) {
  // Transform points in the residue
  int j, m;
  // Isoleucine is the one misbehaving residue since its order is:
  // CA, CB, CG1, CG2, CD1.
  // Torsion bond: CA - CB - CG1 - CD1 should affect only CD1, but CG2 is in the
  // way. So just transfomr CD1 directly, and return.
  if (residue->resID == mapResidue((char *)"ILE") && angleType == AA_CHI2) {
    Atom *currAtom = findAtom(residue, " CD1");
    if (currAtom == 0) {
      return false;
    }
    Vector v(currAtom->m_Position[0], currAtom->m_Position[1],
             currAtom->m_Position[2], 1);
    Vector v2;
    v2.set((*mat) * v);
    currAtom->m_Position[0] = v2[0];
    currAtom->m_Position[1] = v2[1];
    currAtom->m_Position[2] = v2[2];
    return true;
  }
  // Iterate:
  // If the pivot angle is detected, set the flag to indicate we should apply
  // the transformation starting from the NEXT atom
  bool applyMat = false;
  m = (int)residue->m_Atoms.size();
  for (j = 0; j < m; j++) {
    Atom *currAtom = residue->m_Atoms.at(j);
    if (applyMat) {
      Vector v(currAtom->m_Position[0], currAtom->m_Position[1],
               currAtom->m_Position[2], 1);
      Vector v2;
      v2.set((*mat) * v);
      currAtom->m_Position[0] = v2[0];
      currAtom->m_Position[1] = v2[1];
      currAtom->m_Position[2] = v2[2];
    }
    if (applyMat == false &&
        (strcmpCaseInsensitive(pivotAtomName, currAtom->getName()))) {
      applyMat = true;
    }
  }
  // we should have transformed atleast 1 atom
  if (applyMat == false) {
    return false;
  }
  return true;
}

bool Residues::getRotationVector(GroupOfAtoms *residue,
                                 const char *vecTailAtomName,
                                 const char *pivotAtomName, double *vecx,
                                 double *vecy, double *vecz) {
  Atom *pAtom = findAtom(residue, pivotAtomName);
  Atom *tAtom = findAtom(residue, vecTailAtomName);
  if (pAtom == 0 || tAtom == 0) {
    return false;
  }
  *vecx = pAtom->m_Position[0] - tAtom->m_Position[0];
  *vecy = pAtom->m_Position[1] - tAtom->m_Position[1];
  *vecz = pAtom->m_Position[2] - tAtom->m_Position[2];
  return true;
}

bool Residues::getRotationOrigin(GroupOfAtoms *residue,
                                 const char *pivotAtomName, double *origx,
                                 double *origy, double *origz) {
  Atom *pAtom = findAtom(residue, pivotAtomName);
  if (pAtom == 0) {
    return false;
  }
  *origx = pAtom->m_Position[0];
  *origy = pAtom->m_Position[1];
  *origz = pAtom->m_Position[2];
  return true;
}

bool Residues::transform(GroupOfAtoms *residue, int angleType, double newAngle,
                         const char *vecTailAtomName,
                         const char *pivotAtomName) {
  // Find the deltas we need
  double old_chi = 0;
  if (!getAngle(angleType, residue, &old_chi)) {
    return false;
  }
  double d_chi = newAngle - old_chi;
  double vecx = 1, vecy = 0, vecz = 0;
  if (!getRotationVector(residue, vecTailAtomName, pivotAtomName, &vecx, &vecy,
                         &vecz)) {
    return false;
  }
  double origx = 0, origy = 0, origz = 0;
  if (!getRotationOrigin(residue, pivotAtomName, &origx, &origy, &origz)) {
    return false;
  }
  Matrix ntrans = Matrix::translation(-origx, -origy, -origz);
  Matrix rot = Quaternion::rotation(d_chi, vecx, vecy, vecz).buildMatrix();
  Matrix ptrans = Matrix::translation(origx, origy, origz);
  Matrix mat = ptrans * rot * ntrans;
  return transformFrom(residue, angleType, &mat, pivotAtomName);
}
