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
// GOAMolGeometry.cpp: Calculates geometric properties of a molecule based on
// its structure.  For instance, this file contains methods to find the
// torsion angles along protein and RNA chain backbones.

#include <PDBParser/GOAMolGeometry.h>

using namespace PDBParser;

// findTorsion
// -----------
// Given a reference bond (specified by two atoms) and the two other atoms
// adjacent to the bond that are used to measure the torsion about the reference
// bond, finds the torsion angle.  The angle given is the clockwise rotation
// needed to go from the near axis atom to the far axis atom if you are looking
// down the reference bond with the near reference atom closest to you.  If any
// arguments are NULL, returns 0.
//
// Arguments:
// nx - near axis atom, this is the atom bonded to the near reference atom and
// used as one axis of the torsion nr - near reference atom, this is the atom on
// the reference bond closer to the viewpoint if one is looking
//      down the reference bond
// fr - far reference atom, this is the atom on the reference bond further
// fromthe viewpoint if one is looking
//      down the reference bond
// fx - far axis atom, this is the atom bonded to the far reference atom and
// used as one axis of the torsion
//
// Pictorially, the relative location of the atoms looks like this.
//
// nx    fr
//   \  /  \
//    nr    fx
//
// The center bond is the reference bond, and the nr-nx and fr-fx bonds are used
// to find the torsion around the center bond.
//
// Mathmatically, nx and fx are projected to a plane which is perpendicular to
// the reference bond.  The torsion angle found will be the counterclockwise
// rotation of the projection of nx to the projection of fx (i.e.
// counterclockwise rotation is positive, clockwise rotation is negative).
double PDBParser::findTorsion(Atom *nx, Atom *nr, Atom *fr, Atom *fx) {
  if (!nx || !nr || !fr || !fx) {
    return 0;
  }
  // The three bonds connecting nx, nr, fr, and fx.  All of these vectors are
  // oriented so that they point towards the atom closer to nx.
  GOAVector *vprev = 0; // The nr-nx bond.
  GOAVector *vnext = 0; // The fx-fr bond.
  GOAVector *vbond = 0; // The reference bond (fr-nr).
  GOAVector *plane1n =
      0; // The normal to the plane containing the first three atoms
  GOAVector *plane2n =
      0; // The normal to the plane containing the last three atoms
  GOAVector *temp_v1 = 0; // used for orientation, always points in the positive
                          // side for orientation
  double temp_scaler = 0; // The final angle is stored here
  if (nx == 0 || nr == 0 || fr == 0 || fx == 0) {
    return 0;
  }
  vprev = subtractAtoms(nx, nr);
  vbond = subtractAtoms(nr, fr);
  vnext = subtractAtoms(fr, fx);
  // Algorithm: to find the angle between vprev and vbond after both are
  // projected into a plane, find the angle between the normals of plane1 and
  // plane2.
  if (vbond == 0 || vprev == 0) {
    return 0;
  }
  plane1n = crossProduct(
      vbond, vprev); // find the plan containing the first three atoms
  if (vbond == 0 || vnext == 0) {
    return 0;
  }
  plane2n = crossProduct(
      vnext, vbond); // find the plan containing the last three atoms
  // Make unit length
  plane1n->norm();
  plane2n->norm();
  // Get the angle without orientation (aka positive or negative sign)
  temp_scaler = dotProductAngle(plane1n, plane2n);
  // Get the orientation (aka positive or negative sign) of the angle
  if (vbond == 0 || plane1n == 0) {
    return 0;
  }
  temp_v1 = crossProduct(plane1n, vbond);
  if (dotProduct(temp_v1, plane2n) < 0) {
    temp_scaler *= -1;
  }
  return temp_scaler;
}

/*  Calculate psi for the previous residue and phi for the current residue */
/*  Either of both can be NULL. */
void PDBParser::calculatePhiAndPsi(GroupOfAtoms *previousResidue,
                                   GroupOfAtoms *currentResidue) {
  // get currentResidues N, CA, C atoms and previousResidues C atom
  // get phi, psi from these 4 atoms.
  // previous residue could be null.
  // psi angle (i-1) = [ N(i-1) - CA(i-1) - C'(i-1) - N(i) ]
  if (!previousResidue) {
    // nothing to be done
  } else if (!currentResidue) {
    if (!previousResidue->m_ResidueData) {
      previousResidue->m_ResidueData = new ResidueData();
    }
    previousResidue->m_ResidueData->m_Psi = 0;
  } else {
    if (!previousResidue->m_ResidueData) {
      previousResidue->m_ResidueData = new ResidueData();
    }
    Atom *Nprev = findAtom(previousResidue, " N  ");
    Atom *CAprev = findAtom(previousResidue, " CA ");
    Atom *Cprev = findAtom(previousResidue, " C  ");
    Atom *N = findAtom(currentResidue, " N  ");
    previousResidue->m_ResidueData->m_Psi =
        findTorsion(Nprev, CAprev, Cprev, N);
  }
  // phi angle (i) = [ C'(i-1) - N(i) - CA(i) - C'(i) ]
  if (!currentResidue) {
    // nothing to be done
  } else if (!previousResidue) {
    if (!currentResidue->m_ResidueData) {
      currentResidue->m_ResidueData = new ResidueData();
    }
    currentResidue->m_ResidueData->m_Phi = 0;
  } else {
    if (!currentResidue->m_ResidueData) {
      currentResidue->m_ResidueData = new ResidueData();
    }
    Atom *Cprev = findAtom(previousResidue, " C  ");
    Atom *N = findAtom(currentResidue, " N  ");
    Atom *CA = findAtom(currentResidue, " CA ");
    Atom *C = findAtom(currentResidue, " C  ");
    currentResidue->m_ResidueData->m_Phi = findTorsion(Cprev, N, CA, C);
  }
}

void PDBParser::calculateProteinChainTorsionAngles(GroupOfAtoms *currChain) {
  if (!currChain) {
    return;
  }
  if (currChain->m_SubGroups.size() < 1) {
    fprintf(stdout, "AminoAcidTorsion: warning: no secondary structure, no "
                    "angles calculated\n");
    return;
  }
  GroupOfAtoms *previousResidue = 0;
  int i;
  for (i = 0; i < currChain->m_SubGroups.size(); i++) {
    GroupOfAtoms *currentSS = currChain->m_SubGroups.at(i);
    if (!currentSS) {
      continue;
    }
    // go through each secondary structure
    int j;
    for (j = 0; j < currentSS->m_SubGroups.size(); j++) {
      GroupOfAtoms *currentResidue = currentSS->m_SubGroups.at(j);
      if (!currentResidue) {
        continue;
      }
      if (previousResidue) {
        // we have a previous and a current residue, calculate some angles !
        // We calculate the psi of the previous residue and the
        // phi of the current residue
        calculatePhiAndPsi(previousResidue, currentResidue);
      } else {
        // this is the first residue I am looking at. Calculate phi for it.
        calculatePhiAndPsi(0, currentResidue); // this simply sets phi to 0.
      }
      previousResidue =
          currentResidue; // is this the right way ? This is one way.
    }
  }
  // need to call with (previousResidue, 0) to set previousResidue's psi to 0.
  calculatePhiAndPsi(previousResidue, 0); // this simply sets psi to 0.
}

// Case insensitive string compare.  returns 0 if the two input strings are
// case-insensitive identical (in length and content) and returns non-zero
// otherwise.  Returns non-zero if either string is NULL.
int PDBParser::strcmp_ci2(const char *str1, const char *str2) {
  unsigned int i, n;
  if (str1 == 0 || str2 == 0) {
    return 1;
  }
  n = strlen(str1);
  if (strlen(str2) != n) {
    return 1; // unequal length
  }
  for (i = 0; i < n; i++) {
    if (islower(str2[i])) {
      if (tolower(str1[i]) != str2[i]) {
        return 1; // mismatch
      }
    } else if (isupper(str2[i])) {
      if (toupper(str1[i]) != str2[i]) {
        return 1; // mismatch
      }
    } else {
      if (str1[i] != str2[i]) {
        return 1; // mismatch
      }
    }
  }
  // All chars have matched
  return 0;
}

// Globals that contain the residue ID's for the nucleotide bases.
int resID_A = RES_TYPE_DEFAULT;
int resID_G = RES_TYPE_DEFAULT;
int resID_C = RES_TYPE_DEFAULT;
int resID_T = RES_TYPE_DEFAULT;
int resID_U = RES_TYPE_DEFAULT;

// Uses the residue ID's from "ResidueInformation.h"
bool PDBParser::isPurine(GroupOfAtoms *nuc) {
  if (!nuc) {
    return false;
  }
  int id = nuc->resID;
  if (id == resID_A || id == resID_G) {
    return true;
  }
  return false;
}

// Uses the residue ID's from "ResidueInformation.h"
bool PDBParser::isPyrimidine(GroupOfAtoms *nuc) {
  if (!nuc) {
    return false;
  }
  int id = nuc->resID;
  if (id == resID_U || id == resID_C || id == resID_T) {
    return true;
  }
  return false;
}

// Given a residue and an atom name, tries to find the named atom.  Returns
// NULL if the atom is not present.
Atom *PDBParser::findAtom(GroupOfAtoms *res, const char *name) {
  if (!res) {
    return 0;
  }
  if (!name) {
    return 0;
  }
  int i, n;
  Atom *curr = 0;
  n = (int)res->m_Atoms.size();
  for (i = 0; i < n; i++) {
    curr = res->m_Atoms.at(i);
    if (strcmp_ci2(curr->getName(), name) == 0) {
      // match found
      return curr;
    }
  }
  // no match found
  if (name != 0 && res != 0) {
    // arand, this message gives a lot of output...
    //        so I commented it
    // printf("findAtom: warning, requested atom %s not found in residue w/ 1st
    // atom: %d\n", name, res->resPos);
  }
  return 0;
}

// Given a molecule and a residue name, tries to find the named residue.
// It returns the first instance found.
// Returns 0 if the residue is not present.
GroupOfAtoms *PDBParser::findResidue(GroupOfAtoms *mol, int reqResID) {
  if (!mol) {
    return 0;
  }
  // Base case
  if (mol->type == RESIDUE_TYPE) {
    if (mol->resID == reqResID) {
      return mol;
    } else {
      return 0;
    }
  }
  // Recurse
  GroupOfAtoms *ret = 0;
  int i;
  int n = (int)mol->m_SubGroups.size();
  for (i = 0; i < n; i++) {
    GroupOfAtoms *currSubGroup = mol->m_SubGroups.at(i);
    if (currSubGroup) {
      ret = findResidue(currSubGroup, reqResID);
      if (ret != 0) {
        break;
      }
    }
  }
  return ret;
}

// This function processes a given nucleotide's torsion angles.
// The calculated angles are stored in data members inside of "nuc".
//
// Assumption: Each phosphate group is assumed to be grouped with the residue
// that has its O5' bonded with the phosphate group.  The algorithm should be
// changed if this is not the case.
//
// Thus, to calculate all torsion angles for a nucleotide, the previous
// nucleotide's O3' needs to be passed in.  If there is no previous O3', then
// NULL should be passed in.
//
// Likewise, the next nucleotide's P and O5' need to be passed in.  If there is
// no next P or O5', then NULL should be passed in to the corresponding
// argument.
void PDBParser::calculateNucleotideTorsionAngles(Atom *prevO3p,
                                                 GroupOfAtoms *nuc, Atom *nextP,
                                                 Atom *nextO5p, Atom *prevC4p,
                                                 Atom *nextC4p) {
  // check to ensure we are dealing with a residue
  if (nuc->type != PDBParser::RESIDUE_TYPE) {
    return;
  }
  // create ResidueData if its not already there. We should init all pointers to
  // 0
  if (!(nuc->m_ResidueData)) {
    nuc->m_ResidueData = new ResidueData();
  }
  // A 'p' character stands for prime in the atom names
  Atom *P = findAtom(nuc, " P  ");
  Atom *O5p = findAtom(nuc, " O5*");
  Atom *C5p = findAtom(nuc, " C5*");
  Atom *C4p = findAtom(nuc, " C4*");
  Atom *C3p = findAtom(nuc, " C3*");
  Atom *O3p = findAtom(nuc, " O3*");
  Atom *O4p = findAtom(nuc, " O4*");
  Atom *C1p = findAtom(nuc, " C1*");
  Atom *C2p = findAtom(nuc, " C2*");
  Atom *baseN = 0; // The N atom in the base needed to calculate chi
  Atom *baseC = 0; // The C atom in the base needed to calculate chi
  if (isPurine(nuc)) {
    baseN = findAtom(nuc, " N9 ");
    baseC = findAtom(nuc, " C4 ");
    if (!baseN || !baseC) {
      return;
    }
  } else if (isPyrimidine(nuc)) {
    baseN = findAtom(nuc, " N1 ");
    baseC = findAtom(nuc, " C2 ");
    if (!baseN || !baseC) {
      return;
    }
  } else {
    fprintf(stderr, "calculate nuc torsion: warning, nucleotide is neither a "
                    "purine nor a pyrimidine!\n");
  }
  // If any of the atom pointers are NULL, the particular torsion angle will be
  // set to 0
  nuc->m_ResidueData->m_alpha = findTorsion(prevO3p, P, O5p, C5p);
  nuc->m_ResidueData->m_beta = findTorsion(P, O5p, C5p, C4p);
  nuc->m_ResidueData->m_gamma = findTorsion(O5p, C5p, C4p, C3p);
  nuc->m_ResidueData->m_delta = findTorsion(C5p, C4p, C3p, O3p);
  nuc->m_ResidueData->m_epsilon = findTorsion(C4p, C3p, O3p, nextP);
  nuc->m_ResidueData->m_zeta = findTorsion(C3p, O3p, nextP, nextO5p);
  nuc->m_ResidueData->m_nu0 = findTorsion(C4p, O4p, C1p, C2p);
  nuc->m_ResidueData->m_nu1 = findTorsion(O4p, C1p, C2p, C3p);
  nuc->m_ResidueData->m_nu2 = findTorsion(C1p, C2p, C3p, C4p);
  nuc->m_ResidueData->m_nu3 = findTorsion(C2p, C3p, C4p, O4p);
  nuc->m_ResidueData->m_nu4 = findTorsion(C3p, C4p, O4p, C1p);
  nuc->m_ResidueData->m_chi = findTorsion(O4p, C1p, baseN, baseC);
  nuc->m_ResidueData->m_Eta = findTorsion(prevC4p, P, C4p, nextP);
  nuc->m_ResidueData->m_Theta = findTorsion(P, C4p, nextP, nextC4p);
  if (prevC4p == 0) {
    fprintf(stderr, "prevC4p is null\n");
  }
  if (P == 0) {
    fprintf(stderr, "P is null\n");
  }
  if (C4p == 0) {
    fprintf(stderr, "C4p is null\n");
  }
  if (nextP == 0) {
    fprintf(stderr, "nextP is null\n");
  }
  // Puckering
  GOAVector *v1 = 0;
  GOAVector *v2 = 0;
  GOAVector *plane1n =
      0; // The normal to the plane containing C4p, O4p, and C1p
  GOAVector *plane2n =
      0; // The normal to the plane containing O4p, C1p, and C2p
  GOAVector *plane3n =
      0;               // The normal to the plane containing C3p, C4p, and O4p
  double temp2, temp3; // Holds dot product temporaries
  // Find plane1n
  v1 = subtractAtoms(C4p, O4p);
  v2 = subtractAtoms(C1p, O4p);
  if (v1 == 0 || v2 == 0) {
    return;
  }
  plane1n = crossProduct(v1, v2);
  if (!plane1n) {
    return;
  }
  plane1n->norm(); // Normalize
  // pcd debug: better algorithm that needs to be implemented -- get the vectors
  // from O4p to C2p and C3p, respectively then dot product these with the
  // normal to plane1 Find plane2
  v1 = subtractAtoms(O4p, C1p);
  v2 = subtractAtoms(C2p, C1p);
  if (v1 == 0 || v2 == 0) {
    return;
  }
  plane2n = crossProduct(v1, v2);
  if (!plane2n) {
    return;
  }
  plane2n->norm(); // Normalize
  // Find plane3
  v1 = subtractAtoms(C3p, C4p);
  v2 = subtractAtoms(O4p, C4p);
  if (v1 == 0 || v2 == 0) {
    return;
  }
  plane3n = crossProduct(v1, v2);
  if (!plane3n) {
    return;
  }
  plane3n->norm(); // Normalize
  // Find whether plane2n or plane3n is closer to plane1n
  temp2 = dotProduct(plane1n, plane2n);
  temp3 = dotProduct(plane1n, plane3n);
  //     printf("2: %f,    3: %f\n", temp2, temp3);
}

// Main Nucleic Acid Torsion Angle Routine
//
// Goes along the backbone of an nucleic acid chain and calculates the various
// bond angles for each nucleotide.
//
// Assumption: as we step through the atoms in order from the pdb, we are going
// from the 5' end to the 3' end.  Also, each phosphate group is assumed to be
// grouped with the residue that has its O3' bonded with the phosphate group.
// The algorithm should be changed if these assumptions are not true.
void PDBParser::calculateNucleicAcidChainTorsionAngles(
    GroupOfAtoms *currChain) {
  if (!currChain) {
    return;
  }
  int i, n;
  // Initialize the values of the residue ID's for the nucleotide types
  for (i = 1; i < TOTAL_RES_TYPE;
       i++) // 0th residue is a default, dont compare with that
  {
    if (strcmp_ci2("  A", residueTable[i].residueName) == 0) {
      resID_A = i;
    }
    if (strcmp_ci2("  G", residueTable[i].residueName) == 0) {
      resID_G = i;
    }
    if (strcmp_ci2("  C", residueTable[i].residueName) == 0) {
      resID_C = i;
    }
    if (strcmp_ci2("  T", residueTable[i].residueName) == 0) {
      resID_T = i;
    }
    if (strcmp_ci2("  U", residueTable[i].residueName) == 0) {
      resID_U = i;
    }
  }
  // Three pointers that walk down the chain
  GroupOfAtoms *prevNuc = 0;
  GroupOfAtoms *currNuc = 0;
  GroupOfAtoms *nextNuc = 0;
  // Pointers to atoms we need from other residues.  The 'p' character stands
  // for prime.
  Atom *prevO3p = 0;
  Atom *nextP = 0;
  Atom *nextO5p = 0;
  Atom *prevC4p = 0;
  Atom *nextC4p = 0;
  // Names to strings, so they do not have to be reallocated in the loop body
  char *strO3p = (char *)" O3*";
  char *strP = (char *)" P  ";
  char *strO5p = (char *)" O5*";
  char *strC4p = (char *)" C4*";
  // Since a chain has a single SS_NULL secondary structure GOA between it and
  // the residues, we need a pointer to it
  GroupOfAtoms *ssNull;
  if (currChain->m_SubGroups.size() == 1) {
    ssNull = currChain->m_SubGroups.at(0);
    n = ssNull->m_SubGroups.size();
    if (n <= 0) {
      fprintf(stderr, "NucleicAcidTorsion: warning: no nucleotides, no angles "
                      "calculated\n");
      return;
    }
  } else {
    fprintf(
        stderr,
        "NucleicAcidTorsion: warning: empty or improperly formatted chain\n");
    return;
  }
  prevNuc = 0; // no previous nucleotide for the first one in the chain
  currNuc = ssNull->m_SubGroups.at(0);
  if (!currNuc) {
    return;
  }
  for (i = 1; i < n; i++) {
    nextNuc = ssNull->m_SubGroups.at(i);
    // Find atoms we need from neighboring nucleotides
    prevO3p = findAtom(prevNuc, strO3p);
    nextP = findAtom(nextNuc, strP);
    nextO5p = findAtom(nextNuc, strO5p);
    nextC4p = findAtom(nextNuc, strC4p);
    prevC4p = findAtom(prevNuc, strC4p);
    if (nextC4p == 0) {
      fprintf(stderr, "Next C4p was null\n");
    }
    // Find the angles
    calculateNucleotideTorsionAngles(prevO3p, currNuc, nextP, nextO5p, prevC4p,
                                     nextC4p);
    // Update the pointers
    prevNuc = currNuc;
    currNuc = nextNuc;
  }
  // The loop calculates the torsion angles for all nucleotides in the chain
  // except for the last one, so process the last one now.
  prevNuc = currNuc;
  currNuc = nextNuc;
  nextNuc = 0;
  prevO3p = findAtom(prevNuc, strO3p);
  prevC4p = findAtom(prevNuc, strC4p);
  if (!prevO3p || !prevNuc || !currNuc) {
    return;
  }
  nextP = 0;   // nextNuc is NULL, so there is no nextP
  nextO5p = 0; // nextNuc is NULL, so there is no nextO5p
  nextC4p = 0; // nextNuc is NULL, so there is no nextC4p
  calculateNucleotideTorsionAngles(prevO3p, currNuc, nextP, nextO5p, prevC4p,
                                   nextC4p);
  return;
}
