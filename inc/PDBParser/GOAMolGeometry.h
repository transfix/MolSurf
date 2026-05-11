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
#ifndef __GOAMOLGEOMETRY_H__
#define __GOAMOLGEOMETRY_H__

#include <Utility/utility.h>
#include <PDBParser/Atom.h>
#include <PDBParser/GroupOfAtoms.h>
#include <PDBParser/ResidueHeader.h>
#include <PDBParser/ResidueInformation.h>
#include <PDBParser/ResidueData.h>
#include <PDBParser/GOAVector.h>

namespace PDBParser {
class GroupOfAtoms;
class Atom;
void calculateNucleotideTorsionAngles(Atom *prevO3p, GroupOfAtoms *nuc,
                                      Atom *nextP, Atom *nextO5p, Atom *prevC4p,
                                      Atom *nextC4p);
void calculateNucleicAcidChainTorsionAngles(GroupOfAtoms *currChain);
double findTorsion(Atom *nx, Atom *nr, Atom *fr, Atom *fx);
int strcmp_ci2(const char *str1, const char *str2);
bool isPurine(GroupOfAtoms *nuc);
bool isPyrimidine(GroupOfAtoms *nuc);
Atom *findAtom(GroupOfAtoms *res, const char *name);
GroupOfAtoms *findResidue(GroupOfAtoms *mol, int reqResID);
void calculatePhiAndPsi(GroupOfAtoms *previousResidue,
                        GroupOfAtoms *currentResidue);
void calculateProteinChainTorsionAngles(GroupOfAtoms *currChain);
}; // namespace PDBParser

#endif
