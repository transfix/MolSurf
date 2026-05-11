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
#ifndef __GOAFILEIO_H__
#define __GOAFILEIO_H__

#include <Utility/utility.h>
#include <UsefulMath/Matrix.h>
#include <PDBParser/GroupOfAtoms.h>
#include <PDBParser/ResidueHeader.h>
#include <PDBParser/ResidueData.h>
#include <PDBParser/CollectionData.h>
#include <PDBParser/ResidueInformation.h>
#include <PDBParser/ElementInformation.h>
#include <PDBParser/Atom.h>
#include <UsefulMath/Vector.h>

namespace PDBParser {
class GroupOfAtoms;
class Atom;
void printAtomPDB(FILE *stream, Atom *atom,
                  CCVOpenGLMath::Matrix *transformation);
void printGOA(GroupOfAtoms *currGOA, int indent);
int writeGOA2PTS(char *filename, GroupOfAtoms *currGOA);
bool writeGOA2PDB(FILE *stream, GroupOfAtoms *currGOA, int type,
                  CCVOpenGLMath::Matrix *transformation);
bool writeGOAtoFile(FILE *stream, GroupOfAtoms *currGOA);
bool writeTorsionAngles(FILE *stream, GroupOfAtoms *currGOA);
}; // namespace PDBParser

void printAtomPDB(FILE *stream, PDBParser::Atom *atom,
                  CCVOpenGLMath::Matrix *transformation);

#endif
