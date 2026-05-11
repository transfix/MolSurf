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
#if !defined(AFX_RESIDUES_H__05FCD023_EE56_4F87_BE83_656B87A25B73__INCLUDED_)
#define AFX_RESIDUES_H__05FCD023_EE56_4F87_BE83_656B87A25B73__INCLUDED_

#include <Utility/utility.h>
#include <UsefulMath/Matrix.h>
#include <PDBParser/GroupOfAtoms.h>
#include <PDBParser/GOAMolGeometry.h>
#include <PDBParser/Atom.h>
#include <PDBParser/ResidueInformation.h>
#include <PDBParser/parserPDBtoGOA.h>
#include <UsefulMath/Vector.h>
#include <UsefulMath/Matrix.h>
#include <UsefulMath/Quaternion.h>

namespace PDBParser {
class GroupOfAtoms;
class Atom;
class Residues {
public:
  Residues();
  virtual ~Residues();
  static int mapResidue(char *resName);
  static char lookupShortResID(char *resName);
  static const char *lookupLongResID(char resShortID);
  static bool getAngle(int angleNumber, GroupOfAtoms *residue, double *angle);
  static Atom *getAtom(GroupOfAtoms *residue, const char *atomName);
  static bool transform(GroupOfAtoms *residue, int angleType, double newAngle,
                        const char *vecTailAtomName, const char *pivotAtomName);
  static bool transformFrom(GroupOfAtoms *residue, int angleType,
                            CCVOpenGLMath::Matrix *mat,
                            const char *pivotAtomName);
  static bool getRotationVector(GroupOfAtoms *residue,
                                const char *vecTailAtomName,
                                const char *pivotAtomName, double *vecx,
                                double *vecy, double *vecz);
  static bool getRotationOrigin(GroupOfAtoms *residue,
                                const char *pivotAtomName, double *origx,
                                double *origy, double *origz);
};
}; // namespace PDBParser

#endif
