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
#ifndef CCV_TEXMOL_MUTATED_ATOM_INFORMATION_H
#define CCV_TEXMOL_MUTATED_ATOM_INFORMATION_H

#include <Utility/utility.h>
#include <PDBParser/ElementInformation.h>
#include <PDBParser/parserPDBtoGOA.h>
#include <PDBParser/ResidueInformation.h>

namespace PDBParser {
/*  This class contains information for non - standard atoms.  */
/*  Standard atoms have an elementIndex in the element         */
/*  information table. The non - standard atoms will have an   */
/*  index of -1 and requires this structure which helps in     */
/*  providing more information.                                */
class MutatedAtomInformation {
public:
  MutatedAtomInformation(float radius, const char *name,
                         unsigned char residueID, float charge);
  MutatedAtomInformation(const MutatedAtomInformation &m);
  virtual ~MutatedAtomInformation();
  unsigned char getResidueId();
  void setResidueId(unsigned char residueId);
  bool getColor(float *red, float *green, float *blue);
  float m_Radius;
  char m_Name[5];
  float m_Charge;
  float color[3];

protected:
  void resetColor();
  unsigned char m_ResidueId;
};

}; // namespace PDBParser

#endif
