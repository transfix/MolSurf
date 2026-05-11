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
#if !defined(_RESIDUEHEADER_H_)
#define _RESIDUEHEADER_H_

#include <Utility/utility.h>

namespace PDBParser {
class ResidueHeader {
public:
  ResidueHeader();
  ResidueHeader(int a_resID, char a_chainID, int a_seqNum, char a_iCode);
  ResidueHeader(const ResidueHeader &rh);
  virtual ~ResidueHeader();
  // These fields are needed to describe a residue
  int resID;    // An int corresponding to the amino acid that the residue falls
                // under
  char chainID; // The parent chain of this residue
  int seqNum;   // The sequence number that is used, along with the insertion
                // code, in the PDB
  char iCode;   // The insertion code that is used, along with the sequence
                // number, in the PDB
  bool compare(ResidueHeader *rh);
  bool compare(int a_resID, char a_chainID, int a_seqNum, char a_iCode);
};

}; // namespace PDBParser

#endif
