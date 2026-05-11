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
#include <PDBParser/ResidueHeader.h>

using PDBParser::ResidueHeader;

ResidueHeader::ResidueHeader() {
  // Initialize to values that indicate uninitialized status
  resID = -1;
  chainID = '\0';
  seqNum = -999999; // Valid sequence numbers are only 4 digits long
  iCode = '\0';
  // ADD to copy constructer also !
}

ResidueHeader::ResidueHeader(int a_resID, char a_chainID, int a_seqNum,
                             char a_iCode) {
  resID = a_resID;
  chainID = a_chainID;
  seqNum = a_seqNum;
  iCode = a_iCode;
}

ResidueHeader::ResidueHeader(const ResidueHeader &rh) {
  resID = rh.resID;
  chainID = rh.chainID;
  seqNum = rh.seqNum;
  iCode = rh.iCode;
}

ResidueHeader::~ResidueHeader() {}

// Comparing two residue headers.  This is done kind of like Java-style string
// comparison, where you call this method in ResidueHeader 1 using ResidueHeader
// 2 as an arg in order to compare ResidueHeaders 1 and 2.  Returns true if the
// two headers are identical, false otherwise.
bool ResidueHeader::compare(ResidueHeader *rh) {
  if (resID == rh->resID && chainID == rh->chainID && seqNum == rh->seqNum &&
      iCode == rh->iCode) {
    return true;
  }
  return false;
}

// This compare function takes the values directly instead of in a residue
// header. Returns true if the two headers are identical, false otherwise.
bool ResidueHeader::compare(int a_resID, char a_chainID, int a_seqNum,
                            char a_iCode) {
  if (resID == a_resID && chainID == a_chainID && seqNum == a_seqNum &&
      iCode == a_iCode) {
    return true;
  }
  return false;
}
