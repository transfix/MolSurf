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
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "LeastRMSD/LeastRMSD.h"
#include <PDBParser/FlattenGOA.h>
#include <GOAFileTypes/GOALoader.h>

using namespace std;

// arand: added MolSurf interface... 3-30-2011

void usageLeastRMSD() {
  cout << endl
       << "Usage: MolSurf -leastRMSD <pdb 1> <pdb 2>" << endl
       << "    Compute RMSD." << endl;
}

bool getLeastRMSD(int argc, char **argv) {
  if (argc != 4) {
    usageLeastRMSD();
    return false;
  }

  // read PDB from f1
  PDBParser::GroupOfAtoms *molecule1 = GOALoader().loadFile(argv[2]);
  if (!molecule1) {
    exit(-1);
  }
  vector<PDBParser::Atom *> atomList1;
  FlattenGOA(molecule1, atomList1, 0, 0, 0, 0,
             PDBParser::GroupOfAtoms::VDW_RADIUS, PDBParser::ATOM_TYPE, false);

  // read PDB from f2
  PDBParser::GroupOfAtoms *molecule2 = GOALoader().loadFile(argv[3]);
  if (!molecule2) {
    exit(-1);
  }
  vector<PDBParser::Atom *> atomList2;
  FlattenGOA(molecule2, atomList2, 0, 0, 0, 0,
             PDBParser::GroupOfAtoms::VDW_RADIUS, PDBParser::ATOM_TYPE, false);

  // convert PDB to vector< CCVOpenGLMath::Vector >
  vector<CCVOpenGLMath::Vector> pts1;
  for (int i = 0; i < atomList1.size(); i++) {
    float *tmp = atomList1[i]->getPosition();
    CCVOpenGLMath::Vector v(tmp[0], tmp[1], tmp[2], 0.0);
    pts1.push_back(v);
  }

  vector<CCVOpenGLMath::Vector> pts2;
  for (int i = 0; i < atomList2.size(); i++) {
    float *tmp = atomList2[i]->getPosition();
    CCVOpenGLMath::Vector v(tmp[0], tmp[1], tmp[2], 0.0);
    pts2.push_back(v);
  }

  // run the code
  CCVOpenGLMath::LeastRMSD L;
  double rmsd = L.getLeastRMSD(pts1, pts2);

  cout << "RMSD: " << rmsd << endl;

  return true;
}
