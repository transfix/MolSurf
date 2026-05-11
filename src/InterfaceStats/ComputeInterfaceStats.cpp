/*
  Copyright 2011 The University of Texas at Austin

        Authors: Muhibur Rasheed <muhibur@ices.utexas.edu>
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
#include <InterfaceStats/ComputeInterfaceStats.h>
#include <PDBParser/GroupOfAtoms.h>
#include <PDBParser/FlattenGOA.h>
#include <GOAFileTypes/GOALoader.h>
#include <iostream>
#include <fstream>

ComputeInterfaceStats::ComputeInterfaceStats(string pqrFileName1,
                                             string pqrFileName2,
                                             string quadFileName1,
                                             string quadFileName2,
                                             vector<Matrix> xformations) {
  cout << "Trying to compute Stats: in ComputeInterfaceStats now" << endl;
  transformations = xformations;
  numInterfaces = transformations.size();

  cout << "PQR file 1: " << pqrFileName1 << endl;
  cout << "PQR file 2: " << pqrFileName2 << endl;
  cout << "Quad file 1: " << quadFileName1 << endl;
  cout << "Quad file 2: " << quadFileName2 << endl;

  readAtomsFromFiles(pqrFileName1, pqrFileName2);

  readQuadPointsFromFiles(quadFileName1, quadFileName2);

  cout << "files read" << endl;

  computeAllStats();
}

ComputeInterfaceStats::~ComputeInterfaceStats() {
  atomList1.clear();
  atomList2.clear();

  qPoints1.clear();
  qPoints2.clear();

  transformations.clear();
  interfaces.clear();
}

void ComputeInterfaceStats::computeAllStats() {
  atomInterface = new AtomsOnInterface(atomList1, atomList2);
  pointInterface = new PointsOnInterface(qPoints1, qPoints2);

  cout << "Trying to compute Stats: octrees built" << endl;

  for (int i = 0; i < numInterfaces; i++) {
    cout << i << endl;

    InterfaceStats *intStat = new InterfaceStats();

    if (!atomInterface->computeInteractions(transformations[i], intStat)) {
      cout << "Could not compute atomic stat for transformation " << i + 1
           << endl;
    }

    cout << "   atomic stats ready" << endl;

    if (!pointInterface->computeInteractions(transformations[i], intStat)) {
      cout << "Could not compute surface stat for transformation " << i + 1
           << endl;
    }

    cout << "   surface stats ready" << endl;

    interfaces.push_back(intStat);

    cout << endl
         << endl
         << "Printing all stats for xform " << i + 1 << endl
         << endl;
    intStat->print();
    cout << endl << endl;
  }

  cout << "Stats computed" << endl;
}

void ComputeInterfaceStats::readAtomsFromFiles(string receptorFileName,
                                               string ligandfileName) {
  PDBParser::GroupOfAtoms *molecule1 = 0;
  GOALoader *gLoader1 = new GOALoader();
  molecule1 = gLoader1->loadFile(receptorFileName.c_str());
  delete gLoader1;
  if (!molecule1) {
    return;
  }

  PDBParser::GroupOfAtoms *molecule2 = 0;
  GOALoader *gLoader2 = new GOALoader();
  molecule2 = gLoader2->loadFile(ligandfileName.c_str());
  delete gLoader2;
  if (!molecule2) {
    return;
  }

  PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType;
  PDBParser::GroupOfAtoms::intToRadiusType(&radiusType, 0);

  PDBParser::CollectionData *collectionData1 = 0;
  if (molecule1->type == PDBParser::COLLECTION_TYPE) {
    collectionData1 = molecule1->m_CollectionData;
  }
  PDBParser::FlattenGOA(molecule1, atomList1, collectionData1, 0, 0, 0,
                        radiusType, PDBParser::ATOM_TYPE, false);

  PDBParser::CollectionData *collectionData2 = 0;
  if (molecule2->type == PDBParser::COLLECTION_TYPE) {
    collectionData2 = molecule2->m_CollectionData;
  }
  PDBParser::FlattenGOA(molecule2, atomList2, collectionData2, 0, 0, 0,
                        radiusType, PDBParser::ATOM_TYPE, false);
}

void ComputeInterfaceStats::readQuadPointsFromFiles(string receptorFileName,
                                                    string ligandfileName) {
  std::ifstream infile;
  double x, y, z;
  double nx, ny, nz;
  double w;

  infile.open(receptorFileName.c_str());

  if (!infile.is_open()) {
    cout << "Failed to open quadrature points file (%s)! " << receptorFileName
         << endl;
    return;
  }

  while (!infile.eof()) {
    infile >> x;
    infile >> y;
    infile >> z;
    infile >> nx;
    infile >> ny;
    infile >> nz;
    infile >> w;

    DPG::WeightedPoint *wp = new DPG::WeightedPoint(x, y, z, w);

    qPoints1.push_back(wp);
  }

  infile.close();

  infile.open(ligandfileName.c_str());

  if (!infile.is_open()) {
    cout << "Failed to open quadrature points file (%s)! " << ligandfileName
         << endl;
    return;
  }

  while (!infile.eof()) {
    infile >> x;
    infile >> y;
    infile >> z;
    infile >> nx;
    infile >> ny;
    infile >> nz;
    infile >> w;

    DPG::WeightedPoint *wp = new DPG::WeightedPoint(x, y, z, w);

    qPoints2.push_back(wp);
  }

  infile.close();
}
