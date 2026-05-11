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
#include <InterfaceStats/CapriCriteria.h>

CapriCriteria::CapriCriteria(string rf, string lf1, string lf2) {
  receptorFileName = rf;
  correctLigandFileName = lf1;
  predictedLigandFileName = lf2;
  mode = true;

  if (!readAtomsFromFiles()) {
    cout << "Error: could not read molecules" << endl;
    state = false;
  } else {
    state = true;
  }

  identity.set(0, 0, 1);
  identity.set(0, 1, 0);
  identity.set(0, 2, 0);
  identity.set(0, 3, 0);

  identity.set(1, 0, 0);
  identity.set(1, 1, 1);
  identity.set(1, 2, 0);
  identity.set(1, 3, 0);

  identity.set(2, 0, 0);
  identity.set(2, 1, 0);
  identity.set(2, 2, 1);
  identity.set(2, 3, 0);

  identity.set(3, 0, 0);
  identity.set(3, 1, 0);
  identity.set(3, 2, 0);
  identity.set(3, 3, 1);

  init();
}

CapriCriteria::CapriCriteria(string rf, string lf1, Matrix m) {
  receptorFileName = rf;
  correctLigandFileName = lf1;
  mtx = m;
  mode = false;

  if (!readAtomsFromFiles()) {
    cout << "Error: could not read molecules" << endl;
    state = false;
  } else {
    state = true;
  }

  identity.set(0, 0, 1);
  identity.set(0, 1, 0);
  identity.set(0, 2, 0);
  identity.set(0, 3, 0);

  identity.set(1, 0, 0);
  identity.set(1, 1, 1);
  identity.set(1, 2, 0);
  identity.set(1, 3, 0);

  identity.set(2, 0, 0);
  identity.set(2, 1, 0);
  identity.set(2, 2, 1);
  identity.set(2, 3, 0);

  identity.set(3, 0, 0);
  identity.set(3, 1, 0);
  identity.set(3, 2, 0);
  identity.set(3, 3, 1);

  init();
}

CapriCriteria::CapriCriteria(string rf, string lf1) {
  receptorFileName = rf;
  correctLigandFileName = lf1;
  mode = false;

  if (!readAtomsFromFiles()) {
    cout << "Error: could not read molecules" << endl;
    state = false;
  } else {
    state = true;
  }

  identity.set(0, 0, 1);
  identity.set(0, 1, 0);
  identity.set(0, 2, 0);
  identity.set(0, 3, 0);

  identity.set(1, 0, 0);
  identity.set(1, 1, 1);
  identity.set(1, 2, 0);
  identity.set(1, 3, 0);

  identity.set(2, 0, 0);
  identity.set(2, 1, 0);
  identity.set(2, 2, 1);
  identity.set(2, 3, 0);

  identity.set(3, 0, 0);
  identity.set(3, 1, 0);
  identity.set(3, 2, 0);
  identity.set(3, 3, 1);

  mtx = identity;

  init();
}

CapriCriteria::~CapriCriteria() {
  receptorAtoms.clear();
  correctLigandAtoms.clear();
  predictedLigandAtoms.clear();
  correctResiduePairs.clear();
  predictedResiduePairs.clear();
  interfaceAtoms.clear();
}

void CapriCriteria::init() {
  InterfaceAtomDetector *interfaceAtomDetector = new InterfaceAtomDetector(
      receptor, correctLigand, "temp.xyz", 10.0, "backbone", radiusType);
  interfaceAtomDetector->detect();
  interfaceAtoms = interfaceAtomDetector->getInterfaceAtoms();
  cout << "Ligand has " << interfaceAtoms.size()
       << " backbone atoms on the interface " << endl;

  getCorrectResiduePairs();
}

bool CapriCriteria::readAtomsFromFiles() {
  GOALoader *gLoader1 = new GOALoader();
  receptor = gLoader1->loadFile(receptorFileName.c_str());
  delete gLoader1;
  if (!receptor) {
    return false;
  }

  GOALoader *gLoader2 = new GOALoader();
  correctLigand = gLoader2->loadFile(correctLigandFileName.c_str());
  delete gLoader2;
  if (!correctLigand) {
    return false;
  }

  if (mode) {
    GOALoader *gLoader3 = new GOALoader();
    predictedLigand = gLoader3->loadFile(predictedLigandFileName.c_str());
    delete gLoader3;
    if (!predictedLigand) {
      return false;
    }
  }

  PDBParser::GroupOfAtoms::intToRadiusType(&radiusType, 0);

  PDBParser::CollectionData *collectionData1 = 0;
  if (receptor->type == PDBParser::COLLECTION_TYPE) {
    collectionData1 = receptor->m_CollectionData;
  }
  PDBParser::FlattenGOA(receptor, receptorAtoms, collectionData1, 0, 0, 0,
                        radiusType, PDBParser::ATOM_TYPE, false);

  cout << "Receptor has " << receptorAtoms.size() << " atoms" << endl;

  PDBParser::CollectionData *collectionData2 = 0;
  if (correctLigand->type == PDBParser::COLLECTION_TYPE) {
    collectionData2 = correctLigand->m_CollectionData;
  }
  PDBParser::FlattenGOA(correctLigand, correctLigandAtoms, collectionData2, 0,
                        0, 0, radiusType, PDBParser::ATOM_TYPE, false);

  cout << "Ligand has " << correctLigandAtoms.size() << " atoms" << endl;

  if (mode) {
    PDBParser::CollectionData *collectionData3 = 0;
    if (predictedLigand->type == PDBParser::COLLECTION_TYPE) {
      collectionData3 = predictedLigand->m_CollectionData;
    }
    PDBParser::FlattenGOA(predictedLigand, predictedLigandAtoms,
                          collectionData3, 0, 0, 0, radiusType,
                          PDBParser::ATOM_TYPE, false);

    cout << "Predicted Ligand has " << predictedLigandAtoms.size() << " atoms"
         << endl;
  } else {
    predictedLigandAtoms = correctLigandAtoms;
  }

  return true;
}

bool CapriCriteria::isBackBoneAtom(PDBParser::Atom *atom) {
  string name(atom->getName());

  return (name.compare(" N  ") == 0) || (name.compare(" O  ") == 0) ||
         (name.compare(" C  ") == 0) || (name.compare(" CA ") == 0);
}

bool CapriCriteria::getCorrectResiduePairs() {
  AtomsOnInterface *atomInterface =
      new AtomsOnInterface(receptorAtoms, correctLigandAtoms);
  InterfaceStats *intStat = new InterfaceStats();

  if (!atomInterface->computeInteractions(identity, intStat)) {
    cout << "Could not find residue contacts" << endl;
    return false;
  }

  correctResiduePairs = atomInterface->getAllPairs();

  cout << "Correct Ligand has " << correctResiduePairs.size()
       << " pairs of interface residues" << endl;

  return true;
}

bool CapriCriteria::getPredictedResiduePairs() {
  AtomsOnInterface *atomInterface =
      new AtomsOnInterface(receptorAtoms, predictedLigandAtoms);
  InterfaceStats *intStat = new InterfaceStats();
  predictedResiduePairs.clear();

  if (mode) {
    if (!atomInterface->computeInteractions(identity, intStat)) {
      cout << "Could not find residue contacts" << endl;
      return false;
    }
  } else {
    if (!atomInterface->computeInteractions(mtx, intStat)) {
      cout << "Could not find residue contacts" << endl;
      return false;
    }
  }

  predictedResiduePairs = atomInterface->getAllPairs();

  cout << "Predicted Ligand has " << predictedResiduePairs.size()
       << " pairs of interface residues" << endl;

  return true;
}

bool CapriCriteria::computeIRMSD() {
  double distsq;
  double sumdistsq = 0.0;
  int numInterfaceAtoms = interfaceAtoms.size();
  double x, y, z;
  double dx, dy, dz;

  for (int i = 0; i < numInterfaceAtoms; i++) {
    PDBParser::Atom *a1 = correctLigandAtoms[interfaceAtoms[i]];
    PDBParser::Atom *a2 = predictedLigandAtoms[interfaceAtoms[i]];

    if (mode) {
      x = a2->getX();
      y = a2->getY();
      z = a2->getZ();
    } else {
      Vector oldPos(a2->getX(), a2->getY(), a2->getZ(), 1.0);
      Vector newPos = mtx * oldPos;
      x = newPos[0];
      y = newPos[1];
      z = newPos[2];
    }

    dx = x - a1->getX();
    dy = y - a1->getY();
    dz = z - a1->getZ();
    distsq = dx * dx + dy * dy + dz * dz;

    sumdistsq += distsq;
  }
  irmsd = sqrt(sumdistsq / (double)numInterfaceAtoms);

  cout << "I-RMSD is " << irmsd << endl;

  return true;
}

bool CapriCriteria::computeLRMSD() {
  vector<int> correctLigandBackBoneAtomIndex;
  vector<int> predictedLigandBackBoneAtomIndex;

  for (int i = 0; i < correctLigandAtoms.size(); i++)
    if (isBackBoneAtom(correctLigandAtoms[i]))
      correctLigandBackBoneAtomIndex.push_back(i);

  for (int i = 0; i < predictedLigandAtoms.size(); i++)
    if (isBackBoneAtom(predictedLigandAtoms[i]))
      predictedLigandBackBoneAtomIndex.push_back(i);

  int numCorrectBackBoneAtoms = correctLigandBackBoneAtomIndex.size();
  int numPredictedBackBoneAtoms = predictedLigandBackBoneAtomIndex.size();

  cout << "Ligand has " << numCorrectBackBoneAtoms
       << " backbone atoms in total " << endl;

  if (numCorrectBackBoneAtoms != numPredictedBackBoneAtoms) {
    cout << "Error: correct and predicted ligands must have same number of "
            "backbone atoms"
         << endl;
    return false;
  }

  else {
    double distsq;
    double sumdistsq = 0.0;
    double x, y, z;
    double dx, dy, dz;

    for (int i = 0; i < numCorrectBackBoneAtoms; i++) {
      int correctBBAtomIndex = correctLigandBackBoneAtomIndex[i];
      int predictedBBAtomIndex = predictedLigandBackBoneAtomIndex[i];

      PDBParser::Atom *a1 = correctLigandAtoms[correctBBAtomIndex];
      PDBParser::Atom *a2 = predictedLigandAtoms[predictedBBAtomIndex];

      DPG::Point *temp;

      if (mode) {
        x = a2->getX();
        y = a2->getY();
        z = a2->getZ();
      } else {
        Vector oldPos(a2->getX(), a2->getY(), a2->getZ(), 1.0);
        Vector newPos = mtx * oldPos;
        x = newPos[0];
        y = newPos[1];
        z = newPos[2];
      }

      dx = x - a1->getX();
      dy = y - a1->getY();
      dz = z - a1->getZ();
      distsq = dx * dx + dy * dy + dz * dz;

      sumdistsq += distsq;
    }

    lrmsd = sqrt(sumdistsq / (double)numCorrectBackBoneAtoms);
  }

  cout << "L-RMSD is " << lrmsd << endl;

  return true;
}

bool CapriCriteria::computeNats() {
  numNat = 0;
  numPredMinusNat = 0;

  if (!state)
    return state;

  //      if(!getCorrectResiduePairs())
  //              return false;

  if (!getPredictedResiduePairs())
    return false;

  numNat = correctResiduePairs.size();
  int numPred = predictedResiduePairs.size();
  int numCommon = 0;

  for (int i = 0; i < numNat; i++) {
    int aNatPair = correctResiduePairs[i];

    for (int j = 0; j < numPred; j++) {
      int aPredPair = predictedResiduePairs[j];
      if (aNatPair == aPredPair)
        numCommon++;
    }
  }

  numPredMinusNat = numPred - numCommon;

  fnat = (double)numCommon / (double)numNat;
  fnonnat = (double)numPredMinusNat / (double)numNat;

  cout << "fnat " << fnat << endl;
  cout << "fnonnat " << fnonnat << endl;

  return true;
}

bool CapriCriteria::computeCapriScore() {
  if (!state)
    return state;

  if (!computeIRMSD())
    return false;
  if (!computeLRMSD())
    return false;
  if (!computeNats())
    return false;

  verdict = 0;

  if (fnat < 0.1) {
    verdict = 0; // INCORRECT
  } else if (fnat < 0.3) {
    if (irmsd <= 4.0 || lrmsd <= 10.0)
      verdict = 1; // ACCEPTABLE
    else
      verdict = 0;
  } else if (fnat < 0.5) {
    if (irmsd <= 2.0 || lrmsd <= 5.0)
      verdict = 2; // MEDIUM
    else if (irmsd <= 4.0 || lrmsd <= 10.0)
      verdict = 1; // ACCEPTABLE
    else
      verdict = 0; // INCORRECT
  } else {
    if (irmsd <= 1.0 || lrmsd <= 1.0)
      verdict = 3; // HIGH
    else if (irmsd <= 2.0 || lrmsd <= 5.0)
      verdict = 2; // MEDIUM
    else if (irmsd <= 4.0 || lrmsd <= 10.0)
      verdict = 1; // ACCEPTABLE
    else
      verdict = 0; // INCORRECT
  }

  return true;
}

bool CapriCriteria::computeCapriScore(Matrix m) {
  mtx = m;
  mode = false;
  return computeCapriScore();
}
