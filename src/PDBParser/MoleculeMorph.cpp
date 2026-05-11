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
#include <PDBParser/MoleculeMorph.h>

using namespace PDBParser;

MoleculeMorph::MoleculeMorph() {
  m_Molecule1 = 0;
  m_Molecule2 = 0;
  m_OutputFileNamePrefix = 0;
  m_Resolution = 2;
}

MoleculeMorph::~MoleculeMorph() {
  m_Molecule1 = 0;
  m_Molecule2 = 0;
  m_OutputFileNamePrefix = 0;
  m_Resolution = 2;
}

double MoleculeMorph::getSize() {
  double size;
  {
    double minx, miny, minz, maxx, maxy, maxz;
    m_Molecule1->getMinMax(&minx, &miny, &minz, &maxx, &maxy, &maxz);
    double dimx = maxx - minx;
    double dimy = maxy - miny;
    double dimz = maxz - minz;
    if (dimx > dimy) {
      size = dimx;
    } else {
      size = dimy;
    }
    if (dimz > size) {
      size = dimz;
    }
  }
  {
    double minx, miny, minz, maxx, maxy, maxz;
    m_Molecule2->getMinMax(&minx, &miny, &minz, &maxx, &maxy, &maxz);
    double dimx = maxx - minx;
    double dimy = maxy - miny;
    double dimz = maxz - minz;
    if (dimx > size) {
      size = dimx;
    }
    if (dimy > size) {
      size = dimy;
    }
    if (dimz > size) {
      size = dimz;
    }
  }
  return size;
}

bool MoleculeMorph::torsionAnglesPresent(PDBParser::GroupOfAtoms *molecule) {
  if (!molecule) {
    return false;
  }
  if (molecule->type == RESIDUE_TYPE) {
    if (!molecule->m_ResidueData) {
      return false;
    }
    return true;
  }
  int i, n;
  // Iterate through the subgroups array and descend down the tree.
  n = (int)molecule->m_SubGroups.size();
  for (i = 0; i < n; i++) {
    GroupOfAtoms *currSubGroup = molecule->m_SubGroups.at(i);
    if (!torsionAnglesPresent(currSubGroup)) {
      return false;
    }
  }
  return true;
}

//  Works only for AA - SKVINAY
bool MoleculeMorph::getTorsionAngles(PDBParser::GroupOfAtoms *molecule,
                                     double *torsionAngles, int *count) {
  if (!molecule || !torsionAngles || !count) {
    return false;
  }
  if (molecule->type == RESIDUE_TYPE) {
    if (!molecule->m_ResidueData) {
      return false;
    }
    torsionAngles[(*count) * 2 + 0] = molecule->m_ResidueData->m_Phi;
    torsionAngles[(*count) * 2 + 1] = molecule->m_ResidueData->m_Psi;
    (*count)++;
    return true;
  }
  int i, n;
  // Iterate through the subgroups array and descend down the tree.
  n = (int)molecule->m_SubGroups.size();
  for (i = 0; i < n; i++) {
    GroupOfAtoms *currSubGroup = molecule->m_SubGroups.at(i);
    if (!getTorsionAngles(currSubGroup, torsionAngles, count)) {
      return false;
    }
  }
  return true;
}

bool MoleculeMorph::getResidueSequenceIDs(PDBParser::GroupOfAtoms *molecule,
                                          int *residueSequenceNumbers,
                                          int *count) {
  if (!molecule || !residueSequenceNumbers || !count) {
    return false;
  }
  if (molecule->type == RESIDUE_TYPE) {
    residueSequenceNumbers[(*count)] = molecule->m_ResSeq;
    (*count)++;
    return true;
  }
  int i, n;
  // Iterate through the subgroups array and descend down the tree.
  n = (int)molecule->m_SubGroups.size();
  for (i = 0; i < n; i++) {
    GroupOfAtoms *currSubGroup = molecule->m_SubGroups.at(i);
    if (!getResidueSequenceIDs(currSubGroup, residueSequenceNumbers, count)) {
      return false;
    }
  }
  return true;
}

//  Morphs one molecule to another. If resolution is h, then
//  we need 2 pi r / h rotation steps.
//
//  We do not check if residues are the same. As long as there
//  are the same number, we are ok.
bool MoleculeMorph::morph(PDBParser::GroupOfAtoms *molecule1,
                          PDBParser::GroupOfAtoms *molecule2,
                          const char *outputFileNamePrefix, double resolution,
                          int maxSteps) {
  if (!molecule1 || !molecule2 || !outputFileNamePrefix ||
      (resolution < 0.0001)) {
    return false;
  }
  m_Molecule1 = molecule1;
  m_Molecule2 = molecule2;
  m_OutputFileNamePrefix = outputFileNamePrefix;
  m_Resolution = resolution;
  if (!torsionAnglesPresent(molecule1)) {
    return false;
  }
  if (!torsionAnglesPresent(molecule2)) {
    return false;
  }
  int numRes1 = molecule1->getNumberOfGOAs(PDBParser::GroupOfAtoms::RESIDUE);
  int numRes2 = molecule1->getNumberOfGOAs(PDBParser::GroupOfAtoms::RESIDUE);
  if (numRes1 != numRes2) {
    return false;
  }
  if (numRes1 < 1) {
    return false;
  }
  // n = getNumberOfTorsionAngles();
  int n = numRes1 * 2; // SKVINAY works only for AA
  double *torsionAngles1 = new double[n];
  double *torsionAngles2 = new double[n];
  int *residueSequenceNumbers = new int[numRes1];
  int count = 0;
  if (!getTorsionAngles(molecule1, torsionAngles1, &count)) {
    delete[] torsionAngles1;
    delete[] torsionAngles2;
    delete[] residueSequenceNumbers;
    return false;
  }
  count = 0;
  if (!getTorsionAngles(molecule2, torsionAngles2, &count)) {
    delete[] torsionAngles1;
    delete[] torsionAngles2;
    delete[] residueSequenceNumbers;
    return false;
  }
  count = 0;
  if (!getResidueSequenceIDs(molecule1, residueSequenceNumbers, &count)) {
    delete[] torsionAngles1;
    delete[] torsionAngles2;
    delete[] residueSequenceNumbers;
    return false;
  }
  int numModels = 0;
  {
    double size = getSize();
    double radius = size / 2.0;
    numModels = 2 * 3.14159 * radius / resolution;
    if (numModels > maxSteps) {
      numModels = maxSteps;
    }
    if (numModels < 2) {
      return false;
    }
  }
  // write out evole scripts! this is a shortcut. no time - SKVINAY
  {
    FILE *fpBatch = fopen("../Dataset/flexDock/evolve.bat", "w");
    if (!fpBatch) {
      return false;
    }
    int i;
    for (i = 1; i < numModels - 1; i++) // first and last are already there
    {
      char goaFilename[256];
      char evolveFileName[256];
      sprintf(goaFilename, "%s%06d.pdb", outputFileNamePrefix, i);
      sprintf(evolveFileName, "../Dataset/flexDock/backboneRotations%06d.txt",
              i);
      FILE *fp = fopen(evolveFileName, "w");
      if (!fp) {
        return false;
      }
      fprintf(fpBatch,
              "..\\..\\MoleculeViz\\Debug\\MoleculeViz -evolve 1CLL.pdb "
              "backboneRotations%06d.txt %s\n",
              i, goaFilename);
      int j;
      for (j = 0; j < numRes1; j++) {
        double new_phi =
            torsionAngles1[2 * j + 0] +
            (torsionAngles2[2 * j + 0] - torsionAngles1[2 * j + 0]) * i /
                (numModels - 1);
        double new_psi =
            torsionAngles1[2 * j + 1] +
            (torsionAngles2[2 * j + 1] - torsionAngles1[2 * j + 1]) * i /
                (numModels - 1);
        // <0> <sequence number> <new phi angle> <new psi angle>
        fprintf(fp, "0 %d %f %f\n", residueSequenceNumbers[j], new_phi,
                new_psi);
      }
      fclose(fp);
    }
    fclose(fpBatch);
  }
  return true;
}
