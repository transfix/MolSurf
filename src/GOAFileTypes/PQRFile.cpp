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
#include <assert.h>
#include <Utility/utility.h>
// #include <Blurmaps/BlurMapsDataManager.h>
#include <GOAFileTypes/PQRFile.h>
#include <PDBParser/Atom.h>
#include <PDBParser/CollectionData.h>
#include <PDBParser/GroupOfAtoms.h>
#include <PDBParser/parserPDBtoGOA.h>
#include <PDBParser/ResidueHeader.h>
#include <PDBParser/ResidueInformation.h>
#include <PDBParser/Residues.h>
#include <PDBParser/FlattenGOA.h>
#include <sys/stat.h>
#include <sys/types.h>

PQRFile PQRFile::ms_PQRFileRepresentative;

PQRFile::PQRFile() {}

PQRFile::~PQRFile() {}

void PQRFile::setChainColors(PDBParser::GroupOfAtoms *chain, int n) {
  srand(3.14159 / 2.0 - n + 65 + 19.2); // why not ?
  double blue = rand() / (RAND_MAX + 1.0f);
  double red = rand() / (RAND_MAX + 1.0f);
  double green = rand() / (RAND_MAX + 1.0f);
  chain->setGOAColor(red, green, blue);
}

bool isAtomOrHetAtm(char *line) {
  if (!line) {
    return false;
  }
  if (strlen(line) < 6) {
    return false;
  }
  if (!((line[0] == 'A' && line[1] == 'T' && line[2] == 'O' &&
         line[3] == 'M') ||
        (line[0] == 'H' && line[1] == 'E' && line[2] == 'T' && line[3] == 'A' &&
         line[4] == 'T' && line[5] == 'M'))) {
    return false;
  }
  return true;
}

PDBParser::GroupOfAtoms *PQRFile::loadFile(const string &fileName,
                                           bool deleteWater) {
  FILE *fp = fopen(fileName.c_str(), "r");
  // check to make sure the file exists
  if (!fp) {
    printf("Error: could not open file %s\n", fileName.c_str());
    return 0;
  }
  PDBParser::GroupOfAtoms *molecule =
      PDBParser::GroupOfAtoms::getEmptyHierarchy();
  PDBParser::GroupOfAtoms *curProtein = molecule->m_SubGroups[0];
  PDBParser::GroupOfAtoms *curResidue =
      molecule->m_SubGroups[0]->m_SubGroups[0]->m_SubGroups[0]->m_SubGroups[0];
  PDBParser::GroupOfAtoms *curSS =
      molecule->m_SubGroups[0]->m_SubGroups[0]->m_SubGroups[0];
  molecule->type = PDBParser::COLLECTION_TYPE;
  molecule->m_CollectionData = new PDBParser::CollectionData();
  int numAtoms = 0;
  // problem: apbs deletes chain IDs, so we assign them.
  // We start with 'a'. Then we keep incrementing when the current residue
  // number becomes less than the prev residue number
  char cur_chainId = 'A';
  bool fileHasChainId = false;
  int cur_residueNum = -99999999;
  int numberOfResiduesInCurrentChain = 0;
  char line[82];
  memset(line, '\0',
         sizeof(line)); // Clear out the buffer in case of short lines
  while (fgets(line, 81, fp) != 0) {
    // this file format description is from APBS website.
    //
    // Field_name Atom_number Atom_name Residue_name Residue_number X Y Z Charge
    // Radius I have seen pqr files which do not follow this format ! SKVINAY
    //
    char ATOM_or_HETATM[256];
    int atomIndex;
    char atomName[256];
    char residueName[256];
    int residueNumber;
    float x, y, z, charge, radius;
    if (strlen(line) < 6) {
      continue;
    }
    if (!isAtomOrHetAtm(line)) {
      continue;
    }

    if (line[21] == ' ') {
      fileHasChainId = false;
      getString(line, ATOM_or_HETATM, 0, 6);
      getInt(line, &atomIndex, 6, 5);
      getString(line, atomName, 12, 4);
      getString(line, residueName, 17, 3);
      getInt(line, &residueNumber, 22, 4);
      getFloat(line, &x, 30, 8);
      getFloat(line, &y, 38, 8);
      getFloat(line, &z, 46, 8);
      getFloat(line, &charge, 54, 8);
      getFloat(line, &radius, 62, 7);
    } else // works only if there is a chainID column in the pqr file and the
           // pqr file is in a very standard format! --Wenqi
    {
      fileHasChainId = true;
      getString(line, ATOM_or_HETATM, 0, 6);
      getInt(line, &atomIndex, 6, 5);
      getString(line, atomName, 12, 4);
      getString(line, residueName, 17, 3);
      getChar(line, &cur_chainId, 21);
      getInt(line, &residueNumber, 22, 4);
      getFloat(line, &x, 30, 8);
      getFloat(line, &y, 38, 8);
      getFloat(line, &z, 46, 8);
      getFloat(line, &charge, 54, 8);
      getFloat(line, &radius, 62, 7);
    }

    // improve atom name to mimic pdb file format!
    {
      int len = strlen(atomName);
      if (len < 4) {
        char dup[256];
        strcpy(dup, atomName);
        sprintf(atomName, " %s", dup); // add leading space
      }
      len = strlen(atomName);
      if (len < 4) {
        int i;
        for (i = len; i < 4; i++) {
          atomName[i] = ' '; // add trailing spaces
        }
        atomName[4] = '\0';
      }
    }
    if (deleteWater) {
      if (residueName[0] == 'H' && residueName[1] == 'O' &&
          residueName[2] == 'H') {
        break;
      }
    }
    // check if new residue, first residue, new chain etc:
    int currResID = PDBParser::Residues::mapResidue(
        residueName); // The int identifier for the current residue
    // first residue has already been allocated, so use it.
    if (cur_residueNum == -99999999) {
      cur_residueNum = residueNumber;
      PDBParser::ResidueHeader *currHeader =
          new PDBParser::ResidueHeader(currResID, cur_chainId, 0, ' ');
      {
        PDBParser::parserPDBtoGOA *p = new PDBParser::parserPDBtoGOA();
        p->initResidue(curResidue, currHeader, currResID, atomIndex);
        delete p;
      }
    }
    {
      if (cur_residueNum == residueNumber) {
        // same residue
      } else {
        if (residueNumber < cur_residueNum) // not best condition, but ok, ==
                                            // condition is tough to add here
        {
          // we got a new chain.

          // add number of residues read for previous chain, if there was any!
          if (numberOfResiduesInCurrentChain != 0) {
            // arand, 9-6-2011: this list contains the chain id followed by the
            //                  number of residues...
            molecule->m_CollectionData->m_ChainIds.push_back(cur_chainId);

            molecule->m_CollectionData->m_ChainIds.push_back(
                numberOfResiduesInCurrentChain);
            numberOfResiduesInCurrentChain = 0;
          }
          // molecule->m_CollectionData->m_ChainIds.push_back((int)(cur_chainId));

          if (fileHasChainId == false) {
            cur_chainId++;
          }

          PDBParser::GroupOfAtoms *temp_chain = new PDBParser::GroupOfAtoms();

          temp_chain->resID = (int)cur_chainId;
          setChainColors(temp_chain, cur_chainId);

          temp_chain->type = PDBParser::CHAIN_TYPE;
          curProtein->m_SubGroups.push_back(temp_chain);

          if (currResID >= PDBParser::RES_TYPE_FIRST_NUCLEOTIDE) {
            temp_chain->resPos = PDBParser::CHAIN_NUCLEIC;
          } else {
            temp_chain->resPos = PDBParser::CHAIN_PROTEIN;
          }
          PDBParser::GroupOfAtoms *temp_ss = new PDBParser::GroupOfAtoms();
          temp_ss->type = PDBParser::SS_TYPE;
          temp_chain->m_SubGroups.push_back(temp_ss);
          curSS = temp_ss;
        }
        // we got a new residue
        PDBParser::GroupOfAtoms *temp_residue = new PDBParser::GroupOfAtoms();
        temp_residue->type = PDBParser::RESIDUE_TYPE;
        curSS->m_SubGroups.push_back(temp_residue);
        curResidue = temp_residue;
        numberOfResiduesInCurrentChain++;
        cur_residueNum = residueNumber;
        PDBParser::ResidueHeader *currHeader =
            new PDBParser::ResidueHeader(currResID, cur_chainId, 0, ' ');
        {
          PDBParser::parserPDBtoGOA *p = new PDBParser::parserPDBtoGOA();
          p->initResidue(curResidue, currHeader, currResID, atomIndex);
          delete p;
        }
        // arand, 9-6-2011: adding this info to match the PDB parser
        //                  ... needed for colormap editor...
        molecule->m_CollectionData->m_ResidueIds.push_back(
            (int)(PDBParser::Residues::lookupShortResID(residueName)));
        molecule->m_CollectionData->m_ResidueIds.push_back(cur_residueNum);
        molecule->m_CollectionData->m_ResidueIds.push_back(
            ' '); // not sure what to do with this...
      }
    }
    PDBParser::Atom *atom = new PDBParser::Atom();
    atom->m_Position[0] = x;
    atom->m_Position[1] = y;
    atom->m_Position[2] = z;
    {
      PDBParser::parserPDBtoGOA *p = new PDBParser::parserPDBtoGOA();
      // get the element index. If not found, get residue id and set it
      int elementIndex = p->mapPQRElementAndResidue(atomName, residueName);
      atom->setElementIndex(elementIndex);
      if (elementIndex == -1) {
        atom->setResidueId(PDBParser::Residues::mapResidue(residueName));
      }
      delete p;
    }
    atom->setRadius(radius);
    atom->setName(atomName);
    atom->setResidueNum(residueNumber);
    atom->setCharge(charge);
    atom->setChain(cur_chainId);
    curResidue->addAtom(
        atom->getAtomNum(), atom->getElementIndex(), atom->getResidueNum(),
        (char *)atom->getName(), atom->getResidueId(), atom->getChain(),
        atom->m_Position[0], atom->m_Position[1], atom->m_Position[2],
        atom->getRadius(), atom->getCharge(), false);
    delete atom;
    numAtoms++;
  }

  if (numberOfResiduesInCurrentChain != 0) {
    // arand, 9-6-2011:  don't forget the last chain...
    molecule->m_CollectionData->m_ChainIds.push_back(cur_chainId);

    molecule->m_CollectionData->m_ChainIds.push_back(
        numberOfResiduesInCurrentChain);
    numberOfResiduesInCurrentChain = 0;
  }

  double xmin = 0, ymin = 0, zmin = 0, xmax = 1, ymax = 1, zmax = 1;
  molecule->getMinMax(&xmin, &ymin, &zmin, &xmax, &ymax, &zmax);
  molecule->buildHierarchicalInformation(
      (xmin + xmax) / 2.0, (ymin + ymax) / 2.0, (zmin + zmax) / 2.0);
  {
    int count = 0;
    molecule->getNumberOfAtomsRecursive(&count);
  }
  if (numAtoms == 0) {
    delete molecule;
    molecule = 0;
  }
  fclose(fp);
  return molecule;
}

bool PQRFile::checkType(const string &fileName) { return false; }

bool PQRFile::saveFile(PDBParser::GroupOfAtoms *molecule,
                       const string &fileName, unsigned int level,
                       CCVOpenGLMath::Matrix *transformation) {
  FILE *fp = fopen(fileName.c_str(), "w");
  // failed to open the file? bail.
  if (!fp) {
    printf("Could not open file %s\n", fileName.c_str());
    return false;
  }
  // make 2 loops. First loop, get number of atoms. second loop, write them out.
  vector<PDBParser::Atom *> m_AtomList;
  PDBParser::CollectionData *collectionData = 0;
  if (molecule->type == PDBParser::COLLECTION_TYPE) {
    collectionData = molecule->m_CollectionData;
  }
  FlattenGOA(molecule, m_AtomList, collectionData, 0, 0, 0,
             PDBParser::GroupOfAtoms::VDW_RADIUS, level, false);
  int n = m_AtomList.size();
  {

    int i;
    for (i = 0; i < n; i++) {

      // arand: added transformation if available...
      double x, y, z;
      x = m_AtomList[i]->m_Position[0];
      y = m_AtomList[i]->m_Position[1];
      z = m_AtomList[i]->m_Position[2];
      if (transformation) {
        CCVOpenGLMath::Vector newPoint =
            (*transformation) * CCVOpenGLMath::Vector(x, y, z, 1);
        x = newPoint[0];
        y = newPoint[1];
        z = newPoint[2];
      }

      // Field_name Atom_number Atom_name Residue_name Residue_number X Y Z
      // Charge Radius
      //  ARAND: change to get colums in the right place?
      fprintf(fp, "ATOM %6d %4s %3s %5d    %8.3f%8.3f%8.3f %7.4f %6.4f\n",
              // fprintf(fp, "ATOM %6d %4s %3s %5d    %8.3f %8.3f %8.3f %7.4f
              // %6.4f\n",
              i + 1, m_AtomList[i]->getName(), m_AtomList[i]->getResidueName(),
              m_AtomList[i]->getResidueNum(), x, y, z,
              m_AtomList[i]->getCharge(), m_AtomList[i]->getRadius());
    }
  }
  // clean up
  fclose(fp);
  return true;
}

GOAFileType *PQRFile::getRepresentative() { return &ms_PQRFileRepresentative; }

void PQRFile::getString(char *line, char *str, int pos, int len) {
  char buffer[80];
  assert(len < 80);
  strcpy(str, "");
  strncpy(buffer, &line[pos], len);
  buffer[len] = '\0';
  for (int i = 0; i < len; i++) {
    str[i] = buffer[i];
  }
  str[len] = '\0';
}

void PQRFile::getChar(char *line, char *str, int pos) {
  char buffer[80];
  *str = ' ';
  strncpy(buffer, &line[pos], 1);
  buffer[1] = '\0';
  sscanf(buffer, "%c", str);
}

void PQRFile::getInt(char *line, int *i, int pos, int len) {
  char buffer[80];
  *i = -1;
  strncpy(buffer, &line[pos], len);
  buffer[len] = '\0';
  sscanf(buffer, "%d", i);
}

void PQRFile::getFloat(char *line, float *d, int pos, int len) {
  char buffer[80];
  *d = -1;
  strncpy(buffer, &line[pos], len);
  buffer[len] = '\0';
  sscanf(buffer, "%f", d);
}
