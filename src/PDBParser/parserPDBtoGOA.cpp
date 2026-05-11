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
#include <PDBParser/parserPDBtoGOA.h>
// For circular definitions
#include <PDBParser/Residues.h>
#include <PDBParser/GOAVector.h>
#include <PDBParser/Atom.h>
#include <PDBParser/GroupOfAtoms.h>
#include <PDBParser/GOAFileIO.h>
#include <PDBParser/GOAMolGeometry.h>

using namespace CVCUtility;
using namespace PDBParser;

// Initialize member data
parserPDBtoGOA::parserPDBtoGOA() {
  int maxSize = 1000000;
  xmin = maxSize;
  ymin = maxSize;
  zmin = maxSize;
  xmax = -maxSize;
  ymax = -maxSize;
  zmax = -maxSize;
  residueArray = 0;
  residueAssigned = 0;
  residueCheckCount = 0;
  chainIDLast = '\0';
  residueSeqLast = -999999; // Sequence numbers can be negative, but only 4
                            // digits, so -999999 means uninitialized
  iCodeLast = '\0';
  currResidue = 0;
  currSS = 0;
  currChain = 0;
  helixGroupNext = 0;
  sheetGroupNext = 0;
  turnGroupNext = 0;
  sasThreshold = 1;
  m_FileType = UNKNOWN_FILE;
  memset(chainTable, -1, sizeof(chainTable));
  m_NumberOfResiduesInCurrentChain = 0;
}

parserPDBtoGOA::~parserPDBtoGOA() {}

// Craig: these four should be moved to util or should just use std lib tools
void parserPDBtoGOA::getString(char *line, char *str, int pos, int len) {
  /*
          char buffer[80];
          assert(len < 80);
          strcpy(str, "");
          strncpy(buffer, &line[pos], len);
          buffer[len] = '\0';
          for(int i=0; i<len; i++)
          {
                  str[i] = buffer[i];
          }
          str[len] = '\0';
  */
  // CHA: make this straightforward!
  // Assume, pos + len is not exceed length of line array.
  for (int i = 0; i < len; i++)
    str[i] = line[pos + i];
  str[len] = '\0';
}

void parserPDBtoGOA::getChar(char *line, char *str, int pos) {
  /*
          char buffer[80];
          *str = ' ';
          strncpy(buffer, &line[pos], 1);
          buffer[1] = '\0';
          sscanf(buffer, "%c", str);
  */
  // CHA: make this straightforward!
  // Assume, pos is not exceed length of line array.
  *str = line[pos];
}

void parserPDBtoGOA::getInt(char *line, int *i, int pos, int len) {
  /*
  char buffer[80];
  *i = -1;
  strncpy(buffer, &line[pos], len);
  buffer[len] = '\0';
  sscanf(buffer, "%d", i);
  */
  // CHA: make this straightforward! using standard string function.
  // Assume, pos is not exceed length of line array.
  char buf[80];
  assert(len < 80);
  for (int k = 0; k < len; k++)
    buf[k] = line[pos + k];
  buf[len] = '\0';
  *i = atoi(buf);
}

void parserPDBtoGOA::getDouble(char *line, double *d, int pos, int len) {
  /*
          char buffer[80];
          *d = -1;
          strncpy(buffer, &line[pos], len);
          buffer[len] = '\0';
          sscanf(buffer, "%lf", d);
  */
  // CHA: make this straightforward! using standard string function.
  // Assume, pos is not exceed length of line array.
  char buf[80];
  assert(len < 80);
  for (int k = 0; k < len; k++)
    buf[k] = line[pos + k];
  buf[len] = '\0';
  *d = atof(buf);
}

char *parserPDBtoGOA::trim(char *s) {
  int i;

  while (isspace(*s))
    s++; // skip left side white spaces
  for (i = strlen(s) - 1; (isspace(s[i])); i--)
    ; // skip right side white spaces
  s[i + 1] = '\0';
  return s;
}

// Given the current sequence number and inserion code,
// returns true if they are different from the last sequence number insertion
// code.
bool parserPDBtoGOA::residueHasChanged(char chainID, int residueSeq,
                                       char iCode) {
  return (chainIDLast != chainID || residueSeqLast != residueSeq ||
          iCodeLast != iCode);
}

// Allocate new space for a molecule if space has not already been allocated.
// The one case where space might already be allocated is for the first MODEL in
// a file, since we always allocate a default GOA molecule on startup (for files
// that have no MODELs).
void parserPDBtoGOA::parseMODEL(GroupOfAtoms *molecule, char *line) {
  /*      1 -  6       Record name    "MODEL "
          11 - 14       Integer        serial        Model serial number.*/
  char model[10];
  getString(line, model, 0, 6);
  int serial;
  getInt(line, &serial, 10, 4);
  if (!molecule) {
    // Create a new molecule to be filled in (note the global reference is used
    // here)
    currMolecule = new GroupOfAtoms();
    // debug: what about nucleotides, are they GENERAL_TYPE?
    currMolecule->type = PROTEIN_TYPE;
    currMolecule->setMoleculeColor();
  }
}

// Finalize the existing molecule
void parserPDBtoGOA::parseENDMDL(GroupOfAtoms *molecule, char *line) {
  /*      1 -  6         Record name      "ENDMDL"*/
  char endmdl[10];
  getString(line, endmdl, 0, 6);
  if (molecule) {
    finalizeParsing(molecule);
    molModels.add(molecule);
  } else {
    // Something is fishy if there is no prior molecule
    error("parseENDMDL: no previous molecule");
  }
  // Clear out the current molecule (note the global reference is used here)
  currMolecule = NULL;
}

void parserPDBtoGOA::parseSEQADV(GroupOfAtoms *molecule, char *line) {
  /*      1 -  6        Record name     "SEQADV"
           8 - 11        IDcode          idCode         ID code of this entry.
          13 - 15        Residue name    resName        Name of the PDB residue
     in conflict. 17             Character       chainID        PDB chain
     identifier. 19 - 22        Integer         seqNum         PDB sequence
     number. 23             AChar           iCode          PDB insertion code.
          25 - 28        LString         database       Sequence database name.
          30 - 38        LString         dbIdCode       Sequence database
     accession number. 40 - 42        Residue name    dbRes          Sequence
     database residue name. 44 - 48        Integer         dbSeq Sequence
     database sequence number. 50 - 70        LString         conflict Conflict
     comment.*/
  char seqadv[10];
  char idCode[10];
  char resName[10];
  char chainID;
  int seqNum;
  char iCode;
  char database[10];
  char dbIdCode[10];
  char dbRes[10];
  int dbSeq;
  char conflict[22];
  getString(line, seqadv, 0, 6);
  getString(line, idCode, 7, 4);
  getString(line, resName, 12, 3);
  getChar(line, &chainID, 16);
  getInt(line, &seqNum, 18, 4);
  getChar(line, &iCode, 22);
  getString(line, database, 24, 4);
  getString(line, dbIdCode, 29, 9);
  getString(line, dbRes, 39, 3);
  getInt(line, &dbSeq, 43, 5);
  getString(line, conflict, 49, 21);
}

// Craig: move to util
// Case insensitive string compare.  returns 0 if the two input strings are
// case-insensitive identical (in length and content) and returns non-zero
// otherwise.  Returns non-zero if either string is NULL.
bool PDBParser::strcmpCaseInsensitive(const char *str1, const char *str2) {
  if (str1 == NULL || str2 == NULL)
    return false; // either null
  int len1 = strlen(str1);
  int len2 = strlen(str2);
  if (len1 != len2)
    return false; // unequal length
  for (int i = 0; i < len1; i++) {
    if (tolower(str1[i]) != tolower(str2[i]))
      return false; // mismatch
  }
  return true; // all chars matched
}

// When a residue is created, it needs to be initialized with the type of
// residue it is (residueID) and the order in which it appeared in the file
// (resNum), as well as a standard ResidueHeader.
void parserPDBtoGOA::initResidue(GroupOfAtoms *residue,
                                 ResidueHeader *resHeader, int residueID,
                                 int resNum) {
  residue->type = RESIDUE_TYPE;
  residue->resHeader = resHeader;
  residue->resID = residueID;
  residue->resPos = resNum;
  // Set the color based on the residue type
  residue->setResidueColor(residueID);
}

void parserPDBtoGOA::parseSEQRES(GroupOfAtoms *molecule, char *line) {
  static char curChainId;
  static int curLength;
  /*1 -  6       Record name     "SEQRES"
  9 - 10         Integer         serNum        Serial number of the SEQRES
  record for the current chain. Starts at 1 and increments by one each line.
  Reset to 1 for each chain. 12             Character       chainID       Chain
  identifier.  This may be any single legal character, including a blank which
  is used if there is only one chain. 14 - 17        Integer         numRes
  Number of residues in the chain. This value is repeated on every record. 20 -
  22        Residue name    resName       Residue name. 24 - 26        Residue
  name    resName       Residue name. 28 - 30        Residue name    resName
  Residue name. 32 - 34        Residue name    resName       Residue name. 36 -
  38        Residue name    resName       Residue name. 40 - 42        Residue
  name    resName       Residue name. 44 - 46        Residue name    resName
  Residue name. 48 - 50        Residue name    resName       Residue name. 52 -
  54        Residue name    resName       Residue name. 56 - 58        Residue
  name    resName       Residue name. 60 - 62        Residue name    resName
  Residue name. 64 - 66        Residue name    resName       Residue name. 68 -
  70        Residue name    resName       Residue name.*/
  char seqres[10];
  int serNum;
  char chainID;
  int numRes;
  getString(line, seqres, 0, 6);
  getInt(line, &serNum, 8, 2);
  getChar(line, &chainID, 11);
  getInt(line, &numRes, 13, 4);
  if (!residueArray) {
    curChainId = chainID;
    curLength = numRes;
    residueArray = new int[numRes];
    chainTable[chainID] = 0;
  }
  // Grow the arrays when we reach a new chain
  if (curChainId != chainID) {
    int *newResidueArray = new int[numRes + curLength];
    chainTable[chainID] = curLength;
    for (int i = 0; i < curLength; i++) {
      newResidueArray[i] = residueArray[i];
    }
    delete[] residueArray;
    residueArray = 0; // delete the old arrays
    curLength += numRes;
    residueArray = newResidueArray;
    curChainId = chainID;
  }
  int pos = 19;
  for (int i = 0; i < 13; i++) {
    char resName[10];
    getString(line, resName, pos, 3);
    if ((strlen(resName) != 0) && (isalpha(resName[2]))) {
      int residue = Residues::mapResidue(resName);
      residueArray[residueAssigned] = residue;
      residueAssigned++;
    }
    pos += 4;
  }
}

void parserPDBtoGOA::parseHELIX(GroupOfAtoms *molecule, char *line) {
  /*      1 -  6        Record name     "HELIX "
           8 - 10        Integer         serNum          Serial number of the
     helix.   This starts at 1 and increases incrementally. 12 - 14 LString(3)
     helixID         Helix identifier. In addition to a serial number, each
     helix is given a char helix identifier. 16 - 18        Residue name
     initResName     Name of the initial residue. 20             Character
     initChainID     Chain identifier for the chain containing this helix. 22 -
     25        Integer         initSeqNum      Sequence number of the initial
     residue. 26             AChar           initICode       Insertion code of
     the initial residue. 28 - 30        Residue name    endResName      Name of
     the terminal residue of the helix. 32             Character endChainID
     Chain identifier for the chain containing this helix. 34 - 37 Integer
     endSeqNum       Sequence number of the terminal residue. 38 AChar endICode
     Insertion code of the terminal residue. 39 - 40        Integer helixClass
     Helix class (see below). 41 - 70        String          comment Comment
     about this helix. 72 - 76        Integer         length          Length of
     this helix.*/
  char helix[10];
  int serNum;
  char helixID[10];
  char initResName[10];
  char initChainID;
  int initSeqNum;
  char initICode;
  char endResName[10];
  char endChainID;
  int endSeqNum;
  char endICode;
  int helixClass;
  char comment[30];
  int length;
  getString(line, helix, 0, 6);
  getInt(line, &serNum, 7, 3);
  getString(line, helixID, 11, 3);
  getString(line, initResName, 15, 3);
  getChar(line, &initChainID, 19);
  getInt(line, &initSeqNum, 21, 4);
  getChar(line, &initICode, 25);
  getString(line, endResName, 27, 3);
  getChar(line, &endChainID, 31);
  getInt(line, &endSeqNum, 33, 4);
  getChar(line, &endICode, 37);
  getInt(line, &helixClass, 38, 2);
  getString(line, comment, 40, 30);
  getInt(line, &length, 71, 5);
  if (initChainID != endChainID) {
    // debug: need to give a warning somehow
    initChainID = initChainID;
  }
  GroupOfAtoms *helixGroup = new GroupOfAtoms();
  helixGroup->m_SecondaryStructureData = new SecondaryStructureData();
  helixGroup->type = SS_TYPE;
  helixGroup->setSSColor(SS_HELIX);
  helixGroup->resID = SS_HELIX;
  int initResID = Residues::mapResidue(initResName);
  int endResID = Residues::mapResidue(endResName);
  // Record the initial and the ending identifiers
  ResidueHeader *rh;
  rh = new ResidueHeader(initResID, initChainID, initSeqNum, initICode);
  helixGroup->startRes = rh;
  rh = new ResidueHeader(endResID, endChainID, endSeqNum, endICode);
  helixGroup->endRes = rh;
  helixGroups.add(helixGroup);
}

void parserPDBtoGOA::parseSHEET(GroupOfAtoms *molecule, char *line) {
  /*      1 -  6        Record name     "SHEET "
           8 - 10        Integer         strand          Strand number which
     starts at 1 for each strand within a sheet and increases by one. 12 - 14
     LString(3)      sheetID         Sheet identifier. 15 - 16        Integer
     numStrands      Number of strands in sheet. 18 - 20        Residue name
     initResName     Residue name of initial residue. 22             Character
     initChainID     Chain identifier of initial residue in strand. 23 - 26
     Integer         initSeqNum      Sequence number of initial residue in
     strand. 27             AChar           initICode       Insertion code of
     initial residue in strand. 29 - 31        Residue name    endResName
     Residue name of terminal residue. 33             Character       endChainID
     Chain identifier of terminal residue. 34 - 37        Integer endSeqNum
     Sequence number of terminal residue. 38             AChar endICode
     Insertion code of terminal residue. 39 - 40        Integer         sense
     Sense of strand with respect to previous strand in the sheet. 0 if first
     strand, 1 if parallel, -1 if anti-parallel. 42 - 45        Atom curAtom
     Registration. Atom name in current strand. 46 - 48        Residue name
     curResName      Registration. Residue name in current strand. 50 Character
     curChainId      Registration. Chain identifier in current strand. 51 - 54
     Integer         curResSeq       Registration. Residue sequence number in
     current strand. 55             AChar           curICode Registration.
     Insertion code in current strand. 57 - 60        Atom            prevAtom
     Registration. Atom name in previous strand. 61 - 63        Residue name
     prevResName     Registration. Residue name in previous strand. 65 Character
     prevChainId     Registration. Chain identifier in previous strand. 66 - 69
     Integer         prevResSeq      Registration. Residue sequence number in
     previous strand. 70             AChar           prevICode Registration.
     Insertion code in previous strand.*/
  char sheet[10];
  int strand;
  char sheetID[10];
  int numStrands;
  char initResName[10];
  char initChainID;
  int initSeqNum;
  char initICode;
  char endResName[10];
  char endChainID;
  int endSeqNum;
  char endICode;
  int sense;
  char curAtom[10];
  char curResName[10];
  char curChainId;
  int curResSeq;
  char curICode;
  char prevAtom[10];
  char prevResName[10];
  char prevChainId;
  int prevResSeq;
  char prevICode;
  getString(line, sheet, 0, 6);
  getInt(line, &strand, 7, 3);
  getString(line, sheetID, 11, 3);
  getInt(line, &numStrands, 14, 2);
  getString(line, initResName, 17, 3);
  getChar(line, &initChainID, 21);
  getInt(line, &initSeqNum, 22, 4);
  getChar(line, &initICode, 26);
  getString(line, endResName, 28, 3);
  getChar(line, &endChainID, 32);
  getInt(line, &endSeqNum, 33, 5);
  getChar(line, &endICode, 37);
  getInt(line, &sense, 38, 2);
  getString(line, curAtom, 41, 4);
  getString(line, curResName, 45, 3);
  getChar(line, &curChainId, 49);
  getInt(line, &curResSeq, 50, 4);
  getChar(line, &curICode, 54);
  getString(line, prevAtom, 56, 4);
  getString(line, prevResName, 60, 3);
  getChar(line, &prevChainId, 64);
  getInt(line, &prevResSeq, 65, 4);
  getChar(line, &prevICode, 69);
  if (initChainID != endChainID) {
    // debug: need to give a warning somehow
    initChainID = initChainID;
  }
  // we need to create a new groupOfAtoms for a new sheet only if we get a
  // strand with id == 1 else, we add the new strand to the existing sheet If we
  // get a bad case when the strand id is not 1 and a. we dont have any previous
  // GOA,
  //             we curse the pdb folk and create and add a new sheet with this
  //             strand as the first
  // b. we get a prev sheet, but with wrong number of strand,
  //             we curse the pdb folk and add this strand to prev sheet
  // Update the resSequenceNumberPairs array when necessary
  int numberOfSheets = sheetGroups.getNumberOfObjects();
  if ((strand == 1) || (numberOfSheets == 0)) {
    GroupOfAtoms *sheetGroup = new GroupOfAtoms();
    sheetGroup->m_SecondaryStructureData = new SecondaryStructureData();
    sheetGroup->type = SS_TYPE;
    sheetGroup->setSSColor(SS_SHEET);
    sheetGroup->resID = SS_SHEET;
    sheetGroup->m_SecondaryStructureData->m_ResSequenceNumberPairList.push_back(
        initSeqNum);
    sheetGroup->m_SecondaryStructureData->m_ResSequenceNumberPairList.push_back(
        endSeqNum);
    sheetGroups.add(sheetGroup);
    int initResID = Residues::mapResidue(initResName);
    int endResID = Residues::mapResidue(endResName);
    // Record the initial and the ending identifiers
    ResidueHeader *rh;
    rh = new ResidueHeader(initResID, initChainID, initSeqNum, initICode);
    sheetGroup->startRes = rh;
    rh = new ResidueHeader(endResID, endChainID, endSeqNum, endICode);
    sheetGroup->endRes = rh;
  } else {
    GroupOfAtoms *lastSheet = sheetGroups.get(numberOfSheets - 1);
    if (lastSheet->m_SecondaryStructureData == 0) // what to do ? major error!!
    {
      lastSheet->m_SecondaryStructureData =
          new SecondaryStructureData(); // :( bad workaround
    }
    int numberOfPreviousStrands = lastSheet->m_SecondaryStructureData
                                      ->m_ResSequenceNumberPairList.size() /
                                  2;
    if (numberOfPreviousStrands == strand - 1) {
      // ok, add new strand and update its array
    } else {
      // curse pdb
    }
    lastSheet->m_SecondaryStructureData->m_ResSequenceNumberPairList.push_back(
        initSeqNum);
    lastSheet->m_SecondaryStructureData->m_ResSequenceNumberPairList.push_back(
        endSeqNum);
    int endResID = Residues::mapResidue(endResName);
    // Record the the ending identifiers only ! we need to delete prev and
    // update
    ResidueHeader *rh;
    rh = new ResidueHeader(endResID, endChainID, endSeqNum, endICode);
    delete lastSheet->endRes;
    lastSheet->endRes = rh;
  }
}

void parserPDBtoGOA::parseTURN(GroupOfAtoms *molecule, char *line) {
  /*      1 -  6        Record name    "TURN  "
           8 - 10        Integer        seq            Turn number; starts with
     1 and increments by one. 12 - 14        LString(3)     turnId         Turn
     identifier 16 - 18        Residue name   initResName    Residue name of
     initial residue in turn. 20             Character      initChainId    Chain
     identifier for the chain containing this turn. 21 - 24        Integer
     initSeqNum     Sequence number of initial residue in turn. 25 AChar
     initICode      Insertion code of initial residue in turn. 27 - 29 Residue
     name   endResName     Residue name of terminal residue of turn. 31
     Character      endChainId     Chain identifier for the chain containing
     this turn. 32 - 35        Integer        endSeqNum      Sequence number of
     terminal residue of turn. 36             AChar          endICode Insertion
     code of terminal residue of turn. 41 - 70        String         comment
     Associated comment.*/
  char turn[10];
  int seq;
  char turnId[10];
  char initResName[10];
  char initChainId;
  int initSeqNum;
  char initICode;
  char endResName[10];
  char endChainId;
  int endSeqNum;
  char endICode;
  char comment[35];
  getString(line, turn, 0, 6);
  getInt(line, &seq, 7, 3);
  getString(line, turnId, 11, 3);
  getString(line, initResName, 15, 3);
  getChar(line, &initChainId, 19);
  getInt(line, &initSeqNum, 20, 4);
  getChar(line, &initICode, 24);
  getString(line, endResName, 26, 3);
  getChar(line, &endChainId, 30);
  getInt(line, &endSeqNum, 31, 4);
  getChar(line, &endICode, 35);
  getString(line, comment, 40, 30);
  if (initChainId != endChainId) {
    // debug: need to give a warning somehow
    initChainId = initChainId;
  }
  GroupOfAtoms *turnGroup = new GroupOfAtoms();
  turnGroup->m_SecondaryStructureData = new SecondaryStructureData();
  turnGroup->type = SS_TYPE;
  turnGroup->setSSColor(SS_TURN);
  turnGroup->resID = SS_TURN;
  int initResID = Residues::mapResidue(initResName);
  int endResID = Residues::mapResidue(endResName);
  // Record the initial and the ending identifiers
  ResidueHeader *rh;
  rh = new ResidueHeader(initResID, initChainId, initSeqNum, initICode);
  turnGroup->startRes = rh;
  rh = new ResidueHeader(endResID, endChainId, endSeqNum, endICode);
  turnGroup->endRes = rh;
  turnGroups.add(turnGroup);
}

// SKVINAY: I added another parameter, fileType to indicate whether this is a
// PDB file or a PQR file. There is another simple routine called parsepqr which
// does not build a hierarchy. I think this function is better as it builds a
// dummy hierarchy, allowing pqr files to be dealt similar with pdb ones. We
// need this parameter because the atom records have different formats in the
// two files. why ?
void parserPDBtoGOA::parseATOM(GroupOfAtoms *molecule, char *line,
                               bool boundaryAtom,
                               CollectionData *collectionData,
                               bool deleteWater) {

  /*      1 -  6        Record name     "ATOM  "
           7 - 11        Integer         serial        Atom serial number.
          13 - 16        Atom            name          Atom name.
          17             Character       altLoc        Alternate location
     indicator. 18 - 20        Residue name    resName       Residue name. 22
     Character       chainID       Chain identifier. 23 - 26        Integer
     resSeq        Residue sequence number. 27             AChar           iCode
     Code for insertion of residues. 31 - 38        Real(8.3)       x Orthogonal
     coordinates for X in Angstroms. 39 - 46        Real(8.3)       y Orthogonal
     coordinates for Y in Angstroms. 47 - 54        Real(8.3)       z Orthogonal
     coordinates for Z in Angstroms. 55 - 60        Real(6.2)       occupancy
     Occupancy. 61 - 66        Real(6.2)       tempFactor    Temperature factor.
          73 - 76        LString(4)      segID         Segment identifier,
     left-justified. 77 - 78        LString(2)      element       Element
     symbol, right-justified.
          // This 'element' field is not always there! It should never be used.
     Use 'Atom name' instead. 79 - 80        LString(2)      charge Charge on
     the atom.*/
  char atomname[10] = "ATOM";
  int serial = 0;
  char name[10] = "XXXX";
  char altLoc = 0;
  char resName[10] = "XXX";
  char chainID = 0;
  int resSeq = 0;
  char iCode = '0';
  double x = 0;
  double y = 0;
  double z = 0;
  double occupancy = 1;
  double tempFactor = 1;
  char segID[10] = "X";
  char element[10] = "X";
  char charge[10] = "X";
  bool chainChanged = false;
  bool ssChanged = false;
  bool residueChanged = false;
  getString(line, atomname, 0, 6);
  getInt(line, &serial, 6, 5);
  getString(line, name, 12, 4);
  getChar(line, &altLoc, 16);
  getString(line, resName, 17, 3);
  getChar(line, &chainID, 21);
  getInt(line, &resSeq, 22, 4);
  getChar(line, &iCode, 26);
  getDouble(line, &x, 30, 8);
  getDouble(line, &y, 38, 8);
  getDouble(line, &z, 46, 8);
  getDouble(line, &occupancy, 54, 6);
  getDouble(line, &tempFactor, 60, 6);
  getString(line, segID, 72, 4);
  getString(line, element, 76, 2);
  getString(line, charge, 78, 2);
  if (deleteWater) {
    if (strcmp(resName, "HOH") == 0) {
      return; // dont want water
    }
  }
  xmin = minimum(x, xmin);
  ymin = minimum(y, ymin);
  zmin = minimum(z, zmin);
  xmax = maximum(x, xmin);
  ymax = maximum(y, ymin);
  zmax = maximum(z, zmin);
  int currResID = Residues::mapResidue(
      resName); // The int identifier for the current residue
  // Update the chain
  if (chainID != chainIDLast) {
    // add number of residues read for previous chain, if there was any!
    if (m_NumberOfResiduesInCurrentChain) {
      collectionData->m_ChainIds.push_back(m_NumberOfResiduesInCurrentChain);
      m_NumberOfResiduesInCurrentChain = 0;
    }
    collectionData->m_ChainIds.push_back((int)(chainID));
    currChain = new GroupOfAtoms();
    currChain->type = CHAIN_TYPE;
    static int chainNum = 0;
    // Random color
    srand(3.14159 / 2.0 - chainID + 65 + 19.2);
    double blue = rand() / (RAND_MAX + 1.0f);
    double red = rand() / (RAND_MAX + 1.0f);
    double green = rand() / (RAND_MAX + 1.0f);
    currChain->setGOAColor(red, green, blue);
    chainNum++;
    if (chainNum == 4) {
      chainNum = 0;
    }
    currChain->resID = (int)chainID;
    if (currResID >= RES_TYPE_FIRST_NUCLEOTIDE) {
      currChain->resPos = CHAIN_NUCLEIC;
    } else {
      currChain->resPos = CHAIN_PROTEIN;
    }
    // Add to the set of all chains
    chainGroups.add(currChain);
    chainIDLast = chainID; // update
    chainChanged = true;
    // pcd debug: on 1J5E.pdb, we need to unhook the current secondary structure
    // manually when the chain changes, need to look into why the normal
    // unhooking code (unhook if the terminating residue of an SS is seen) is
    // not working, this noticably occurs for the last chains T and V
    currSS = 0;
  }
  // End Chain update section
  // Begin chain / residue / ss update section
  // Check if this atom is the start of a new residue.
  // If so, make a new GOA for the new residue and add it to the appropriate
  // array. Also, error check.  Since it is possible for some PDB's do not
  // contain any preliminary SEQRES information before the ATOM fields, a
  // special exception is made in these error checks for when no prior SEQRES
  // information has been seen.
  if (residueHasChanged(chainID, resSeq, iCode) == true) {
    // The number of residues from the ATOM fields (the check count) cannot
    // exceed the number of residues declared from the SEQRES fields.  The '>'
    // comparator is used here because on the final residue, all atoms after the
    // first one will have a residueCheckCount equal to the value of
    // residueAssigned. Update the variables that represent the current residue
    residueSeqLast = resSeq;
    iCodeLast = iCode;
    // Create a new residue and add it to the appropriate array.
    currResidue = new GroupOfAtoms();
    ResidueHeader *currHeader =
        new ResidueHeader(currResID, chainID, resSeq, iCode);
    //            initResidue(currResidue, currHeader, currResID,
    //            residueCheckCount);
    // Pass in the serial number of the first atom in the residue, this way,
    // residues in data structures can be quickly linked back to the part of the
    // ATOMS that created them
    initResidue(currResidue, currHeader, currResID, serial);
    residueChanged = true;
    // add to collection data's residue information.
    collectionData->m_ResidueIds.push_back(
        (int)(Residues::lookupShortResID(resName)));
    collectionData->m_ResidueIds.push_back(resSeq);
    collectionData->m_ResidueIds.push_back((int)iCode);
    m_NumberOfResiduesInCurrentChain++;
    // Begin Residue / SS update section
    // Assumption: a residue can belong to at most one SS.  Code will need to be
    // changed if this is not true. Assumption: helices and turns will appear in
    // order ascending residue sequence number in a file, and sheet entries can
    // appear in any order. Unlike chains and residues, all SS's are finalized
    // and moved to the appropriate location when the last element is parsed
    // because the SS information gathered earlier indicates the last residue in
    // the SS.  For chains and residues, we only know the chain or residue is
    // done when we transition to a new chain or residue or if we reach the end
    // of the file. Check if we have a real secondary structure
    //  - currSS == NULL on startup and when the last residue before the current
    //  one was at the
    //              end of an SS
    //  - currSS->resID == SS_NULL when we are in the middle of processing a
    //  stream
    //                     of residues that do not belong to a real secondary
    //                     structure
    //
    // Therefore, when both of these conditions are not true, we in the middle
    // of parsing an actual secondary strucutre, and we can add the new residue
    // directly to the SS. Ensure that if we are in a sheet, we allow other SS's
    // to be within its strands!
    if (currSS != 0 && currSS->resID != SS_NULL && currSS->resID != SS_SHEET) {
      // If there is a current SS, place the residue there and check if it is
      // the terminating residue to the SS
      currSS->m_SubGroups.push_back(currResidue);
      if (currSS->endRes->compare(currHeader) == true) {
        // Remove focus from the current SS, allowing a new SS to take its place
        // next trime
        currSS = 0;
      }
    } else {
      // Check if the start residue of any of the next available SS's match the
      // current residue If so, switch to that SS.
      GroupOfAtoms *tempSS;
      // Check the next helix
      {
        if (helixGroupNext < (int)helixGroups.getNumberOfObjects()) {
          tempSS = helixGroups.get(helixGroupNext);

          // If the current residue matches the starting residue, switch
          if (tempSS->startRes->compare(currHeader) == true) {
            currSS = tempSS;
            helixGroupNext++;
            // The newly found SS needs to be added to the chain
            currChain->m_SubGroups.push_back(currSS);
            ssChanged = true;
          }
        }
      }
      // Check the next turn
      {
        if (turnGroupNext < (int)turnGroups.getNumberOfObjects()) {
          tempSS = turnGroups.get(turnGroupNext);
          // If the current residue matches the starting residue, switch
          if (tempSS->startRes->compare(currHeader) == true) {
            currSS = tempSS;
            turnGroupNext++;
            // The newly found SS needs to be added to the chain
            currChain->m_SubGroups.push_back(currSS);
            ssChanged = true;
          }
        }
      }
      // Check all of the sheets, since sheets are unordered in the input file
      // in each sheet, we need to go through and see if the residue is part of
      // some strand A Sheet is not a continuous list of residues, but a
      // collection of strands.
      bool belongsToSheet = false;
      {
        // check to make sure this residue is not part of a helix!
        // for some idiotic reason, thats quite common!!
        bool screwedUpDuplicateEntry = false;
        {
          int n = (int)helixGroups.getNumberOfObjects();
          for (int i = 0; i < n; i++) {
            GroupOfAtoms *tempHelix = helixGroups.get(i);
            if ((tempHelix->startRes->chainID == currHeader->chainID) &&
                (tempHelix->startRes->seqNum <= currHeader->seqNum) &&
                (tempHelix->endRes->seqNum >= currHeader->seqNum)) {
              screwedUpDuplicateEntry = true;
            }
          }
        }
        if (!screwedUpDuplicateEntry) {
          int n = (int)sheetGroups.getNumberOfObjects();
          for (sheetGroupNext = 0; sheetGroupNext < n; sheetGroupNext++) {
            tempSS = sheetGroups.get(sheetGroupNext);
            // If the current residue matches the starting residue, switch
            // if (tempSS->startRes->compare(currHeader) == true)
            if (tempSS->m_SecondaryStructureData != 0) {
              // now check to see if the residue is within any strand of the
              // sheet only then add to sheet! the SS data should be allocated,
              // otherwise we cry
              int numberOfStrands = tempSS->m_SecondaryStructureData
                                        ->m_ResSequenceNumberPairList.size() /
                                    2;
              if (numberOfStrands > 0) {
                for (int i = 0; i < numberOfStrands; i++) {
                  int startResSeqId =
                      tempSS->m_SecondaryStructureData
                          ->m_ResSequenceNumberPairList[i * 2 + 0];
                  int endResSeqId =
                      tempSS->m_SecondaryStructureData
                          ->m_ResSequenceNumberPairList[i * 2 + 1];
                  int ssChainId = tempSS->startRes->chainID;
                  if ((resSeq >= startResSeqId) && (resSeq <= endResSeqId) &&
                      (ssChainId == chainID)) {
                    belongsToSheet = true;
                  }
                }
              }
              if (belongsToSheet) {
                currSS = tempSS;
                // The newly found SS needs to be added to the chain, but only
                // once ! go through chain and see if its already in it, if not,
                // add it
                int c = currChain->m_SubGroups.size();
                bool sheetInChain = false;
                for (int i = 0; i < c; i++) {
                  if (currChain->m_SubGroups.at(i) == currSS) {
                    sheetInChain = true;
                  }
                }
                if (!sheetInChain) {
                  currChain->m_SubGroups.push_back(currSS);
                }
                ssChanged = true;
                break;
              }
            }
          }
        }
      }
      // WARNING: at this point, we have already seen all SS's and determined
      // that its not part of anything or added it to something ! This implies
      // that currSS could point to a prev found Sheet, but the residue may not
      // belong to it. Hence, if either the currSS is null, or if it was sheet,
      // but was not found, create a NULL and add to it! Check if we have found
      // an SS that starts with the current residue.  If we have not, create a
      // null secondary structure that will represent the absense of a secondary
      // structure.
      if ((currSS == 0) || ((currSS->resID == SS_SHEET) && (!belongsToSheet))) {
        // Add the residue to a null secondary structure, which means all
        // component residues do not belong to a real secondary structure (i.e.
        // helix, sheet, or turn).
        currSS = new GroupOfAtoms();
        currSS->m_SecondaryStructureData = new SecondaryStructureData();
        currSS->type = SS_TYPE;
        currSS->setSSColor(SS_NULL);
        currSS->resID = SS_NULL;
        // The null SS needs to be added to the chain so that the hierarchy is
        // connected.
        currChain->m_SubGroups.push_back(currSS);
      }
      // Add the new residue to the current secondary structure, thus linking
      // up the hierarchy.
      currSS->m_SubGroups.push_back(currResidue);
      // If we have a non-NULL secondary structure, we need to check if the
      // current residue is the last residue in the secondary structure.
      if (currSS->resID != SS_NULL &&
          currSS->endRes->compare(currHeader) == true) {
        // Remove focus from the current SS, allowing a new SS to take its place
        // next time
        currSS = 0;
      }
    }
    // End Residue / SS update section
    residueCheckCount++; // Increase the number of residues seen while error
                         // checking
  }
  // End chain / residue / ss update section
  // get the value for m_Element of atoms. This should be a table lookup from
  // 'name'
  int elementIndex = 0;
  // if( fileType == 0 ) // pdb
  //      elementIndex = mapElement( name );
  // else
  //      elementIndex = mapPQRElement( name );
  elementIndex = mapElementAndResidue(name, resName);
  // if( elementIndex >= 0 )
  currResidue->addAtom(serial, elementIndex, resSeq, name, currResID, chainID,
                       x, y, z, occupancy, 0, boundaryAtom);
  // Error check the flags residueChanged, ssChanged, and chainChanged
  // The following cases are produce warnings:
  //
  // 1.  The chain changed but the ss did not change
  // 2.  The chain changed but the residue did not change
  // 3.  The ss changed but the residue did not change
  // debug: add checks that compare these flags to each other
}

void parserPDBtoGOA::parseHBOND(GroupOfAtoms *molecule, const char *fname) {
  FILE *fp = fopen(fname, "r");
  if (!fp) {
    printf("Hydrogen-bond file open fail: %s\n", fname);
    return;
  }
  char line[1025];
  // read total number of atoms
  int natomsR = 0, natomsL = 0, natoms = 0;
  fgets(line, 1024, fp);
  sscanf(line, "%d %d", &natomsR, &natomsL);
  fprintf(stderr, "receptor: %d, ligand: %d\n", natomsR, natomsL);
  natoms = natomsR + natomsL;
  molecule->setNumberOfReceptorAtoms(natomsR);
  molecule->setNumberOfLigandAtoms(natomsL);

  Atom **atomlist = new Atom *[natoms + 1];
  flattenAtom(molecule, atomlist, natoms);
  fprintf(stderr, "atom flatten\n");

  std::vector<HBOND> *_hbond = molecule->getHBond();

  while (NULL != fgets(line, 1024, fp)) {
    int sid, eid;
    float e;
    sscanf(line, "%d %d %f", &sid, &eid, &e);
    HBOND hb;
    if (sid > natoms || eid > natoms) {
      fprintf(stderr, "hbond id: %d %d\n", sid, eid);
      continue;
    }
    hb.d = atomlist[sid];
    hb.a = atomlist[eid];
    hb.energy = e;

    if (((sid <= natomsR) && (eid <= natomsR)) ||
        ((sid >= natomsR) && (eid >= natomsR)))
      hb.internal = true;
    else
      hb.internal = false;

    hb.d->m_donorOrAcceptor = true;
    hb.a->m_donorOrAcceptor = true;

    int isInternal = hb.internal ? 2 : 1;

    hb.d->m_internalHbond = hb.d->m_internalHbond | isInternal;
    hb.a->m_internalHbond = hb.a->m_internalHbond | isInternal;

    _hbond->push_back(hb);
  }
  fprintf(stderr, "%d hbond have been parsed\n", _hbond->size());
  delete[] atomlist;
}

void parserPDBtoGOA::flattenAtom(PDBParser::GroupOfAtoms *molecule,
                                 Atom **atomlist, int maxIdx) {
  // CHA: atomlist must have a length of total-number-of-atoms + 1
  if (!molecule) {
    return;
  }

  // If the proper level of detail is set, render the group of atoms in
  // multiresolution form but not the component subgroups
  if (molecule->type == PDBParser::RESIDUE_TYPE) {
    for (int i = 0; i < molecule->m_Atoms.size(); i++) {
      PDBParser::Atom *atom = molecule->m_Atoms.at(i);
      if (atom->getAtomNum() > maxIdx)
        fprintf(stderr, "atom id: %d\n", atom->getAtomNum());
      else
        atomlist[atom->getAtomNum()] = atom;
    }
  }

  // If control has reached this point, iterate through the subgroup array and
  // draw all subgroups
  if (molecule->type == PDBParser::COLLECTION_TYPE) {
    PDBParser::GroupOfAtoms *currSubGroup =
        molecule->m_SubGroups.at(0); // render only the first model

    if (currSubGroup) {
      flattenAtom(currSubGroup, atomlist, maxIdx);
    }
  } else {
    int n = (int)molecule->m_SubGroups.size();
    for (int i = 0; i < n; i++) {
      PDBParser::GroupOfAtoms *currSubGroup = molecule->m_SubGroups.at(i);

      if (currSubGroup) {
        flattenAtom(currSubGroup, atomlist, maxIdx);
      }
    }
  }
}

// Transfers the subgroups to the final molecule, does some postprocessing, and
// then readies the parser for parsing a new molecule.  Call this routine once
// the entire PDB file has been parsed.
void parserPDBtoGOA::finalizeParsing(GroupOfAtoms *molecule) {
  int n = (int)chainGroups.getNumberOfObjects();

  for (int i = 0; i < n; i++) {
    GroupOfAtoms *currChain = chainGroups.get(i);

    molecule->m_SubGroups.push_back(currChain);
    // Calculate either protein or chain torsion angles
    if (currChain->resPos == CHAIN_PROTEIN) {
      calculateProteinChainTorsionAngles(currChain);
    } else if (currChain->resPos == CHAIN_NUCLEIC) {
      calculateNucleicAcidChainTorsionAngles(currChain);
    }
  }
  chainGroups.nonDestructiveClear();
  helixGroups.nonDestructiveClear();
  sheetGroups.nonDestructiveClear();
  turnGroups.nonDestructiveClear();
  helixGroupNext = 0;
  sheetGroupNext = 0;
  turnGroupNext = 0;
  residueCheckCount = 0;
  chainIDLast = '\0';
  residueSeqLast = -999999; // Sequence numbers can be negative, but only 4
                            // digits, so -999999 means uninitialized
  iCodeLast = '\0';
  currResidue = 0;
  currSS = 0;
  currChain = 0;

  // showGOAhierarchy(molecule, 0);
}

void parserPDBtoGOA::showGOAhierarchy(GroupOfAtoms *molecule, int level) {
  printf("Level: %d, Type: ", level);
  switch (molecule->type) {
  case GENERAL_TYPE:
    printf("GENERAL_TYPE");
    break;
  case ATOM_TYPE:
    printf("ATOM_TYPE");
    break;
  case RESIDUE_TYPE:
    printf("RESIDUE_TYPE");
    break;
  case SS_TYPE:
    printf("SS_TYPE");
    break;
  case PROTEIN_TYPE:
    printf("PROTEIN_TYPE");
    break;
  case CHAIN_TYPE:
    printf("CHAIN_TYPE");
    break;
  case COLLECTION_TYPE:
    printf("COLLECTION_TYPE\n");
    break;
  default:
    printf("NONE");
  }
  printf(", resID: ");
  switch (molecule->resID) {
  case SS_HELIX:
    printf("SS_HELIX\n");
    break;
  case SS_SHEET:
    printf("SS_SHEET\n");
    break;
  case SS_TURN:
    printf("SS_TURN\n");
    break;
  case SS_NULL:
    printf("SS_NULL\n");
    break;
  case TOTAL_SS_TYPE:
    printf("TOTAL_SS_TYPE\n");
    break;
  default:
    printf("NONE\n");
  }

  for (int i = 0; i < molecule->m_SubGroups.size(); i++)
    showGOAhierarchy(molecule->m_SubGroups[i], level + 1);
}

// Move to utility, deprecate
bool beginsWithCStyle(char *string, char *substring) {
  char *temp;
  if (temp = strstr(string, substring)) {
    int index = temp - string + 1;
    if ((temp != 0) && (index == 1)) {
      return true;
    }
  }
  return false;
}

// This is an internal local function. DONT call it.
// The collection data structure holds the header of the pdb and other global
// information about it. Functionality has been subsumed by parsePDBModels.
// Returns the first molecule in the PDB file.
GroupOfAtoms *parserPDBtoGOA::parsePDB(string filename,
                                       CollectionData *collectionData,
                                       bool deleteWater) {
  // Main file
  FILE *fp = fileRead(filename);

  // Boundary file
  string bdyFileName = filename + ".bdy";
  FILE *bdyfp = fopen(bdyFileName.c_str(), "r");
  bool boundaryPresent = bdyfp;

  // Hbond file name
  std::vector<std::string> hbond_files;

  // Current molecule
  currMolecule = new GroupOfAtoms();
  currMolecule->type = PROTEIN_TYPE;
  currMolecule->setMoleculeColor();
  // Parse the lines
  char line[129];
  memset(line, '\0',
         sizeof(line)); // Clear out the buffer in case of short lines
  while (fgets(line, 128, fp)) {
    char *str;
    int index;
    // Assumption: Records come in order as specified by PDB. No state
    // information needs to be stored. Exception:  Know the last chain and
    // residue from the prior line to know when a new chain or residue is
    // reached. Only the parsing routines for MODEL and ENDMDL will mutate
    // currMolecule.
    if ((str = strstr(line, "HEADER"))) {
      index = str - line + 1;
      if ((str != 0) && (index == 1)) {
        parseHEADER(line, collectionData);
      }
    } else if ((str = strstr(line, "TITLE"))) {
      index = str - line + 1;
      if ((str != 0) && (index == 1)) {
        parseTITLE(line, collectionData);
      }
    } else if ((str = strstr(line, "COMPND"))) {
      index = str - line + 1;
      if ((str != 0) && (index == 1)) {
        parseCOMPND(line, collectionData);
      }
    } else if ((str = strstr(line, "SOURCE"))) {
      index = str - line + 1;
      if ((str != 0) && (index == 1)) {
        parseSOURCE(line, collectionData);
      }
    } else if ((str = strstr(line, "KEYWDS"))) {
      index = str - line + 1;
      if ((str != 0) && (index == 1)) {
        parseKEYWDS(line, collectionData);
      }
    } else if ((str = strstr(line, "EXPDTA"))) {
      index = str - line + 1;
      if ((str != 0) && (index == 1)) {
        parseEXPDTA(line, collectionData);
      }
    } else if ((str = strstr(line, "AUTHOR"))) {
      index = str - line + 1;
      if ((str != 0) && (index == 1)) {
        parseAUTHOR(line, collectionData);
      }
    } else if ((str = strstr(line, "JRNL"))) {
      index = str - line + 1;
      if ((str != 0) && (index == 1)) {
        parseJRNL(line, collectionData);
      }
    } else if ((str = strstr(line, "REMARK"))) {
      index = str - line + 1;
      if ((str != 0) && (index == 1)) {
        parseREMARK(line, collectionData);
      }
    } else if ((str = strstr(line, "MODEL"))) {
      index = str - line + 1;
      if ((str != 0) && (index == 1)) {
        parseMODEL(currMolecule, line);
      }
    } else if ((str = strstr(line, "ENDMDL"))) {
      index = str - line + 1;
      if ((str != 0) && (index == 1)) {
        parseENDMDL(currMolecule, line);
      }
    } else if ((str = strstr(line, "SEQADV"))) {
      index = str - line + 1;
      if ((str != 0) && (index == 1)) {
        parseSEQADV(currMolecule, line);
      }
    } else if ((str = strstr(line, "SEQRES"))) {
      index = str - line + 1;
      if ((str != 0) && (index == 1)) {
        parseSEQRES(currMolecule, line);
      }
    } else if ((str = strstr(line, "HELIX"))) {
      index = str - line + 1;
      if ((str != 0) && (index == 1)) {
        parseHELIX(currMolecule, line);
      }
    } else if ((str = strstr(line, "SHEET"))) {
      index = str - line + 1;
      if ((str != 0) && (index == 1)) {
        parseSHEET(currMolecule, line);
      }
    } else if ((str = strstr(line, "TURN"))) {
      index = str - line + 1;
      if ((str != 0) && (index == 1)) {
        parseTURN(currMolecule, line);
      }
    } else if ((str = strstr(line, "ATOM")) || (str = strstr(line, "HETATM"))) {
      bool boundaryAtom = false;
      if (boundaryPresent) {
        double sasArea = 0.0;
        int temp = fscanf(bdyfp, "%lf\n", &sasArea);
        if (temp != 1)
          error("Error reading sasArea from " + string(filename) + " on line " +
                string(line));
        assert(temp != EOF);
        if (sasArea > sasThreshold) {
          boundaryAtom = true;
        }
      }
      index = str - line + 1;
      if ((str != 0) && (index == 1)) {
        parseATOM(currMolecule, line, boundaryAtom, collectionData,
                  deleteWater);
      }
    } else if ((str = strstr(line, "HBOND"))) {
      fprintf(stderr, "%s", line);
      char fname[1024];
      sscanf(line, "HBOND %s", fname);
      hbond_files.push_back(std::string(fname));
    }
    memset(line, '\0',
           sizeof(line)); // Clear out the buffer in case of short lines
  }
  /*
                       if(beginsWith(line, "HEADER"))     parseHEADER(line,
     collectionData); else if(beginsWith(line, "TITLE" ))     parseTITLE (line,
     collectionData); else if(beginsWith(line, "COMPND"))     parseCOMPND(line,
     collectionData); else if(beginsWith(line, "SOURCE"))     parseSOURCE(line,
     collectionData); else if(beginsWith(line, "KEYWDS"))     parseKEYWDS(line,
     collectionData); else if(beginsWith(line, "EXPDTA"))     parseEXPDTA(line,
     collectionData); else if(beginsWith(line, "AUTHOR"))     parseAUTHOR(line,
     collectionData); else if(beginsWith(line, "JRNL"  ))     parseJRNL  (line,
     collectionData); else if(beginsWith(line, "REMARK"))     parseREMARK(line,
     collectionData); else if(beginsWith(line, "MODEL" ))     parseMODEL
     (currMolecule, line); else if(beginsWith(line, "ENDMDL"))
     parseENDMDL(currMolecule, line); else if(beginsWith(line, "SEQADV"))
     parseSEQADV(currMolecule, line); else if(beginsWith(line, "SEQRES"))
     parseSEQRES(currMolecule, line); else if(beginsWith(line, "HELIX" ))
     parseHELIX (currMolecule, line); else if(beginsWith(line, "SHEET" ))
     parseSHEET (currMolecule, line); else if(beginsWith(line, "TURN"  ))
     parseTURN  (currMolecule, line); else if(beginsWith(line, "ATOM") ||
     beginsWith(line, "HETATM"))
  */
  // showGOAhierarchy( chainGroups.get(0), 0 );

  // Finalize the current molecule
  if (currMolecule) {
    finalizeParsing(currMolecule);
    molModels.add(currMolecule);
  }
  // Update all molecules seen with the global bounding box information.
  // Also, create the multiresolution representation for each molecule.
  // For speed, comment out the creation of multiresolution if it is not needed.
  int n = molModels.getNumberOfObjects();

  for (int i = 0; i < n; i++) {
    currMolecule = molModels.get(i);
    currMolecule->buildHierarchicalInformation(
        (xmin + xmax) / 2.0, (ymin + ymax) / 2.0, (zmin + zmax) / 2.0);
  }
  if (boundaryPresent) {
    fclose(bdyfp);
  }

  // CHA: parse HBOND file, if any. And construct the hbond structure
  // currently, we assume only one file name was appeared in the file.
  if (hbond_files.size() == 1)
    parseHBOND(molModels.get(0), hbond_files[0].c_str());

  fprintf(stderr, "n hbond: %d\n", molModels.get(0)->getHBond()->size());
  // Return the first molecule that was seen
  return molModels.get(0);
}

// Function returns an index in the table 'elementTable'
// given the 4 character long atom name and 3 character long residue name
// if the atom is not found, it returns -1
int parserPDBtoGOA::mapElementAndResidue(char *atomName, char *residueName) {
  if (!atomName || !residueName) {
    return -1;
  }
  int len = strlen(atomName);
  if (len < 4) {
    return -1;
  }
  len = strlen(residueName);
  if (len < 3) {
    return -1;
  }
  for (int i = 0; i < MAX_BIOCHEM_ELEMENTS; i++) {
    ElementInformation eInfo = elementTable[i];
    if (strcmp(trim(eInfo.atomName), trim(atomName)) == 0 &&
        strcmp(trim(eInfo.residueName), trim(residueName)) == 0) {
      return i;
    }
  }

  return -1;
}

// WARNING: check what pqr is all about and change this.
// this messes up with hydrogens anyway, cause it says unknown residue and maps
// it to XXX
int parserPDBtoGOA::mapPQRElementAndResidue(const char *const atomName,
                                            const char *const residueName) {
  if (!atomName || !residueName) {
    return -1;
  }
  char tableAtomName[20];
  char tableResidueName[20];
  int i;
  for (i = 0; i < MAX_BIOCHEM_ELEMENTS; i++) {
    ElementInformation eInfo = elementTable[i];
    sscanf(eInfo.atomName, "%s", tableAtomName);
    sscanf(eInfo.residueName, "%s", tableResidueName);
    if (!strcmp(tableAtomName, atomName) &&
        !strcmp(tableResidueName, residueName)) {
      return i;
    }
  }
  return -1;
}

void parserPDBtoGOA::getChainIds(GroupOfAtoms *molecule,
                                 std::vector<char> *chainIds) {
  if (!molecule || !chainIds) {
    return;
  }
  if (molecule->type == CHAIN_TYPE) {
    chainIds->push_back(molecule->resID); // this variable is overloaded
  }
  int i, m;
  m = (int)molecule->m_SubGroups.size();
  for (i = 0; i < m; i++) {
    GroupOfAtoms *currSubGroup = molecule->m_SubGroups.at(i);
    if (currSubGroup) {
      getChainIds(currSubGroup, chainIds);
    }
  }
}

// If this is a pdb1 file type, we add the 60 matrices.
// If any matrices are present already, we just return.
// We also use the chains only from the first model.
bool parserPDBtoGOA::addTransformationsIfNeeded(GroupOfAtoms *molecule) {
  if (!(m_FileType == PDB1_FILE)) {
    return true;
  }
  if (!molecule) {
    return false;
  }
  if (molecule->type != COLLECTION_TYPE) {
    return false;
  }
  // parse chains
  std::vector<char> chainIds;
  getChainIds(molecule, &chainIds);
  if (chainIds.size() < 1) {
    return false;
  }
  if (!molecule->m_CollectionData) {
    molecule->m_CollectionData = new CollectionData();
  }
  BiomoleculeTransformations *b = new BiomoleculeTransformations();
  b->addChains(chainIds);
  int i;
  for (i = 0; i < 60 * 3; i++) {
    b->rotation.push_back(virusTransformations[i * 3 + 0]);
    b->rotation.push_back(virusTransformations[i * 3 + 1]);
    b->rotation.push_back(virusTransformations[i * 3 + 2]);
    b->translation.push_back(0);
  }
  molecule->m_CollectionData->m_BiomoleculeTransformations.push_back(b);
  return true;
}

//  Sets and return the file type.
// The function returns false if it cannot determine the file type.
// The variable m_FileType is also set.
bool parserPDBtoGOA::getFileType(const char *const fileName) {
  m_FileType = UNKNOWN_FILE;
  int len = strlen(fileName);
  if (len < 5) {
    return false;
  }
  if ((fileName[len - 3] == 'p' || fileName[len - 3] == 'P') &&
      (fileName[len - 2] == 'd' || fileName[len - 2] == 'D') &&
      (fileName[len - 1] == 'b' || fileName[len - 1] == 'B')) {
    m_FileType = PDB_FILE;
  } else if ((fileName[len - 4] == 'p' || fileName[len - 4] == 'P') &&
             (fileName[len - 3] == 'd' || fileName[len - 3] == 'D') &&
             (fileName[len - 2] == 'b' || fileName[len - 2] == 'B') &&
             (fileName[len - 1] == '1')) {
    if (len == 5) {
      return false;
    }
    m_FileType = PDB1_FILE;
  } else {
    return false;
  }
  return true;
}

// This is the correct function to call to parse both a pdb and a pqr file.
// It returns a collection type.
// Returns a GOA with all of the models (aka molecules) in the subgroups array
// of the return value GOA. We should make sure that it works even for molecules
// with some symmetry information etc. If there are no models in the file (i.e.
// if the file only contains a single molecule), then the single molecule will
// be in the first slot in the returned GOA.
GroupOfAtoms *parserPDBtoGOA::parsePDBModels(const char *fileName,
                                             bool deleteWater) {
  if (!getFileType(fileName)) {
    return 0;
  }
  GroupOfAtoms *allModels = new GroupOfAtoms();
  allModels->type = COLLECTION_TYPE;
  allModels->m_CollectionData = new CollectionData();
  // Parse the file.  The models will be placed in the molModels class member.
  if (!parsePDB(fileName, allModels->m_CollectionData, deleteWater)) {
    delete allModels;
    return 0;
  }
  if (m_NumberOfResiduesInCurrentChain > 0) {
    // need to do this in the end, but the file formats do not have an end
    // marker all the time. So doing it here.
    allModels->m_CollectionData->m_ChainIds.push_back(
        m_NumberOfResiduesInCurrentChain);
    m_NumberOfResiduesInCurrentChain = 0;
  }

  // Transfer all of the models from parsing into the subgroups array of the
  // return GOA
  int n = molModels.getNumberOfObjects();
  for (int i = 0; i < n; i++) {
    allModels->m_SubGroups.push_back(molModels.get(i));
  }
  if (m_FileType == PDB1_FILE) {
    addTransformationsIfNeeded(allModels);
  }
  return allModels;
}
