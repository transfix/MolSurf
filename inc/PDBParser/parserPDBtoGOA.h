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
#if !defined(                                                                  \
    AFX_PARSERPDBTOGOA_H__06E4D08D_156E_46F5_9FA0_2B2FE9CCA790__INCLUDED_)
#define AFX_PARSERPDBTOGOA_H__06E4D08D_156E_46F5_9FA0_2B2FE9CCA790__INCLUDED_

#include <Utility/utility.h>
#include <PDBParser/ExpandableArray.h>

namespace PDBParser {
class Atom;
class GroupOfAtoms;
class ResidueHeader;
class CollectionData;
const int UNKNOWN_FILE = 0;
const int PDB_FILE = 1;
const int PDB1_FILE = 2;
bool strcmpCaseInsensitive(const char *str1, const char *str2);
class parserPDBtoGOA {
public:
  parserPDBtoGOA();
  virtual ~parserPDBtoGOA();
  GroupOfAtoms *parsePDBModels(const char *fileName, bool deleteWater = false);
  int mapElementAndResidue(char *atomName, char *residueName);
  int mapPQRElementAndResidue(const char *constatomName,
                              const char *const residueName);
  void initResidue(GroupOfAtoms *residue, ResidueHeader *resHeader,
                   int residueID, int resNum);

protected:
  GroupOfAtoms *parsePDB(string filename, CollectionData *collectionData,
                         bool deleteWater = false);
  void getString(char *line, char *str, int pos, int len);
  void getChar(char *line, char *str, int pos);
  void getInt(char *line, int *i, int pos, int len);
  void getDouble(char *line, double *d, int pos, int len);
  char *trim(char *s);
  bool residueHasChanged(char chainID, int residueSeq, char iCode);
  int mapResidue(char *resName);

  // for parsing the header
  void parseHEADER(char *line, CollectionData *collectionData);
  void parseTITLE(char *line, CollectionData *collectionData);
  void parseCOMPND(char *line, CollectionData *collectionData);
  void parseSOURCE(char *line, CollectionData *collectionData);
  void parseKEYWDS(char *line, CollectionData *collectionData);
  void parseEXPDTA(char *line, CollectionData *collectionData);
  void parseAUTHOR(char *line, CollectionData *collectionData);
  void parseJRNL(char *line, CollectionData *collectionData);
  void parseREMARK(char *line, CollectionData *collectionData);
  void parseREMARK350(char *line, CollectionData *collectionData);
  void getChainIds(char *line, std::vector<char> *chainIds);

  // for parsing the rest of the molecule
  void parseMODEL(GroupOfAtoms *molecule, char *line);
  void parseENDMDL(GroupOfAtoms *molecule, char *line);
  void parseSEQADV(GroupOfAtoms *molecule, char *line);
  void parseSEQRES(GroupOfAtoms *molecule, char *line);
  void parseHELIX(GroupOfAtoms *molecule, char *line);
  void parseSHEET(GroupOfAtoms *molecule, char *line);
  void parseTURN(GroupOfAtoms *molecule, char *line);
  void parseATOM(GroupOfAtoms *molecule, char *line, bool boundaryAtom,
                 CollectionData *collectionData, bool deleteWater = false);
  void finalizeParsing(GroupOfAtoms *molecule);
  // CHA: this is to check PDB parser
  void showGOAhierarchy(GroupOfAtoms *molecule, int level);
  // CHA: parsing hbond file
  void parseHBOND(GroupOfAtoms *molecule, const char *fname);
  void flattenAtom(PDBParser::GroupOfAtoms *molecule, Atom **atomlist,
                   int maxIdx);

  bool addTransformationsIfNeeded(GroupOfAtoms *molecule);
  bool getFileType(const char *const fileName);
  void getChainIds(GroupOfAtoms *molecule, std::vector<char> *chainIds);
  double xmin, ymin, zmin, xmax, ymax, zmax;
  int *residueArray;
  int residueAssigned;
  int chainTable[256]; // maps chain Id to the first residue number in that
                       // chain

  // Used to store all of the helices, sheets, and turns
  ExpandableArray<GroupOfAtoms> chainGroups;
  ExpandableArray<GroupOfAtoms> helixGroups;
  ExpandableArray<GroupOfAtoms> sheetGroups;
  ExpandableArray<GroupOfAtoms> turnGroups;
  int helixGroupNext;
  int sheetGroupNext;
  int turnGroupNext;

  // All residues that are not in an SS are placed in here
  ExpandableArray<GroupOfAtoms> residueGroups;

  // Temporaries when building the multiresolution tree
  GroupOfAtoms *currResidue;
  GroupOfAtoms *currSS;
  GroupOfAtoms *currChain;
  int residueCheckCount;
  char chainIDLast;   // The last chain ID seen
  int residueSeqLast; // The last residue sequence number seen
  char iCodeLast;     // The last insertion code seen
  double sasThreshold;
  int m_FileType;
  GroupOfAtoms
      *currMolecule; // The current molecule being built from the PDB file
  int m_NumberOfResiduesInCurrentChain;

  // Keeps track of all models (aka entire molecules) within the PDB file
  ExpandableArray<GroupOfAtoms> molModels;
};

}; // namespace PDBParser

#endif
