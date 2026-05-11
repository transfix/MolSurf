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
/*   This file contains some extra parsing routines, dealing exclusively   */
/*   with the PDB header. Since it is quite long and complicated, I have   */
/*   divided the parser into 2 c++ files.                                  */
/*   See parserPDBtoGOA.cpp for the main part of the parser                */

#include <PDBParser/parserPDBtoGOA.h>
#include <PDBParser/CollectionData.h>
#include <PDBParser/BiomoleculeTransformations.h>

using namespace PDBParser;

void parserPDBtoGOA::parseHEADER(char *line, CollectionData *collectionData) {
  if (!collectionData || !line) {
    return;
  }
  /*   1 -  6        Record name     "HEADER"

          11 - 50        String(40)      classification  Classifies the
     molecule(s)

          51 - 59        Date            depDate         Deposition date.  This
     is the date the coordinates were received by the PDB

          63 - 66        IDcode          idCode          This identifier is
     unique within PDB
  */
  getString(line, collectionData->m_Classification, 10, 40);
  getString(line, collectionData->m_DepDate, 50, 9);
  getString(line, collectionData->m_PDBID, 62, 4);
}

void parserPDBtoGOA::parseTITLE(char *line, CollectionData *collectionData) {
  if (!collectionData || !line) {
    return;
  }
  /*
           1 -  6        Record name     "TITLE "
           9 - 10        Continuation    continuation   Allows concatenation of
     multiple records. 11 - 70        String          title          Title of
     the experiment.
  */
  if (!collectionData->m_Title) {
    collectionData->m_Title = new char[61];
    getString(line, collectionData->m_Title, 10, 60);
  } else {
    int curLen = strlen(collectionData->m_Title);
    char *temp = new char[curLen + 1];
    strcpy(temp, collectionData->m_Title);
    delete[] collectionData->m_Title;
    collectionData->m_Title = new char[curLen + 62];
    char extension[61];
    getString(line, extension, 10, 60);
    strcpy(collectionData->m_Title, temp);
    strcat(collectionData->m_Title, "\n");
    strcat(collectionData->m_Title, extension);
    delete[] temp;
  }
}

void parserPDBtoGOA::parseCOMPND(char *line, CollectionData *collectionData) {
  if (!collectionData || !line) {
    return;
  }
  /*
           1 -  6        Record name       "COMPND"
           9 - 10        Continuation      continuation   Allows concatenation
     of multiple records. 11 - 70        Specification     compound Description
     of the molecular components.
  */
  if (!collectionData->m_MoleculeDescription) {
    collectionData->m_MoleculeDescription = new char[61];
    getString(line, collectionData->m_MoleculeDescription, 10, 60);
  } else {
    int curLen = strlen(collectionData->m_MoleculeDescription);
    char *temp = new char[curLen + 1];
    strcpy(temp, collectionData->m_MoleculeDescription);
    delete[] collectionData->m_MoleculeDescription;
    collectionData->m_MoleculeDescription = new char[curLen + 62];
    char extension[61];
    getString(line, extension, 10, 60);
    strcpy(collectionData->m_MoleculeDescription, temp);
    strcat(collectionData->m_MoleculeDescription, "\n");
    strcat(collectionData->m_MoleculeDescription, extension);
    delete[] temp;
  }
}

void parserPDBtoGOA::parseSOURCE(char *line, CollectionData *collectionData) {
  if (!collectionData || !line) {
    return;
  }
  /*
           1 -  6        Record name       "SOURCE"
           9 - 10        Continuation      continuation   Allows concatenation
     of multiple records. 11 - 70   Specification list     srcName Identifies
     the source of the macromolecule in a token:
  */
  if (!collectionData->m_Source) {
    collectionData->m_Source = new char[61];
    getString(line, collectionData->m_Source, 10, 60);
  } else {
    int curLen = strlen(collectionData->m_Source);
    char *temp = new char[curLen + 1];
    strcpy(temp, collectionData->m_Source);
    delete[] collectionData->m_Source;
    collectionData->m_Source = new char[curLen + 62];
    char extension[61];
    getString(line, extension, 10, 60);
    strcpy(collectionData->m_Source, temp);
    strcat(collectionData->m_Source, "\n");
    strcat(collectionData->m_Source, extension);
    delete[] temp;
  }
}

void parserPDBtoGOA::parseKEYWDS(char *line, CollectionData *collectionData) {
  if (!collectionData || !line) {
    return;
  }
  /*
           1 -  6        Record name     "KEYWDS"
           9 - 10        Continuation    continuation   Allows concatenation of
     records if necessary. 11 - 70        List            keywds Comma-separated
     list of keywords relevant to the entry.
  */
  if (!collectionData->m_Keywords) {
    collectionData->m_Keywords = new char[61];
    getString(line, collectionData->m_Keywords, 10, 60);
  } else {
    int curLen = strlen(collectionData->m_Keywords);
    char *temp = new char[curLen + 1];
    strcpy(temp, collectionData->m_Keywords);
    delete[] collectionData->m_Keywords;
    collectionData->m_Keywords = new char[curLen + 62];
    char extension[61];
    getString(line, extension, 10, 60);
    strcpy(collectionData->m_Keywords, temp);
    strcat(collectionData->m_Keywords, "\n");
    strcat(collectionData->m_Keywords, extension);
    delete[] temp;
  }
}

void parserPDBtoGOA::parseEXPDTA(char *line, CollectionData *collectionData) {
  if (!collectionData || !line) {
    return;
  }
  /*
           1 -  6       Record name    "EXPDTA"
           9 - 10       Continuation   continuation  Allows concatenation of
     multiple records. 11 - 70       SList          technique     The
     experimental technique(s) with optional comment describing the sample or
     experiment.
  */
  if (!collectionData->m_ExperimentType) {
    collectionData->m_ExperimentType = new char[61];
    getString(line, collectionData->m_ExperimentType, 10, 60);
  } else {
    int curLen = strlen(collectionData->m_ExperimentType);
    char *temp = new char[curLen + 1];
    strcpy(temp, collectionData->m_ExperimentType);
    delete[] collectionData->m_ExperimentType;
    collectionData->m_ExperimentType = new char[curLen + 62];
    char extension[61];
    getString(line, extension, 10, 60);
    strcpy(collectionData->m_ExperimentType, temp);
    strcat(collectionData->m_ExperimentType, "\n");
    strcat(collectionData->m_ExperimentType, extension);
    delete[] temp;
  }
}

void parserPDBtoGOA::parseAUTHOR(char *line, CollectionData *collectionData) {
  if (!collectionData || !line) {
    return;
  }
  /*
           1 -  6       Record name    "AUTHOR"
           9 - 10       Continuation   continuation  Allows concatenation of
     multiple records. 11 - 70       List           authorList    List of the
     author names, separated by commas.
  */
  if (!collectionData->m_Authors) {
    collectionData->m_Authors = new char[61];
    getString(line, collectionData->m_Authors, 10, 60);
  } else {
    int curLen = strlen(collectionData->m_Authors);
    char *temp = new char[curLen + 1];
    strcpy(temp, collectionData->m_Authors);
    delete[] collectionData->m_Authors;
    collectionData->m_Authors = new char[curLen + 62];
    char extension[61];
    getString(line, extension, 10, 60);
    strcpy(collectionData->m_Authors, temp);
    strcat(collectionData->m_Authors, "\n");
    strcat(collectionData->m_Authors, extension);
    delete[] temp;
  }
}

void parserPDBtoGOA::parseJRNL(char *line, CollectionData *collectionData) {
  if (!collectionData || !line) {
    return;
  }
  /*
           1 -  6    Record name    "JRNL  "
          13 - 70    LString        text      See Details below.
  */
  if (!collectionData->m_Citation) {
    collectionData->m_Citation = new char[59];
    getString(line, collectionData->m_Citation, 12, 58);
  } else {
    int curLen = strlen(collectionData->m_Citation);
    char *temp = new char[curLen + 1];
    strcpy(temp, collectionData->m_Citation);
    delete[] collectionData->m_Citation;
    collectionData->m_Citation = new char[curLen + 62];
    char extension[59];
    getString(line, extension, 12, 58);
    strcpy(collectionData->m_Citation, temp);
    strcat(collectionData->m_Citation, "\n");
    strcat(collectionData->m_Citation, extension);
    delete[] temp;
  }
}

void parserPDBtoGOA::parseREMARK(char *line, CollectionData *collectionData) {
  if (!collectionData || !line) {
    return;
  }
  /*
           1 -  6      Record name     "REMARK"
           8 - 10      Integer         remarkNum      Remark number. It is not
     an error for remark n to exist in an entry when remark n-1 does not. 12 -
     70      LString         empty          Left as white space in first line of
                                                              each new remark.
  */
  int remarkNumber = -1;
  getInt(line, &remarkNumber, 7, 3);
  char remark[80];
  getString(line, remark, 11,
            69); // remark 350 can take characters till column 80
  if (remarkNumber == 350) {
    parseREMARK350(remark, collectionData);
  }
}

// The input is a sequence of characters separated by commas.
// That is about as specific as the PDB file format gets.
void parserPDBtoGOA::getChainIds(char *line, std::vector<char> *chainIds) {
  int len = strlen(line);
  int i;
  bool addedSomething = false;
  for (i = 0; i < len - 1; i++) {
    if (line[i] != ' ' &&
        line[i] != ',') // not checking for alternating commas and charecters.
    {
      chainIds->push_back(line[i]);
      addedSomething = true;
    }
  }
  // what else to do ? Some idiot wrote the file format.
  if (!addedSomething) {
    chainIds->push_back(' ');
  }
}

// We only send the required fields. Namely the string from line 12 - 70 ( index
// starting from 1 ), since the previous stuff is simply REMARK <Remark number>.
void parserPDBtoGOA::parseREMARK350(char *line,
                                    CollectionData *collectionData) {
  if (!collectionData || !line) {
    return;
  }
  /*
                           1         2         3         4         5         6 7
          1234567890123456789012345678901234567890123456789012345678901234567890
          REMARK 350
          REMARK 350 GENERATING THE BIOMOLECULE
          REMARK 350 COORDINATES FOR A COMPLETE MULTIMER REPRESENTING THE KNOWN
          REMARK 350 BIOLOGICALLY SIGNIFICANT OLIGOMERIZATION STATE OF THE
          REMARK 350 MOLECULE CAN BE GENERATED BY APPLYING BIOMT TRANSFORMATIONS
          REMARK 350 GIVEN BELOW.  BOTH NON-CRYSTALLOGRAPHIC AND
          REMARK 350 CRYSTALLOGRAPHIC OPERATIONS ARE GIVEN.
          REMARK 350
          REMARK 350 APPLY THE FOLLOWING TO CHAINS: ?, ?...
          REMARK 350   BIOMT1   N  N.NNNNNN  N.NNNNNN  N.NNNNNN        N.NNNNN
          REMARK 350   BIOMT2   N  N.NNNNNN  N.NNNNNN  N.NNNNNN        N.NNNNN
          REMARK 350   BIOMT3   N  N.NNNNNN  N.NNNNNN  N.NNNNNN        N.NNNNN
  */
  // We need to keep state information to see which remark we are parsing, and
  // also the transformation matrix row number ( [1..3] ). We could store this
  // information in the collectionData itself. So no ugly static / global
  // variables are needed. Get rid of anything useless.
  if (!strstr(line, "APPLY THE FOLLOWING TO CHAINS:") &&
      !strstr(line, "  BIOMT")) {
    return;
  }
  // This can repeat multiple times, even on successive lines!
  if (strstr(line, "APPLY THE FOLLOWING TO CHAINS:")) {
    // parse chains
    std::vector<char> chainIds;
    getChainIds(&(line[30]), &chainIds);
    if (chainIds.size() < 1) {
      return; // nothing to add.
    }
    // if we dont have a transformation yet, keep adding to existing chains
    // or create a new BiomoleculeTransformations and make this the new
    // list of chains.
    int n = collectionData->m_BiomoleculeTransformations.size();
    BiomoleculeTransformations *b = 0;
    if (n > 0) {
      if (collectionData->m_BiomoleculeTransformations.at(n - 1)) {
        // assume that no transformations have been added.
        if (!((collectionData->m_BiomoleculeTransformations.at(n - 1))
                  ->rotation.size())) {
          b = collectionData->m_BiomoleculeTransformations.at(n - 1);
        }
      } else {
        // there is an empty pointer ! Leave it empty. Should not happen
        // though...
      }
    } else {
      // either no existing transformations, or we have to enter chain
      // information in to a new tuple
      b = new BiomoleculeTransformations();
      collectionData->m_BiomoleculeTransformations.push_back(b);
    }
    if (!b) {
      return;
    }
    b->addChains(chainIds);
  } else {
    int n = collectionData->m_BiomoleculeTransformations.size();
    if (n < 1) {
      return;
    }
    // add transformation to latest BiomoleculeTransformations.
    // if they do wierd things like repeat BIOMT1 for the same
    // transformation, this code will just overwrite the
    // existing stuff.
    // If there is no existing BiomoleculeTransformations, dont
    // do anything.
    // if N is larger than current+1, bad luck.
    //  BIOMTX   N  N.NNNNNN  N.NNNNNN  N.NNNNNN        N.NNNNN
    double r1 = 0, r2 = 0, r3 = 0, t = 0;
    getDouble(line, &r1, 12, 10);
    getDouble(line, &r2, 22, 10);
    getDouble(line, &r3, 32, 10);
    getDouble(line, &t, 42, 15);
    // we ignore the transformationNumber. If there is some transformation
    // missing, im not going to worry. If there are not enough rotations,
    // translations etc in the end, the code using it should make sure things
    // dont crash. I also ignore the BIOMTX value. Whats the use of knowing
    // there is something wrong with the input file ? Might as well shut up and
    // carry on in life. add new values
    BiomoleculeTransformations *b =
        collectionData->m_BiomoleculeTransformations.at(n - 1);
    if (!b) {
      return;
    }
    b->rotation.push_back(r1);
    b->rotation.push_back(r2);
    b->rotation.push_back(r3);
    b->translation.push_back(t);
  }
}
