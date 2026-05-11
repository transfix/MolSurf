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
#include <PDBParser/GOAFileIO.h>

using namespace PDBParser;

static const char NULL_ELEM_CHAR = '-';

// GroupOfAtoms file format

// Given a positive integer INDENT, prints (2 * INDENT) blank spaces.
void printIndent(int indent) {
  for (int i = 0; i < indent; i++) {
    printf("  ");
  }
}

// Recursively print the contents of the current GroupofAtoms.  INDENT is the
// level of indentation.
void PDBParser::printGOA(GroupOfAtoms *currGOA, int indent) {
  static int total_atoms = 0;
  static int total_res = 0;
  static int total_ss = 0;
  static int total_chain = 0;
  int i, n;
  switch (currGOA->type) {
  case GENERAL_TYPE:
    printIndent(indent);
    printf("GENERAL MOLECULE\n");
    break;
  case PROTEIN_TYPE:
    printIndent(indent);
    printf("PROTEIN\n");
    break;
  case CHAIN_TYPE:
    printf("\n\n");
    printIndent(indent);
    printf("CHAIN:  ID: %c\n", (char)currGOA->resID);
    total_chain++;
    break;
  case SS_TYPE:
    printf("\n\n");
    printIndent(indent);
    printf("Secondary Structure: %s\n", ssString[currGOA->resID]);
    if (currGOA->resID != SS_NULL) {
      total_ss++; // SS_NULL is just a placeholder -- it's not a real SS
    }
    break;
  case RESIDUE_TYPE:
    // SKVINAY : print out the angles here.
    printf("\n");
    printIndent(indent);
    if (currGOA->resID > -1) {
      printf("Residue Type: %s   Chain: %c   Seq Num:%4d%c   Position: %d\n",
             residueTable[currGOA->resID].residueName,
             currGOA->resHeader->chainID, currGOA->resHeader->seqNum,
             currGOA->resHeader->iCode, currGOA->resPos);
    } else {
      printf("Residue Type: ???   Chain: %c   Seq Num:%4d%c   Position: %d\n",
             currGOA->resHeader->chainID, currGOA->resHeader->seqNum,
             currGOA->resHeader->iCode, currGOA->resPos);
    }
    total_res++;
    break;
  case COLLECTION_TYPE:
    // Print all of the molecules in the collection (a collection is the
    // hierarchical level right above a molecule)
    n = currGOA->m_SubGroups.size();
    for (i = 0; i < n; i++) {
      printIndent(indent);
      printf("\n********************* Begin Model *********************\n");
      printGOA(currGOA->m_SubGroups.at(i), indent);
      printIndent(indent);
      printf("\n********************** End Model **********************\n\n\n");
    }
    // For a collection, we do not need to print out the subarrays since we have
    // already explicitly iterated through them
    return;
  default:
    printIndent(indent);
    printf("Unknown_type\n");
  }
  // Iterate through the atoms array
  int atom_num = currGOA->m_Atoms.size();
  if (atom_num > 0) {
    printIndent(indent);
    printf("---- Component Atoms: ----\n");
  }
  for (i = 0; i < atom_num; i++) {
    Atom *currAtom = currGOA->m_Atoms.at(i);
    printIndent(indent);
    printf("Elem: %4.4s, Rad: %2.2lf, Loc: (%4.2lf, %4.2lf, %4.2lf)\n",
           currAtom->getName(), currAtom->getRadius(), currAtom->m_Position[0],
           currAtom->m_Position[1], currAtom->m_Position[2]);
    total_atoms++;
  }
  // Iterate through the subgroups array
  n = (int)currGOA->m_SubGroups.size();
  if (n > 0) {
    printIndent(indent);
    printf("---- Component Subgroups: ----\n");
  }
  for (i = 0; i < n; i++) {
    GroupOfAtoms *currSubGroup = currGOA->m_SubGroups.at(i);
    printGOA(currSubGroup, indent + 1);
  }
  // Print overall statistics at the top level of a molecule or model
  if (currGOA->type == GENERAL_TYPE || currGOA->type == PROTEIN_TYPE) {
    printf("\n\n");
    printIndent(indent);
    printf("Note: Statistics are cumulative across all models.\n");
    printIndent(indent);
    printf("Total Chains:               %d\n", total_chain);
    printIndent(indent);
    printf("Total Secondary Structures: %d\n", total_ss);
    printIndent(indent);
    printf("Total Residues:             %d\n", total_res);
    printIndent(indent);
    printf("Total Atoms:                %d\n", total_atoms);
  }
}

// Vorocomp PTS file format

// Count how many atoms are in the tree rooted at the given GroupOfAtoms
int countGOA2PTSRecursive(GroupOfAtoms *currGOA) {
  int num_atoms = 0;
  // Count the number of atoms on this level
  num_atoms += currGOA->m_Atoms.size();
  // Iterate through the subgroups array
  int n = (int)currGOA->m_SubGroups.size();
  for (int i = 0; i < n; i++) {
    GroupOfAtoms *currSubGroup = currGOA->m_SubGroups.at(i);
    num_atoms += countGOA2PTSRecursive(currSubGroup);
  }
  return num_atoms;
}

// The recursive helper call to writing the GOA to the file in the vorocomp .pts
// format
void writeGOA2PTSRecursive(FILE *stream, GroupOfAtoms *currGOA) {
  // Iterate through the atoms array
  int atom_num = currGOA->m_Atoms.size();
  for (int i = 0; i < atom_num; i++) {
    Atom *currAtom = currGOA->m_Atoms.at(i);
    fprintf(stream, "%.3lf %.3lf %.3lf %.3lf\n", currAtom->m_Position[0],
            currAtom->m_Position[1], currAtom->m_Position[2],
            currAtom->getRadius());
  }
  // Iterate through the subgroups array
  int n = (int)currGOA->m_SubGroups.size();
  for (int i = 0; i < n; i++) {
    GroupOfAtoms *currSubGroup = currGOA->m_SubGroups.at(i);
    writeGOA2PTSRecursive(stream, currSubGroup);
  }
}

// Write the GOA to the file in the vorocomp .pts format
// Return zero on success, non-zero on error.
int PDBParser::writeGOA2PTS(char *filename, GroupOfAtoms *currGOA) {
  FILE *stream;
  int num_atoms = 0;
  if ((stream = fopen(filename, "w")) == (FILE *)0) {
    // Handle file open error
    return 1;
  }
  num_atoms = countGOA2PTSRecursive(currGOA);
  fprintf(stream, "%d\n", num_atoms);
  writeGOA2PTSRecursive(stream, currGOA);
  fclose(stream);
  return 0;
}

int atom_serial_num; // Used to keep track of how many atoms have been numbered
int curr_res; // The current residue position (0 for first residue in pdb, 1 for
              // second residue, etc.)
int curr_res_in_chain; // The current residues position in the chain (1 for
                       // first residue in chain, 2 for second residue, etc.)
int curr_ss;           // The current secondary structure
int curr_chain;        // The current chain

// Print an Atom object to the stream in the style of the PDB's ATOM field.
void PDBParser::printAtomPDB(FILE *stream, Atom *atom,
                             CCVOpenGLMath::Matrix *transformation) {
  /*      1 -  6        Record name     "ATOM  "
           7 - 11        Integer         serial        Atom serial number.
     (right-justified) 13 - 16        Atom            name          Atom name.
     (2 char elem sym is right-justified) 17             Character       altLoc
     Alternate location indicator. 18 - 20        Residue name    resName
     Residue name. 22             Character       chainID       Chain
     identifier. 23 - 26        Integer         resSeq        Residue sequence
     number. (right justified) 27             AChar           iCode         Code
     for insertion of residues. (right justified) 31 - 38        Real(8.3) x
     Orthogonal coordinates for X in Angstroms. (right justified) 39 - 46
     Real(8.3)       y             Orthogonal coordinates for Y in Angstroms.
     (right justified) 47 - 54        Real(8.3)       z             Orthogonal
     coordinates for Z in Angstroms. (right justified) 55 - 60        Real(6.2)
     occupancy     Occupancy. (right justified) 61 - 66        Real(6.2)
     tempFactor    Temperature factor. (right justified) 73 - 76 LString(4)
     segID         Segment identifier, left-justified. 77 - 78        LString(2)
     element       Element symbol, right-justified. 79 - 80        LString(2)
     charge        Charge on the atom.*/
  float *loc;
  double x, y, z, r;
  loc = atom->getPosition();
  x = loc[0];
  y = loc[1];
  z = loc[2];
  // transform atom is transformation matrix is not null
  if (transformation) {
    CCVOpenGLMath::Vector newPoint =
        (*transformation) * CCVOpenGLMath::Vector(x, y, z, 1);
    x = newPoint[0];
    y = newPoint[1];
    z = newPoint[2];
  }
  r = atom->getRadius();
  // Atom serial number has to be 5 digits or less
  if (atom_serial_num >= 100000) {
    // fprintf(stderr, "error: atom serial number has overflowed!\n");
    printf("error: atom serial number has overflowed!\n");
    atom_serial_num = 1;
  }
  fprintf(stream, "ATOM"); // "ATOM "
  // debug: we need to use a string to get exactly 5 digits
  float tempf = (float)atom_serial_num;
  fprintf(stream, "%7.0f", tempf); // serial
  fprintf(stream, " ");            // char 12 is blank
  const char *atom_name = atom->getName();
  if (strlen(atom_name) >= 4) {
    fprintf(stream, "%-4.4s", atom->getName());
  } else {
    fprintf(stream, " %-3.3s", atom->getName());
  }
  fprintf(stream, " ");                           // altLoc
  fprintf(stream, "%3s", atom->getResidueName()); // resName
  fprintf(stream, " ");                           // char 21 is blank
  fprintf(stream, "%c", atom->getChain());
  // use this "commented out code" for residue-level lod
  // only works if there is one atom per residue!
  // basically, says each atom is in a residue by itself
  //     fprintf(stream, "%4d", (atom_serial_num%1000)); // resSeq
  // end "commented out code"
  fprintf(stream, "%4d", atom->getResidueNum()); // resSeq
  fprintf(stream, " ");                          // iCode
  fprintf(stream, "   ");                        // chars 28, 29, 30 are blank
  fprintf(stream, "%8.3f", x);                   // x
  fprintf(stream, "%8.3f", y);                   // y
  fprintf(stream, "%8.3f", z);                   // z
  fprintf(stream, "%6.2f", r);                   // occupancy
  fprintf(stream, "      ");                     // tempFactor
  fprintf(stream, "    ");                       // segID
  fprintf(stream, "      ");                     // chars 67-72 are blank
  fprintf(stream, "  ");                         // chars 67-72 are blank
  fprintf(stream, "  ");                         // charge
  fprintf(stream, "\n");
  atom_serial_num++;
}

// The recusive call that traverses the hierarchy.
void writeGOA2PDBRecursive(FILE *stream, GroupOfAtoms *currGOA, int type,
                           CCVOpenGLMath::Matrix *transformation) {

  int i, n;
  // Set the global for which chain we are in
  if (currGOA->type == CHAIN_TYPE) {
    curr_chain = currGOA->resID;
    curr_res_in_chain = 0; // reset the residue counter for this chain
  }
  if (currGOA->type == RESIDUE_TYPE) {
    curr_res = currGOA->resID;
    curr_res_in_chain++;
  }
  // If the current type is equal to the type, then output the mutiresolution
  // molecules and then return so that nothing lower in the tree will be
  // rendered.
  if (type == currGOA->type) {
    // For chains, print out the backbone atoms instead of the bounding balls.
    // It is done later, Hence continue if its CHAIN_TYPE
    if (currGOA->type != CHAIN_TYPE) {
      // Iterate through the multires array
      n = currGOA->m_MultiresAtoms.size();
      for (i = 0; i < n; i++) {
        Atom *currAtom = currGOA->m_MultiresAtoms[i];
        PDBParser::printAtomPDB(stream, currAtom, transformation);
      }
    }
    return;
  }
  // If we've reached a residue, write out the atoms no matter what.
  // Exclude HETATOMs when writing out the chain type, so that no atoms are
  // printed at this level except if control flow reaches a residue.
  if (currGOA->type == RESIDUE_TYPE) {
    n = currGOA->m_Atoms.size();
    for (i = 0; i < n; i++) {
      Atom *currAtom = currGOA->m_Atoms.at(i);
      if (type == CHAIN_TYPE) {
        if (currAtom->isBoundaryAtom()) {
          PDBParser::printAtomPDB(stream, currAtom, transformation);
        }
      } else {
        // pcd debug: don't print hetatoms for now (hetatoms are stored at the
        // molecule level's atom array
        PDBParser::printAtomPDB(stream, currAtom, transformation);
      }
    }
  }
  // Iterate through the subgroups array and descend down the tree.
  n = (int)currGOA->m_SubGroups.size();

  for (i = 0; i < n; i++) {
    GroupOfAtoms *currSubGroup = currGOA->m_SubGroups.at(i);
    writeGOA2PDBRecursive(stream, currSubGroup, type, transformation);
  }
}

// Writes the CURRGOA to the STREAM in PDB format at the resolution of TYPE
// (i.e. ATOM_TYPE to print out atoms, RESIDUE_TYPE to print out residues, etc,
// etc.).
bool PDBParser::writeGOA2PDB(FILE *stream, GroupOfAtoms *currGOA, int type,
                             CCVOpenGLMath::Matrix *transformation) {
  if (!stream || !currGOA) {
    return false;
  }
  // Initialize the counters
  atom_serial_num = 1;
  int init_serial_num = atom_serial_num;
  writeGOA2PDBRecursive(stream, currGOA, type, transformation);
  printf("Number of atoms written: %d\n", atom_serial_num - init_serial_num);
  return true;
}

void printAtomGOA(FILE *stream, Atom *atom) {
  float *loc;
  double x, y, z, r;
  loc = atom->getPosition();
  x = loc[0];
  y = loc[1];
  z = loc[2];
  r = atom->getRadius();
  fprintf(stream, "ATOM  "); // "ATOM "
  // Print the atom name
  const char *atom_name = atom->getName();
  if (strlen(atom_name) >= 4) {
    fprintf(stream, "%-4.4s", atom->getName());
  } else {
    fprintf(stream, " %-3.3s", atom->getName());
  }
  fprintf(stream, " ");        // iCode
  fprintf(stream, "%8.3f", x); // x
  fprintf(stream, "%8.3f", y); // y
  fprintf(stream, "%8.3f", z); // z
  fprintf(stream, "%6.2f", r); // occupancy
  fprintf(stream, "\n");
}

// The recusive call that traverses the hierarchy.
void writeGOARecursive(FILE *stream, GroupOfAtoms *currGOA, int chainType) {
  if (!stream || !currGOA) {
    return;
  }
  int i, n;
  if (currGOA->type == COLLECTION_TYPE) {
    if (currGOA->m_CollectionData) {
      currGOA->m_CollectionData->print(stream);
    }
  }
  // Set the global for which chain we are in
  if (currGOA->type == CHAIN_TYPE) {
    chainType = currGOA->resPos; // this variable is 'uglily' overloaded.
    curr_chain = currGOA->resID;
    curr_res_in_chain = 0; // reset the residue counter for this chain
    fprintf(stream, "Begin new chain %c\n", (char)curr_chain);
  }
  if (currGOA->type == SS_TYPE) {
    switch (currGOA->resID) {
    case SS_HELIX:
      fprintf(stream, "Begin new secondary structure: HELIX\n");
      break;
    case SS_SHEET:
      fprintf(stream, "Begin new secondary structure: SHEET\n");
      break;
    case SS_TURN:
      fprintf(stream, "Begin new secondary structure: TURN\n");
      break;
    case SS_NULL:
      fprintf(stream, "Begin new secondary structure: NULL\n");
      break;
    default:
      fprintf(stream, "Begin new secondary structure: NULL\n");
      break;
    }
  }
  if (currGOA->type == RESIDUE_TYPE) {
    curr_res = currGOA->resID;
    curr_res_in_chain++;
    fprintf(stream, "Begin new residue %3s\n",
            residueTable[curr_res].residueName);
  }
  // first print out the angles!
  if (currGOA->type == RESIDUE_TYPE) {
    if (chainType == CHAIN_PROTEIN) {
      if (!currGOA->m_ResidueData) {
        fprintf(stream, "Amino acid. Torsion angles not found\n");
      } else {
        fprintf(stream, "Amino acid. Torsion angles are : PHI %f, PSI %f\n",
                currGOA->m_ResidueData->m_Phi, currGOA->m_ResidueData->m_Psi);
      }
    } else if (chainType == CHAIN_NUCLEIC) {
      if (!currGOA->m_ResidueData) {
        fprintf(stream, "Nucleic acid. Torsion angles not found\n");
      } else {
        fprintf(
            stream,
            "Nucleic acid. Torsion angles are : ALPHA %f, BETA %f, GAMMA %f, "
            "DELTA %f, EPSILON %f, ZETA %f, NU0 %f, NU1 %f, NU2 %f, NU3 %f, "
            "NU4 %f, CHI %f, ETA %f, THETA %f\n",
            currGOA->m_ResidueData->m_alpha, currGOA->m_ResidueData->m_beta,
            currGOA->m_ResidueData->m_gamma, currGOA->m_ResidueData->m_delta,
            currGOA->m_ResidueData->m_epsilon, currGOA->m_ResidueData->m_zeta,
            currGOA->m_ResidueData->m_nu0, currGOA->m_ResidueData->m_nu1,
            currGOA->m_ResidueData->m_nu2, currGOA->m_ResidueData->m_nu3,
            currGOA->m_ResidueData->m_nu4, currGOA->m_ResidueData->m_chi,
            currGOA->m_ResidueData->m_Eta, currGOA->m_ResidueData->m_Theta);
      }
    } else {
      fprintf(stream, "Unknown residue type\n");
    }
  }
  // If we've reached a residue, write out the atoms no matter what.
  // Exclude HETATOMs when writing out the chain type, so that no atoms are
  // printed at this level except if control flow reaches a residue.
  //     if (currGOA->type == RESIDUE_TYPE || type != CHAIN_TYPE) {
  if (currGOA->type == RESIDUE_TYPE) {
    n = currGOA->m_Atoms.size();
    for (i = 0; i < n; i++) {
      Atom *currAtom = currGOA->m_Atoms.at(i);
      // pcd debug: don't print hetatoms for now (hetatoms are stored at the
      // molecule level's atom array
      printAtomGOA(stream, currAtom);
    }
  }
  // Iterate through the subgroups array and descend down the tree.
  n = (int)currGOA->m_SubGroups.size();
  for (i = 0; i < n; i++) {
    GroupOfAtoms *currSubGroup = currGOA->m_SubGroups.at(i);
    writeGOARecursive(stream, currSubGroup, chainType);
  }
  if (currGOA->type == CHAIN_TYPE) {
    fprintf(stream, "End chain\n");
  }
  if (currGOA->type == SS_TYPE) {
    fprintf(stream, "End secondary structure\n");
  }
  if (currGOA->type == RESIDUE_TYPE) {
    fprintf(stream, "End residue\n");
  }
}

bool PDBParser::writeGOAtoFile(FILE *stream, GroupOfAtoms *currGOA) {
  if (!stream || !currGOA) {
    return false;
  }
  int chainType = CHAIN_NULL;
  writeGOARecursive(stream, currGOA, chainType);
  return true;
}

void writeTorsionAnglesRecursive(FILE *stream, GroupOfAtoms *currGOA,
                                 int chainType) {
  if (!stream || !currGOA) {
    return;
  }
  int i, n;
  // Set the global for which chain we are in
  if (currGOA->type == CHAIN_TYPE) {
    chainType = currGOA->resPos; // this variable is 'uglily' overloaded.
  }
  // first print out the angles!
  if (currGOA->type == RESIDUE_TYPE) {
    if (chainType == CHAIN_PROTEIN) {
      if (!currGOA->m_ResidueData) {
        fprintf(stream, "Amino acid. Torsion angles not found\n");
      } else {
        fprintf(stream,
                "Amino acid. Torsion angles are : PHI % 7.5f, PSI % 7.5f\n",
                currGOA->m_ResidueData->m_Phi, currGOA->m_ResidueData->m_Psi);
      }
    } else if (chainType == CHAIN_NUCLEIC) {
      if (!currGOA->m_ResidueData) {
        fprintf(stream, "Nucleic acid. Torsion angles not found\n");
      } else {
        fprintf(
            stream,
            "Nucleic acid. Torsion angles are : ALPHA %f, BETA %f, GAMMA %f, "
            "DELTA %f, EPSILON %f, ZETA %f, NU0 %f, NU1 %f, NU2 %f, NU3 %f, "
            "NU4 %f, CHI %f, ETA %f, THETA %f\n",
            currGOA->m_ResidueData->m_alpha, currGOA->m_ResidueData->m_beta,
            currGOA->m_ResidueData->m_gamma, currGOA->m_ResidueData->m_delta,
            currGOA->m_ResidueData->m_epsilon, currGOA->m_ResidueData->m_zeta,
            currGOA->m_ResidueData->m_nu0, currGOA->m_ResidueData->m_nu1,
            currGOA->m_ResidueData->m_nu2, currGOA->m_ResidueData->m_nu3,
            currGOA->m_ResidueData->m_nu4, currGOA->m_ResidueData->m_chi,
            currGOA->m_ResidueData->m_Eta, currGOA->m_ResidueData->m_Theta);
      }
    } else {
      fprintf(stream, "Unknown residue type\n");
    }
  }
  // Iterate through the subgroups array and descend down the tree.
  n = (int)currGOA->m_SubGroups.size();
  for (i = 0; i < n; i++) {
    GroupOfAtoms *currSubGroup = currGOA->m_SubGroups.at(i);
    writeTorsionAnglesRecursive(stream, currSubGroup, chainType);
  }
}

bool PDBParser::writeTorsionAngles(FILE *stream, GroupOfAtoms *currGOA) {
  if (!stream || !currGOA) {
    return false;
  }
  int chainType = CHAIN_NULL;
  writeTorsionAnglesRecursive(stream, currGOA, chainType);
  return true;
}
