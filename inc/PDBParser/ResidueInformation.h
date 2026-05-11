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
// ResidueInformation.h: information for residue i/o
#ifndef _RESIDUE_INFORMATION_H_
#define _RESIDUE_INFORMATION_H_

#include <Utility/utility.h>

namespace PDBParser {
static const int TOTAL_RES_TYPE =
    26; // 25 residues total + 1 default , the 0th one
static const int RES_TYPE_DEFAULT = 0; // index of the default residue
static const char RES_TYPE_DEFAULT_SHORT_ID =
    'X';                                        // index of the default residue
static const char *RES_TYPE_DEFAULT_ID = "XXX"; // index of the default residue
static const int RES_TYPE_FIRST_NUCLEOTIDE =
    22; // index of the first nucleotide (i.e. A, C, G, T, or U)
static const int TOTAL_AA_RES_TYPE = 21; // includes the default "XXX"
static const int AA_PSI = 0;
static const int AA_PHI = 1;
static const int AA_CHI1 = 2;
static const int AA_CHI2 = 3;
static const int AA_CHI3 = 4;
static const int AA_CHI4 = 5;
static const int AA_CHI5 = 6;

typedef struct _ResidueInformation {
  char residueName[4];
  char shortResidueName;
  double hydrophobicity;
} ResidueInformation;

static ResidueInformation residueTable[TOTAL_RES_TYPE] = {
    {"XXX", 'X', 0.000},  {"ALA", 'A', -0.232}, {"ARG", 'R', 1.000},
    {"ASN", 'N', 0.528},  {"ASP", 'D', 0.944},  {"CYS", 'C', -0.360},
    {"GLN", 'Q', 0.498},  {"GLU", 'E', 0.914},  {"GLY", 'G', -0.002},
    {"HIS", 'H', 0.670},  {"ILE", 'I', -0.886}, {"LEU", 'L', -0.886},
    {"LYS", 'K', 0.434},  {"MET", 'M', -0.676}, {"PHE", 'F', -1.000},
    {"PRO", 'P', -0.422}, {"SER", 'S', 0.282},  {"THR", 'T', 0.100},
    {"TRP", 'W', -0.756}, {"TYR", 'Y', -0.760}, {"VAL", 'V', -0.650},
    {"  A", 'a', -3.225}, {"  C", 'c', -3.262}, {"  G", 'g', -4.210},
    {"  T", 't', -3.286}, {"  U", 'u', -3.286}};

// Default colors for residues

// The red value for the colors of 1 default, the twenty three residues + 5
// nucleic acid residues
static const double resColorR[TOTAL_RES_TYPE] = {
    /* Unknown */ 0.9,
    /* Amino Acids */
    0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
    0.6, 0.9, 0.9, 0.9, 0.9,
    /* Nucleotides */ 0.6, 0.9, 0.9, 0.6, 0.6};

// The green value for the colors of 1 default, the twenty three residues + 5
// nucleic acid residues
static const double resColorG[TOTAL_RES_TYPE] = {
    /* Unknown */ 0.9,
    /* Amino Acids */
    0.3, 0.3, 0.3, 0.6, 0.6, 0.9, 0.9, 0.3, 0.3, 0.3, 0.6, 0.6, 0.6, 0.9, 0.9,
    0.9, 0.3, 0.3, 0.3, 0.6,
    /* Nucleotides */ 0.6, 0.5, 0.3, 0.9, 0.6};

// The blue value for the colors of  1 default, the twenty three residues + 5
// nucleic acid residues
static const double resColorB[TOTAL_RES_TYPE] = {
    /* Unknown */ 0.9,
    /* Amino Acids */
    0.3, 0.6, 0.9, 0.3, 0.9, 0.3, 0.6, 0.3, 0.6, 0.9, 0.3, 0.6, 0.9, 0.3, 0.6,
    0.9, 0.3, 0.6, 0.9, 0.6,
    /* Nucleotides */ 0.9, 0.3, 0.3, 0.6, 0.6};

// Different atoms used in the definition of Chi angles for amino acids
static const char *const AminoAcidChi1Table[TOTAL_AA_RES_TYPE][5] = {
    {"XXX", "    ", "    ", "    ", "    "},

    {"ALA", "    ", "    ", "    ", "    "},
    {"ARG", " N  ", " CA ", " CB ", " CG "},
    {"ASN", " N  ", " CA ", " CB ", " CG "},
    {"ASP", " N  ", " CA ", " CB ", " CG "},

    {"CYS", " N  ", " CA ", " CB ", " SG "},
    {"GLN", " N  ", " CA ", " CB ", " CG "},
    {"GLU", " N  ", " CA ", " CB ", " CG "},
    {"GLY", "    ", "    ", "    ", "    "},

    {"HIS", " N  ", " CA ", " CB ", " CG "},
    {"ILE", " N  ", " CA ", " CB ", " CG1"},
    {"LEU", " N  ", " CA ", " CB ", " CG "},
    {"LYS", " N  ", " CA ", " CB ", " CG "},

    {"MET", " N  ", " CA ", " CB ", " CG "},
    {"PHE", " N  ", " CA ", " CB ", " CG "},
    {"PRO", " N  ", " CA ", " CB ", " CG "},
    {"SER", " N  ", " CA ", " CB ", " OG "},

    {"THR", " N  ", " CA ", " CB ", " OG1"},
    {"TRP", " N  ", " CA ", " CB ", " CG "},
    {"TYR", " N  ", " CA ", " CB ", " CG "},
    {"VAL", " N  ", " CA ", " CB ", " CG1"}};

static const char *const AminoAcidChi2Table[TOTAL_AA_RES_TYPE][5] = {
    {"XXX", "    ", "    ", "    ", "    "},

    {"ALA", "    ", "    ", "    ", "    "},
    {"ARG", " CA ", " CB ", " CG ", " CD "},
    {"ASN", " CA ", " CB ", " CG ", " OD1"},
    {"ASP", " CA ", " CB ", " CG ", " OD1"},

    {"CYS", "    ", "    ", "    ", "    "},
    {"GLN", " CA ", " CB ", " CG ", " CD "},
    {"GLU", " CA ", " CB ", " CG ", " CD "},
    {"GLY", "    ", "    ", "    ", "    "},

    {"HIS", " CA ", " CB ", " CG ", " ND1"},
    {"ILE", " CA ", " CB ", " CG1", " CD "},
    {"LEU", " CA ", " CB ", " CG ", " CD1"},
    {"LYS", " CA ", " CB ", " CG ", " CD "},

    {"MET", " CA ", " CB ", " CG ", " SD "},
    {"PHE", " CA ", " CB ", " CG ", " CD1"},
    {"PRO", " CA ", " CB ", " CG ", " CD "},
    {"SER", "    ", "    ", "    ", "    "},

    {"THR", "    ", "    ", "    ", "    "},
    {"TRP", " CA ", " CB ", " CG ", " CD1"},
    {"TYR", " CA ", " CB ", " CG ", " CD1"},
    {"VAL", "    ", "    ", "    ", "    "}};

static const char *const AminoAcidChi3Table[TOTAL_AA_RES_TYPE][5] = {
    {"  X", "    ", "    ", "    ", "    "},

    {"ALA", "    ", "    ", "    ", "    "},
    {"ARG", " CB ", " CG ", " CD ", " NE "},
    {"ASN", "    ", "    ", "    ", "    "},
    {"ASP", "    ", "    ", "    ", "    "},

    {"CYS", "    ", "    ", "    ", "    "},
    {"GLN", " CB ", " CG ", " CD ", " OE1"},
    {"GLU", " CB ", " CG ", " CD ", " OE1"},
    {"GLY", "    ", "    ", "    ", "    "},

    {"HIS", "    ", "    ", "    ", "    "},
    {"ILE", "    ", "    ", "    ", "    "},
    {"LEU", "    ", "    ", "    ", "    "},
    {"LYS", " CB ", " CG ", " CD ", " CE "},

    {"MET", " CB ", " CG ", " SD ", " CE "},
    {"PHE", "    ", "    ", "    ", "    "},
    {"PRO", "    ", "    ", "    ", "    "},
    {"SER", "    ", "    ", "    ", "    "},

    {"THR", "    ", "    ", "    ", "    "},
    {"TRP", "    ", "    ", "    ", "    "},
    {"TYR", "    ", "    ", "    ", "    "},
    {"VAL", "    ", "    ", "    ", "    "}};

static const char *const AminoAcidChi4Table[TOTAL_AA_RES_TYPE][5] = {
    {"  X", "    ", "    ", "    ", "    "},

    {"ALA", "    ", "    ", "    ", "    "},
    {"ARG", " CG ", " CD ", " NE ", " CZ "},
    {"ASN", "    ", "    ", "    ", "    "},
    {"ASP", "    ", "    ", "    ", "    "},

    {"CYS", "    ", "    ", "    ", "    "},
    {"GLN", "    ", "    ", "    ", "    "},
    {"GLU", "    ", "    ", "    ", "    "},
    {"GLY", "    ", "    ", "    ", "    "},

    {"HIS", "    ", "    ", "    ", "    "},
    {"ILE", "    ", "    ", "    ", "    "},
    {"LEU", "    ", "    ", "    ", "    "},
    {"LYS", " CG ", " CD ", " CE ", " NZ "},

    {"MET", "    ", "    ", "    ", "    "},
    {"PHE", "    ", "    ", "    ", "    "},
    {"PRO", "    ", "    ", "    ", "    "},
    {"SER", "    ", "    ", "    ", "    "},

    {"THR", "    ", "    ", "    ", "    "},
    {"TRP", "    ", "    ", "    ", "    "},
    {"TYR", "    ", "    ", "    ", "    "},
    {"VAL", "    ", "    ", "    ", "    "}};

static const char *const AminoAcidChi5Table[TOTAL_AA_RES_TYPE][5] = {
    {"  X", "    ", "    ", "    ", "    "},

    {"ALA", "    ", "    ", "    ", "    "},
    {"ARG", " CD ", " NE ", " CZ ", " NH1"},
    {"ASN", "    ", "    ", "    ", "    "},
    {"ASP", "    ", "    ", "    ", "    "},

    {"CYS", "    ", "    ", "    ", "    "},
    {"GLN", "    ", "    ", "    ", "    "},
    {"GLU", "    ", "    ", "    ", "    "},
    {"GLY", "    ", "    ", "    ", "    "},

    {"HIS", "    ", "    ", "    ", "    "},
    {"ILE", "    ", "    ", "    ", "    "},
    {"LEU", "    ", "    ", "    ", "    "},
    {"LYS", "    ", "    ", "    ", "    "},

    {"MET", "    ", "    ", "    ", "    "},
    {"PHE", "    ", "    ", "    ", "    "},
    {"PRO", "    ", "    ", "    ", "    "},
    {"SER", "    ", "    ", "    ", "    "},

    {"THR", "    ", "    ", "    ", "    "},
    {"TRP", "    ", "    ", "    ", "    "},
    {"TYR", "    ", "    ", "    ", "    "},
    {"VAL", "    ", "    ", "    ", "    "}};

static int relevantNumberOfAtoms[TOTAL_AA_RES_TYPE] = {
    4, 5, 11, 8, 8, 6, 9, 9, 4, 10, 8, 8, 9, 8, 11, 7, 6, 7, 14, 12, 7};

static const char
    *const AminoAcidResidueContentsOrdering[TOTAL_AA_RES_TYPE][15] = {
        {"  X", " N  ", " CA ", " C  ", " O  ", "    ", "    ", "    ", "    ",
         "    ", "    ", "    ", "    ", "    ", "    "},

        {"ALA", " N  ", " CA ", " C  ", " O  ", " CB ", "    ", "    ", "    ",
         "    ", "    ", "    ", "    ", "    ", "    "},
        {"ARG", " N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD ", " NE ",
         " CZ ", " NH1", " NH2", "    ", "    ", "    "},
        {"ASN", " N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " OD1", " ND2",
         "    ", "    ", "    ", "    ", "    ", "    "},
        {"ASP", " N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " OD1", " OD2",
         "    ", "    ", "    ", "    ", "    ", "    "},

        {"CYS", " N  ", " CA ", " C  ", " O  ", " CB ", " SG ", "    ", "    ",
         "    ", "    ", "    ", "    ", "    ", "    "},
        {"GLN", " N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD ", " OE1",
         " NE2", "    ", "    ", "    ", "    ", "    "},
        {"GLU", " N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD ", " OE1",
         " OE2", "    ", "    ", "    ", "    ", "    "},
        {"GLY", " N  ", " CA ", " C  ", " O  ", "    ", "    ", "    ", "    ",
         "    ", "    ", "    ", "    ", "    ", "    "},

        {"HIS", " N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " ND1", " CD2",
         " CE1", " NE2", "    ", "    ", "    ", "    "},
        {"ILE", " N  ", " CA ", " C  ", " O  ", " CB ", " CG1", " CG2", " CD1",
         "    ", "    ", "    ", "    ", "    ", "    "},
        {"LEU", " N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD1", " CD2",
         "    ", "    ", "    ", "    ", "    ", "    "},
        {"LYS", " N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD ", " CE ",
         " NZ ", "    ", "    ", "    ", "    ", "    "},

        {"MET", " N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " SD ", " CE ",
         "    ", "    ", "    ", "    ", "    ", "    "},
        {"PHE", " N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD1", " CD2",
         " CE1", " CE2", " CZ ", "    ", "    ", "    "},
        {"PRO", " N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD ", "    ",
         "    ", "    ", "    ", "    ", "    ", "    "},
        {"SER", " N  ", " CA ", " C  ", " O  ", " CB ", " OG ", "    ", "    ",
         "    ", "    ", "    ", "    ", "    ", "    "},

        {"THR", " N  ", " CA ", " C  ", " O  ", " CB ", " OG1", " CG2", "    ",
         "    ", "    ", "    ", "    ", "    ", "    "},
        {"TRP", " N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD1", " CD2",
         " NE1", " CE2", " CE3", " CZ2", " CZ3", " CH2"},
        {"TYR", " N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD1", " CD2",
         " CE1", " CE2", " CZ ", " OH ", "    ", "    "},
        {"VAL", " N  ", " CA ", " C  ", " O  ", " CB ", " CG1", " CG2", "    ",
         "    ", "    ", "    ", "    ", "    ", "    "}};
}; // namespace PDBParser

#endif
