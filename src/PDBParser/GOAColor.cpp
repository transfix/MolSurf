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
#include <PDBParser/GOAColor.h>

using PDBParser::GOAColor;

GOAColor::GOAColor() {}

GOAColor::~GOAColor() {}

void GOAColor::parseColormapString(std::vector<std::string> strings) {
  char line[1024], *ptr = 0;
  unsigned int lineNum = 1;
  int sectionNum = 0;
  // int i;
  ColormapEntry cmapEnt;
  if (strings.size() == 0) {
    return;
  }
  {
    // clear the line buffer
    memset(line, 0, sizeof(line));
    // read the file into cmap
    // while (fgets(line, sizeof(line)-1, fp) != NULL)
    int s;
    for (s = 0; s < strings.size(); s++) {
      strcpy(line, strings[s].c_str());
      // skip over whitespace before checking for the opening quotation mark
      ptr = eatWhitespace(line);
      // determine if this line is a comment or an entry
      if (ptr[0] == '\"') {
        // this line is an entry
        char **tokens = parseEntryLine(line);
        // add the line number to the new entry
        cmapEnt.lineNum = lineNum;
        cmapEnt.transformSection = sectionNum;
        if (genEntryFromTokens(tokens, &cmapEnt)) {
          // add the newly parsed entry to "the list"
          m_EntList.push_back(cmapEnt);
        }
        freeTokens(tokens);
      } else if (ptr[0] == '[') {
        // this line marks the beginning of a transformation's section
        int nSection = parseSectionMarker(line);
        if (nSection > -1) {
          sectionNum = nSection;
        }
      } else if (ptr[0] != '#' && ptr[0] != '\0')
      // be verbose about misunderstood lines
      {
        printf("Warning: ignoring line %d\n", lineNum);
      }
      // else this line is a comment or is blank
      //  zero the line buffer
      memset(line, 0, sizeof(line));
      // increment the line number
      lineNum++;
    }
  }
}

void GOAColor::parseColormapFile(char *filename) {
  char line[1024], *ptr = 0;
  unsigned int lineNum = 1;
  int sectionNum = 0;
  ColormapEntry cmapEnt;
  FILE *fp;
  // open the file
  fp = fopen(filename, "r");
  if (fp != NULL) {
    // clear the line buffer
    memset(line, 0, sizeof(line));
    // read the file into cmap
    while (fgets(line, sizeof(line) - 1, fp) != NULL) {
      // skip over whitespace before checking for the opening quotation mark
      ptr = eatWhitespace(line);
      // determine if this line is a comment or an entry
      if (ptr[0] == '\"') {
        // this line is an entry
        char **tokens = parseEntryLine(line);
        // add the line number to the new entry
        cmapEnt.lineNum = lineNum;
        cmapEnt.transformSection = sectionNum;
        if (genEntryFromTokens(tokens, &cmapEnt)) {
          // add the newly parsed entry to "the list"
          m_EntList.push_back(cmapEnt);
        }
        freeTokens(tokens);
      } else if (ptr[0] == '[') {
        // this line marks the beginning of a transformation's section
        int nSection = parseSectionMarker(line);
        if (nSection > -1) {
          sectionNum = nSection;
        }
      } else if (ptr[0] != '#' && ptr[0] != '\0')
      // be verbose about misunderstood lines
      {
        printf("Warning: ignoring line %d\n", lineNum);
      }
      // else this line is a comment or is blank
      //  zero the line buffer
      memset(line, 0, sizeof(line));
      // increment the line number
      lineNum++;
    }
    // close the file
    fclose(fp);
  }
}

void GOAColor::printColormap() const {
  unsigned int i;
  for (i = 0; i < m_EntList.size(); i++) {
    printEntry(m_EntList[i]);
  }
}

void GOAColor::ApplyColormap(PDBParser::GroupOfAtoms *grp,
                             int transIndex) const {
  unsigned int i;
  printColormap();
  for (i = 0; i < m_EntList.size(); i++)
    if (m_EntList[i].transformSection == transIndex) {
      ApplyColormapEntryToGOA(m_EntList[i], grp);
    }
  // XXX: pseudo-hack for coloring HETATM atoms
  // in the current iteration of the parser, HETATM entries get put
  // in the m_Atoms array at the top level of a PROTEIN_TYPE GOA.
  if (grp->type == PDBParser::COLLECTION_TYPE) {
    grp = grp->m_SubGroups.at(0); // pick the first model. (it's a hack, ok?)
  }
  if (grp->type == PDBParser::PROTEIN_TYPE && grp->m_Atoms.size() > 0) {
    for (i = 0; i < grp->m_Atoms.size(); i++) {
      PDBParser::Atom *at = grp->m_Atoms.at(i);
      float rgb[3];
      float alpha;
      at->getColor(rgb + 0, rgb + 1, rgb + 2, &alpha);
      ApplyColorToAtom(rgb[0], rgb[1], rgb[2], 1.0, at);
    }
  }
}

void GOAColor::ApplyColorsByGOALevel(
    PDBParser::GroupOfAtoms *molecule,
    PDBParser::GroupOfAtoms::GOA_TYPE colorLevel) {
  switch (colorLevel) {
  case PDBParser::GroupOfAtoms::ATOM:
    ApplyElementColors(molecule);
    break;
  case PDBParser::GroupOfAtoms::RESIDUE:
    ApplyResidueColors(molecule);
    break;
  case PDBParser::GroupOfAtoms::SECONDARY_STRUCTURE:
    ApplyElementColors(molecule);
    break;
  case PDBParser::GroupOfAtoms::CHAIN: {
    ColorManager::ColorGenerator *pastelColorGenerator =
        new ColorManager::PastelColorGenerator();
    pastelColorGenerator->reset();
    ApplyChainColors(molecule, pastelColorGenerator);
    delete pastelColorGenerator;
  } break;
  case PDBParser::GroupOfAtoms::PROTEIN: // this is quite useless!
  {
    double red, green, blue, alpha = 1.0;
    srand(molecule->m_Max[0] - molecule->m_Min[0]);
    red = rand() / ((double)(RAND_MAX));
    green = rand() / ((double)(RAND_MAX));
    blue = rand() / ((double)(RAND_MAX));
    ApplyColorToGOA(red, green, blue, alpha, molecule);
  } break;
  case PDBParser::GroupOfAtoms::COLLECTION:
    // Here, we still have just one copy. Hence we cannot do much right now.
    // This needs to be done when we render each
    break;
  }
}

void GOAColor::ApplyElementColors(PDBParser::GroupOfAtoms *grp) const {
  unsigned int i;
  for (i = 0; i < grp->m_SubGroups.size(); i++) {
    ApplyElementColors(grp->m_SubGroups.at(i));
  }
  for (i = 0; i < grp->m_Atoms.size(); i++) {
    PDBParser::Atom *at = grp->m_Atoms.at(i);
    float rgb[3];
    float alpha;
    at->getColor(rgb + 0, rgb + 1, rgb + 2, &alpha);
    ApplyColorToAtom(rgb[0], rgb[1], rgb[2], 1.0, at);
  }
}

void GOAColor::ApplyResidueColors(PDBParser::GroupOfAtoms *grp) const {
  unsigned int i;
  for (i = 0; i < grp->m_SubGroups.size(); i++) {
    ApplyResidueColors(grp->m_SubGroups.at(i));
  }
  if (grp->type == PDBParser::RESIDUE_TYPE)
    ApplyColorToGOA(PDBParser::resColorR[grp->resID],
                    PDBParser::resColorG[grp->resID],
                    PDBParser::resColorB[grp->resID], 1.0, grp);
}

void GOAColor::ApplyChainColors(
    PDBParser::GroupOfAtoms *grp,
    ColorManager::ColorGenerator *colorGenerator) const {
  unsigned int i;
  for (i = 0; i < grp->m_SubGroups.size(); i++) {
    ApplyChainColors(grp->m_SubGroups.at(i), colorGenerator);
  }
  if (grp->type == CHAIN_TYPE) {
    ColorManager::Color color = colorGenerator->getNextColor();
    // arand comment, 9-9-2011
    // printf("%f %f %f %f\n", color[0]*256, color[1]*256, color[2]*256,
    // color[3]);
    ApplyColorToGOA(color[0], color[1], color[2], 1.0, grp);
  }
}

bool GOAColor::UsesDefaultColoring(int transIndex) const {
  bool ret = false;
  for (int i = 0; i < m_EntList.size() && !ret; i++)
    if (m_EntList[i].transformSection == transIndex) {
      ret = true; // this will cause the loop to terminate
    }
  return ret;
}

// Private Functions
void GOAColor::ApplyColormapEntryToGOA(const ColormapEntry &ent,
                                       PDBParser::GroupOfAtoms *grp) const {
  // crawl through the ColormapEntry until you hit an entry that's set to
  // its default value. Then apply the entry's color for every atom below
  // whatever level in the GOA you're at.
  // roughly speaking :)
  PDBParser::GroupOfAtoms *atoms = grp;
  if (atoms->type == PDBParser::COLLECTION_TYPE && ent.molecule > -1 &&
      ent.molecule < (int)atoms->m_SubGroups.size()) {
    atoms = atoms->m_SubGroups.at(ent.molecule);
  } else if (atoms->type == PDBParser::COLLECTION_TYPE) {
    // gripe about the out of range index
    printf("Warning: Colormap entry on line %d contains an ", ent.lineNum);
    printf("out-of-range molecule index (%d)\n", ent.molecule);
    printf("Colormap entry skipped\n");
    return;
  }
  if (atoms->type == PDBParser::PROTEIN_TYPE && ent.chainIndex > -1 &&
      ent.chainIndex < (int)atoms->m_SubGroups.size()) {
    atoms = atoms->m_SubGroups.at(ent.chainIndex);
    // check to see if ent.chainID matches the chain id of
    // atoms (the chain id is a member of the Atom class).
    // gripe if it doesn't match.
    PDBParser::Atom *at = getAnAtom(atoms);
    if (at) {
      if (at->m_Chain != ent.chainID) {
        printf("Warning: Colormap entry on line %d contains ", ent.lineNum);
        printf("a chain index/id pair which does not match the parser's ");
        printf("output.\nColormap: %d/%c\n", ent.chainIndex, ent.chainID);
        printf("GroupOfAtoms: %d/%c\n", ent.chainIndex, at->m_Chain);
      }
    }
  } else if (atoms->type == PDBParser::PROTEIN_TYPE && ent.chainIndex != -1) {
    // gripe about the out of range index
    printf("Warning: Colormap entry on line %d contains an ", ent.lineNum);
    printf("out-of-range chain index (%d)\n", ent.chainIndex);
    printf("Colormap entry skipped\n");
    return;
  }
  // IGNORE THE SS INDEX. Search the entire chain for a residue.
  // ColorMapEntry::resIndex
  if (atoms->type == PDBParser::CHAIN_TYPE && ent.resIndex > -1) {
    unsigned int i, j;
    PDBParser::GroupOfAtoms *ss, *found = NULL;
    for (i = 0; i < atoms->m_SubGroups.size() && found == NULL; i++) {
      ss = atoms->m_SubGroups.at(i);
      for (j = 0; j < ss->m_SubGroups.size(); j++) {
        PDBParser::GroupOfAtoms *tmp = ss->m_SubGroups.at(j);
        if (ent.resIndex == tmp->m_ResSeq) {
          found = tmp;
          break;
        }
      }
    }
    // assign what was found (unless nothing's there)
    if (found != NULL) {
      atoms = found;
    }
    // check to see if ent.resName matches the residue name of atoms.
    // gripe if it doesn't match.
    if (atoms->type == PDBParser::RESIDUE_TYPE &&
        strcmp(ent.resName,
               PDBParser::residueTable[atoms->resID].residueName) != 0) {
      printf("Warning: Colormap entry on line %d contains a ", ent.lineNum);
      printf("residue number (%d) which could not be located\n", ent.resIndex);
      printf("Colormap entry skipped\n");
      return;
    }
  } else if ((atoms->type == PDBParser::SS_TYPE ||
              atoms->type == PDBParser::CHAIN_TYPE) &&
             ent.resIndex != -1) {
    // gripe about the out of range index
    printf("Warning: Colormap entry on line %d contains an ", ent.lineNum);
    printf("out-of-range residue index (%d)\n", ent.resIndex);
    printf("Colormap entry skipped\n");
    return;
  }
  if (atoms->type == PDBParser::RESIDUE_TYPE && ent.atomIndex > -1 &&
      ent.atomIndex < (int)atoms->m_Atoms.size()) {
    unsigned int i = 0;
    if (strncmp(ent.atomName, "XXXX", 4) == 0) {
      i = ent.atomIndex;
    } else {
      for (i = 0; i < atoms->m_Atoms.size(); i++) {
        PDBParser::Atom *at = atoms->m_Atoms.at(i);
        if (strncmp(at->getName(), ent.atomName, 4) == 0) {
          break;
        }
      }
    }
    if (i < atoms->m_Atoms.size()) {
      ApplyColorToAtom(ent.red, ent.green, ent.blue, ent.alpha,
                       atoms->m_Atoms.at(i));
    } else {
      printf("Atom could not be found!, line: %d\n", ent.lineNum);
    }
    // we are only coloring one atom, so we're done
    return;
  }

  ApplyColorToGOA(ent.red, ent.green, ent.blue, ent.alpha, atoms);
}

// this old version uses indices derived from the GOA rather than indices
// derived from the PDB file. Turns out that was a tedious way of doing
// things.
void GOAColor::ApplyColormapEntryToGOA_old(const ColormapEntry &ent,
                                           PDBParser::GroupOfAtoms *grp) const {
  // crawl through the ColormapEntry until you hit an entry that's set to
  // its default value. Then apply the entry's color for every atom below
  // whatever level in the GOA you're at.
  // roughly speaking :)
  PDBParser::GroupOfAtoms *atoms = grp;
  if (atoms->type == PDBParser::COLLECTION_TYPE && ent.molecule > -1 &&
      ent.molecule < (int)atoms->m_SubGroups.size()) {
    atoms = atoms->m_SubGroups.at(ent.molecule);
  } else if (atoms->type == PDBParser::COLLECTION_TYPE) {
    // gripe about the out of range index
    printf("Warning: Colormap entry on line %d contains an ", ent.lineNum);
    printf("out-of-range molecule index (%d)\n", ent.molecule);
    printf("Colormap entry skipped\n");
    return;
  }
  if (atoms->type == PDBParser::PROTEIN_TYPE && ent.chainIndex > -1 &&
      ent.chainIndex < (int)atoms->m_SubGroups.size()) {
    atoms = atoms->m_SubGroups.at(ent.chainIndex);
    // check to see if ent.chainID matches the chain id of
    // atoms (the chain id is a member of the Atom class).
    // gripe if it doesn't match.
    PDBParser::Atom *at = getAnAtom(atoms);
    if (at) {
      if (at->m_Chain != ent.chainID) {
        printf("Warning: Colormap entry on line %d contains ", ent.lineNum);
        printf("a chain index/id pair which does not match the parser's ");
        printf("output.\nColormap: %d/%c\n", ent.chainIndex, ent.chainID);
        printf("GroupOfAtoms: %d/%c\n", ent.chainIndex, at->m_Chain);
      }
    }
  } else if (atoms->type == PDBParser::PROTEIN_TYPE && ent.chainIndex != -1) {
    // gripe about the out of range index
    printf("Warning: Colormap entry on line %d contains an ", ent.lineNum);
    printf("out-of-range chain index (%d)\n", ent.chainIndex);
    printf("Colormap entry skipped\n");
    return;
  }
  if (atoms->type == PDBParser::CHAIN_TYPE && ent.ssIndex > -1 &&
      ent.ssIndex < (int)atoms->m_SubGroups.size()) {
    atoms = atoms->m_SubGroups.at(ent.ssIndex);
  } else if (atoms->type == PDBParser::CHAIN_TYPE && ent.ssIndex != -1) {
    // gripe about the out of range index
    printf("Warning: Colormap entry on line %d contains an ", ent.lineNum);
    printf("out-of-range secondary structure index (%d)\n", ent.ssIndex);
    printf("Colormap entry skipped\n");
    return;
  }
  if ((atoms->type == PDBParser::SS_TYPE ||
       atoms->type == PDBParser::CHAIN_TYPE) &&
      ent.resIndex > -1 && ent.resIndex < (int)atoms->m_SubGroups.size()) {
    atoms = atoms->m_SubGroups.at(ent.resIndex);
    // XXX:
    // check to see if ent.resName matches the residue name of atoms.
    // gripe if it doesn't match.
  } else if ((atoms->type == PDBParser::SS_TYPE ||
              atoms->type == PDBParser::CHAIN_TYPE) &&
             ent.resIndex != -1) {
    // gripe about the out of range index
    printf("Warning: Colormap entry on line %d contains an ", ent.lineNum);
    printf("out-of-range residue index (%d)\n", ent.resIndex);
    printf("Colormap entry skipped\n");
    return;
  }
  if (atoms->type == PDBParser::RESIDUE_TYPE && ent.atomIndex > -1 &&
      ent.atomIndex < (int)atoms->m_Atoms.size()) {
    unsigned int i = 0;
    if (strncmp(ent.atomName, "XXXX", 4) == 0) {
      i = ent.atomIndex;
    } else {
      for (i = 0; i < atoms->m_Atoms.size(); i++) {
        PDBParser::Atom *at = atoms->m_Atoms.at(i);
        if (strncmp(at->getName(), ent.atomName, 4) == 0) {
          break;
        }
      }
    }
    ApplyColorToAtom(ent.red, ent.green, ent.blue, ent.alpha,
                     atoms->m_Atoms.at(i));
    // we are only coloring one atom, so we're done
    return;
  }
  ApplyColorToGOA(ent.red, ent.green, ent.blue, ent.alpha, atoms);
}

void GOAColor::ApplyColorToGOA(double red, double green, double blue,
                               double alpha,
                               PDBParser::GroupOfAtoms *grp) const {
  unsigned int i;
  // multiresolution colors
  grp->multi_r = (float)red;
  grp->multi_g = (float)green;
  grp->multi_b = (float)blue;
  for (i = 0; i < grp->m_SubGroups.size(); i++) {
    ApplyColorToGOA(red, green, blue, alpha, grp->m_SubGroups.at(i));
  }
  for (i = 0; i < grp->m_Atoms.size(); i++) {
    ApplyColorToAtom(red, green, blue, alpha, grp->m_Atoms.at(i));
  }
}

void GOAColor::ApplyColorToAtom(float red, float green, float blue, float alpha,
                                PDBParser::Atom *at) const {
  unsigned int r, g, b, a;
  r = (unsigned int)(red * 255.0);
  g = (unsigned int)(green * 255.0);
  b = (unsigned int)(blue * 255.0);
  a = (unsigned int)(alpha * 255.0);
  at->m_UserDefinedColor = ((r << 24) & 0xff000000) | ((g << 16) & 0xff0000) |
                           ((b << 8) & 0xff00) | (a & 0xff);
}

PDBParser::Atom *GOAColor::getAnAtom(PDBParser::GroupOfAtoms *grp) const {
  if (grp->m_Atoms.size() == 0) {
    if (grp->m_SubGroups.size() > 0) {
      return getAnAtom(grp->m_SubGroups.at(0));
    } else // no subgroups ??!?
    {
      return NULL;
    }
  }
  return grp->m_Atoms.at(0);
}

void GOAColor::printEntry(const ColormapEntry &ent) const {
  printf("ColormapEntry:\n");
  printf("\tLine: %d\n", ent.lineNum);
  printf("\tTransformation Index: %d\n", ent.transformSection);
  printf("\tMolecule Index: %d\n", ent.molecule);
  printf("\tChain Index/ID: %d / %c\n", ent.chainIndex, ent.chainID);
  printf("\tSecondary Structure Index: %d\n", ent.ssIndex);
  printf("\tResidue Index/Name: %d / \'%s\'\n", ent.resIndex, ent.resName);
  printf("\tAtom Name/Index: \'%s\' / %d\n", ent.atomName, ent.atomIndex);
  printf("\tColor: (%lf,%lf,%lf,%lf)\n", ent.red, ent.green, ent.blue,
         ent.alpha);
}

char *GOAColor::eatWhitespace(char *buffer) {
  while (*buffer == ' ' || *buffer == '\t' || *buffer == '\n' ||
         *buffer == '\r') {
    buffer++;
  }
  return buffer;
}

char **GOAColor::parseEntryLine(char *line) {
  char *starts[NumTokens], *ends[NumTokens], **tokens, *ptr;
  int i;
  // allocate space for tokens (the return value)
  tokens = (char **)malloc(NumTokens * sizeof(char *));
  // ptr is line-1 because we want line+0 passed to strchr. (see below)
  ptr = line - 1;
  // zero the pointers
  for (i = 0; i < NumTokens; i++) {
    starts[i] = ends[i] = NULL;
  }
  // break up the line
  for (i = 0; i < NumTokens; i++) {
    // find the opening double-quote
    ptr = strchr(ptr + 1, '\"');
    // assign the quote to the starts array
    if (ptr != NULL) {
      starts[i] = ptr;
    } else {
      break; // bail if we run out of double-quotes
    }
    // find the closing double-quote
    ptr = strchr(ptr + 1, '\"');
    // assign the quote to the ends array
    if (ptr != NULL) {
      ends[i] = ptr;
    } else {
      break; // bail if we run out of double-quotes
    }
  }
  // generate a list of tokens from the broken line
  for (i = 0; i < NumTokens; i++) {
    // calculate the length of the token
    int toklen = ends[i] - starts[i];
    // create the token (if it's not NULL or empty)
    if (starts[i] != NULL && ends[i] != NULL && toklen > 1) {
      // allocate space
      tokens[i] = (char *)malloc(toklen);
      // copy only what is between the quotes
      memcpy(tokens[i], starts[i] + 1, ends[i] - starts[i] - 1);
      // NULL terminate
      tokens[i][toklen - 1] = '\0';
    } else {
      tokens[i] = NULL;
    }
  }
  return tokens;
}

int GOAColor::parseSectionMarker(char *line) {
  int ret = -1;
  char *ptr = line, *start = NULL, *end = NULL;
  // find the opening brace
  ptr = strchr(ptr, '[');
  // make sure the opener was found
  if (ptr != NULL) {
    // note it
    start = ptr + 1;
    // then find the closing brace
    ptr = strchr(start, ']');
    if (ptr != NULL) {
      end = ptr;
    }
  }
  if (start != NULL && end != NULL && isdigit(start[0])) {
    ret = atoi(start);
  }
  return ret;
}

bool GOAColor::genEntryFromTokens(char **tokens, ColormapEntry *ent) {
  int i;
  bool ret = true, blankSeen = false;
  for (i = 0; i < NumTokens && ret; i++) {
    switch (i) {
    case MolIndex:
      if (tokens[MolIndex] != NULL && isdigit(tokens[MolIndex][0])) {
        ent->molecule = atoi(tokens[MolIndex]);
      } else {
        printf("Colormap Parser Error:%d:"
               " Missing molecule index\n",
               ent->lineNum);
        ent->molecule = -1;
        ret = false; // molecule must have a value
      }
      break;
    case ChainIndex:
      if (tokens[ChainIndex] != NULL && isdigit(tokens[ChainIndex][0])) {
        ent->chainIndex = atoi(tokens[ChainIndex]);
      } else if (tokens[ChainIndex] != NULL &&
                 !isdigit(tokens[ChainIndex][0])) {
        printf("Colormap Parser Error:%d:"
               " Chain index should be a number.\n",
               ent->lineNum);
        ret = false; // molecule must have a value
      } else if (tokens[ChainIndex] == NULL) {
        blankSeen = true; // note that we've encountered a blank token
        ent->chainIndex = -1;
      }
      break;
    case ChainID:
      if (tokens[ChainID] != NULL &&
          (isalnum(tokens[ChainID][0]) ||
           (tokens[ChainID][0] ==
            ' '))) // SKVINAY bug? some chains come with space as their Id
      {
        // it is an error to have non-blank items follow blanks
        if (blankSeen) {
          printf("Colormap Parser Error:%d:"
                 "Non-blank entry following a blank entry.\n",
                 ent->lineNum);
          //" Non-blank chain ID entry following a blank"
          //" molecule index or chain index.\n");
          ret = false;
        }
        ent->chainID = tokens[ChainID][0];
      } else {
        // make sure that chainIndex was not specified
        // otherwise it is an error
        if (ent->chainIndex != -1) {
          printf("Colormap Parser Error:%d:"
                 " The chain ID entry cannot be blank if a chain index"
                 " is specified.\n",
                 ent->lineNum);
          ret = false;
        }
        blankSeen = true; // note that we've encountered a blank token
        ent->chainID = '\0';
      }
      break;
    case SSIndex:
      if (tokens[SSIndex] != NULL && isdigit(tokens[SSIndex][0])) {
        // it is an error to have non-blank items follow blanks
        if (blankSeen) {
          printf("Colormap Parser Error:%d:"
                 "Non-blank entry following a blank entry.\n",
                 ent->lineNum);
          ret = false;
        }
        ent->ssIndex = atoi(tokens[SSIndex]);
      } else {
        blankSeen = true; // note that we've encountered a blank token
        ent->ssIndex = -1;
      }
      break;
    case ResidueIndex:
      if (tokens[ResidueIndex] != NULL && isdigit(tokens[ResidueIndex][0])) {
        // it is an error to have non-blank items follow blanks
        if (blankSeen) {
          printf("Colormap Parser Error:%d:"
                 "Non-blank entry following a blank entry.\n",
                 ent->lineNum);
          ret = false;
        }
        ent->resIndex = atoi(tokens[ResidueIndex]);
      } else {
        blankSeen = true; // note that we've encountered a blank token
        ent->resIndex = -1;
      }
      break;
    case ResidueName:
      if (tokens[ResidueName] != NULL && isascii(tokens[ResidueName][0])) {
        // it is an error to have non-blank items follow blanks
        if (blankSeen) {
          printf("Colormap Parser Error:%d:"
                 "Non-blank entry following a blank entry.\n",
                 ent->lineNum);
          ret = false;
        }
        strncpy(ent->resName, tokens[ResidueName], 3);
        ent->resName[3] = '\0';
      } else {
        // make sure that resIndex was not specified
        // otherwise it is an error
        if (ent->resIndex != -1) {
          printf("Colormap Parser Error:%d:"
                 " The residue index entry cannot be blank if a residue name"
                 " is specified.\n",
                 ent->lineNum);
          ret = false;
        }
        blankSeen = true; // note that we've encountered a blank token
        memset(ent->resName, '\0', 4);
      }
      break;
    case AtomName:
      if (tokens[AtomName] != NULL && isascii(tokens[AtomName][0])) {
        // it is an error to have non-blank items follow blanks
        if (blankSeen) {
          printf("Colormap Parser Error:%d:"
                 "Non-blank entry following a blank entry.\n",
                 ent->lineNum);
          ret = false;
        }
        strncpy(ent->atomName, tokens[AtomName], 4);
        ent->atomName[4] = '\0';
      } else {
        blankSeen = true; // note that we've encountered a blank token
        memset(ent->atomName, '\0', 5);
      }
      break;
    case AtomIndex:
      if (tokens[AtomIndex] != NULL && isdigit(tokens[AtomIndex][0])) {
        // it is an error to have non-blank items follow blanks
        if (blankSeen) {
          printf("Colormap Parser Error:%d:"
                 "Non-blank entry following a blank entry.\n",
                 ent->lineNum);
          ret = false;
        }
        ent->atomIndex = atoi(tokens[AtomIndex]);
      } else {
        // make sure that atomName was not specified
        // otherwise it is an error
        if (ent->atomName[0] != '\0') {
          printf("Colormap Parser Error:%d:"
                 " The atom index entry cannot be blank if an atom name"
                 " is specified.\n",
                 ent->lineNum);
          ret = false;
        }
        ent->atomIndex = -1;
      }
      break;
    case Red:
      if (tokens[Red] != NULL && isdigit(tokens[Red][0])) {
        sscanf(tokens[Red], "%lf", &(ent->red));
        if (ent->red > 1.0) {
          ent->red = 1.0;
        } else if (ent->red < 0.0) {
          ent->red = 0.0;
        }
      } else {
        printf("Colormap Parser Error:%d:"
               " Each component of a color must be specified.\n",
               ent->lineNum);
        ret = false; // blank color entries are an error
        ent->red = -1.0;
      }
      break;
    case Green:
      if (tokens[Green] != NULL && isdigit(tokens[Green][0])) {
        sscanf(tokens[Green], "%lf", &(ent->green));
        if (ent->green > 1.0) {
          ent->green = 1.0;
        } else if (ent->green < 0.0) {
          ent->green = 0.0;
        }
      } else {
        printf("Colormap Parser Error:%d:"
               " Each component of a color must be specified.\n",
               ent->lineNum);
        ret = false; // blank color entries are an error
        ent->green = -1.0;
      }
      break;
    case Blue:
      if (tokens[Blue] != NULL && isdigit(tokens[Blue][0])) {
        sscanf(tokens[Blue], "%lf", &(ent->blue));
        if (ent->blue > 1.0) {
          ent->blue = 1.0;
        } else if (ent->blue < 0.0) {
          ent->blue = 0.0;
        }
      } else {
        printf("Colormap Parser Error:%d:"
               " Each component of a color must be specified.\n",
               ent->lineNum);
        ret = false; // blank color entries are an error
        ent->blue = -1.0;
      }
      break;
    case Alpha:
      if (tokens[Alpha] != NULL && isdigit(tokens[Alpha][0])) {
        sscanf(tokens[Alpha], "%lf", &(ent->alpha));
        if (ent->alpha > 1.0) {
          ent->alpha = 1.0;
        } else if (ent->alpha < 0.0) {
          ent->alpha = 0.0;
        }
      } else {
        printf("Colormap Parser Error:%d:"
               " Each component of a color must be specified.\n",
               ent->lineNum);
        ret = false; // blank color entries are an error
        ent->alpha = -1.0;
      }
      break;
    default:
      break;
    }
  }
  return ret;
}

void GOAColor::freeTokens(char **tokens) {
  // individual entries
  for (int i = 0; i < NumTokens; i++) {
    if (tokens[i] != NULL) {
      free(tokens[i]);
    }
  }
  // the container
  free(tokens);
}
