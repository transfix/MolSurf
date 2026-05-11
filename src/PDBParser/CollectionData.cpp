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
#include <PDBParser/CollectionData.h>

using PDBParser::CollectionData;

CollectionData::CollectionData() {
  m_Classification[0] = '\0';
  m_DepDate[10] = '\0';
  m_PDBID[5] = '\0';
  m_Title = 0;
  m_MoleculeDescription = 0;
  m_Source = 0;
  m_Keywords = 0;
  m_ExperimentType = 0;
  m_Authors = 0;
  m_Citation = 0;
}

CollectionData::CollectionData(const CollectionData &collectionData) {
  strcpy(m_Classification, collectionData.m_Classification);
  strcpy(m_DepDate, collectionData.m_DepDate);
  strcpy(m_PDBID, collectionData.m_PDBID);
  if (collectionData.m_Title) {
    m_Title = new char[strlen(collectionData.m_Title) + 1];
    strcpy(m_Title, collectionData.m_Title);
  } else {
    m_Title = 0;
  }
  if (collectionData.m_MoleculeDescription) {
    m_MoleculeDescription =
        new char[strlen(collectionData.m_MoleculeDescription) + 1];
    strcpy(m_MoleculeDescription, collectionData.m_MoleculeDescription);
  } else {
    m_MoleculeDescription = 0;
  }
  if (collectionData.m_Source) {
    m_Source = new char[strlen(collectionData.m_Source) + 1];
    strcpy(m_Source, collectionData.m_Source);
  } else {
    m_Source = 0;
  }
  if (collectionData.m_Keywords) {
    m_Keywords = new char[strlen(collectionData.m_Keywords) + 1];
    strcpy(m_Keywords, collectionData.m_Keywords);
  } else {
    m_Keywords = 0;
  }
  if (collectionData.m_ExperimentType) {
    m_ExperimentType = new char[strlen(collectionData.m_ExperimentType) + 1];
    strcpy(m_ExperimentType, collectionData.m_ExperimentType);
  } else {
    m_ExperimentType = 0;
  }
  if (collectionData.m_Authors) {
    m_Authors = new char[strlen(collectionData.m_Authors) + 1];
    strcpy(m_Authors, collectionData.m_Authors);
  } else {
    m_Authors = 0;
  }
  if (collectionData.m_Citation) {
    m_Citation = new char[strlen(collectionData.m_Citation) + 1];
    strcpy(m_Citation, collectionData.m_Citation);
  } else {
    m_Citation = 0;
  }
  int i, n;
  n = collectionData.m_BiomoleculeTransformations.size();
  for (i = 0; i < n; i++) {
    BiomoleculeTransformations *subCopy = new BiomoleculeTransformations(
        *(collectionData.m_BiomoleculeTransformations.at(i)));
    m_BiomoleculeTransformations.push_back(subCopy);
  }
  n = collectionData.m_ChainIds.size();
  for (i = 0; i < n; i++) {
    m_ChainIds.push_back(collectionData.m_ChainIds.at(i));
  }
  n = collectionData.m_ResidueIds.size();
  for (i = 0; i < n; i++) {
    m_ResidueIds.push_back(collectionData.m_ResidueIds.at(i));
  }
}

CollectionData::~CollectionData() {
  delete[] m_Title;
  m_Title = 0;
  delete[] m_MoleculeDescription;
  m_MoleculeDescription = 0;
  delete[] m_Source;
  m_Source = 0;
  delete[] m_Keywords;
  m_Keywords = 0;
  delete[] m_ExperimentType;
  m_ExperimentType = 0;
  delete[] m_Authors;
  m_Authors = 0;
  delete[] m_Citation;
  m_Citation = 0;
  // delete the transformations
  {
    std::vector<BiomoleculeTransformations *>::iterator
        iter = m_BiomoleculeTransformations.begin(),
        end = m_BiomoleculeTransformations.end();
    for (; iter != end; ++iter) {
      delete *iter;
    }
    m_BiomoleculeTransformations.clear();
  }
}

// print the data in a nicely formatted manner to fp.
bool CollectionData::print(FILE *fp) {
  if (!fp) {
    return false;
  }
  fprintf(fp, "Begin molecule information\n");
  fprintf(fp, "Classification\n");
  fprintf(fp, "%s\n", m_Classification);
  fprintf(fp, "Date of deposition\n");
  fprintf(fp, "%s\n", m_DepDate);
  fprintf(fp, "PDB ID\n");
  fprintf(fp, "%s\n", m_PDBID);
  if (m_Title) {
    fprintf(fp, "Title\n");
    fprintf(fp, "%s\n", m_Title);
  }
  if (m_MoleculeDescription) {
    fprintf(fp, "Description\n");
    fprintf(fp, "%s\n", m_MoleculeDescription);
  }
  if (m_Source) {
    fprintf(fp, "Source\n");
    fprintf(fp, "%s\n", m_Source);
  }
  if (m_Keywords) {
    fprintf(fp, "Keywords\n");
    fprintf(fp, "%s\n", m_Keywords);
  }
  if (m_ExperimentType) {
    fprintf(fp, "Experiment type\n");
    fprintf(fp, "%s\n", m_ExperimentType);
  }
  if (m_Authors) {
    fprintf(fp, "Authors\n");
    fprintf(fp, "%s\n", m_Authors);
  }
  if (m_Citation) {
    fprintf(fp, "Citation\n");
    fprintf(fp, "%s\n", m_Citation);
  }
  int n;
  n = m_BiomoleculeTransformations.size();
  if (n > 0) {
    fprintf(fp, "Transformation matrices:\n");
    int i;
    for (i = 0; i < n; i++) {
      if (m_BiomoleculeTransformations.at(i)) {
        if (!m_BiomoleculeTransformations.at(i)->print(fp)) {
          fprintf(fp,
                  "Something wrong in printing out transformation number %d\n",
                  i);
        }
      } else {
        fprintf(fp,
                "Something wrong in printing out transformation number %d\n",
                i);
      }
    }
  }
  fprintf(fp, "End molecule information\n");
  return true;
}

/*  Go through each element in the list 'm_BiomoleculeTransformations' */
/*  If any element's list of chainIds has a chainId equal to 'chainId' */
/*  return true, else return false.                                    */
bool CollectionData::isTransformationPresent(char chainId) {
  int i, n;
  n = m_BiomoleculeTransformations.size();
  for (i = 0; i < n; i++) {
    if (m_BiomoleculeTransformations.at(i)->isTransformationPresent(chainId)) {
      return true;
    }
  }
  return false;
}

/*  Go through each element in the list 'm_BiomoleculeTransformations' */
/*  If any element's list of chainIds has a chainId equal to 'chainId' */
/*  return the number of transformations associated with it. If we     */
/*  dont find any such chainId, we return -1, signalling an error.     */
int CollectionData::getNumberOfTransformations(char chainId) {
  int i, n;
  n = m_BiomoleculeTransformations.size();
  for (i = 0; i < n; i++) {
    int numTrans =
        m_BiomoleculeTransformations.at(i)->getNumberOfTransformations(chainId);
    if (numTrans != -1) {
      return numTrans;
    }
  }
  return -1;
}

// a mutation of the above, which returns number of transformations for the
// first set of chains
int CollectionData::getNumberOfTransformations() {
  if (m_BiomoleculeTransformations.size()) {
    return m_BiomoleculeTransformations.at(0)->getNumberOfTransformations();
  }
  return -1;
}

/*  Go through each element in the list 'm_BiomoleculeTransformations' */
/*  If any element's list of chainIds has a chainId equal to 'chainId' */
/*  return the rotations associated with it as a 9xN list of doubles   */
/*  and return true. If we dont find any such chainId, we return false */
bool CollectionData::getRotations(char chainId, double *rotations) {
  if (!rotations) {
    return false;
  }
  int i, n;
  n = m_BiomoleculeTransformations.size();
  for (i = 0; i < n; i++) {
    if (m_BiomoleculeTransformations.at(i)->isTransformationPresent(chainId)) {
      return m_BiomoleculeTransformations.at(i)->getRotations(chainId,
                                                              rotations);
    }
  }
  return false;
}

bool CollectionData::getRotations(double *rotations) {
  if (!rotations) {
    return false;
  }
  if (m_BiomoleculeTransformations.size()) {
    return m_BiomoleculeTransformations.at(0)->getRotations(rotations);
  }
  return false;
}

/*  Go through each element in the list 'm_BiomoleculeTransformations'    */
/*  If any element's list of chainIds has a chainId equal to 'chainId'    */
/*  return the translation associated with it as a 9xN list of doubles    */
/*  and return true. If we dont find any such chainId, we return false    */
bool CollectionData::getTranslations(char chainId, double *translation) {
  if (!translation) {
    return false;
  }
  int i, n;
  n = m_BiomoleculeTransformations.size();
  for (i = 0; i < n; i++) {
    if (m_BiomoleculeTransformations.at(i)->isTransformationPresent(chainId)) {
      return m_BiomoleculeTransformations.at(i)->getTranslations(chainId,
                                                                 translation);
    }
  }
  return false;
}

bool CollectionData::getTranslations(double *translation) {
  if (!translation) {
    return false;
  }
  if (m_BiomoleculeTransformations.size()) {
    return m_BiomoleculeTransformations.at(0)->getTranslations(translation);
  }
  return false;
}
