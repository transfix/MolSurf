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
    AFX_COLLECTIONDATA_H__622FD44A_FAFF_4126_82A2_CB5D720AAA6F__INCLUDED_)
#define AFX_COLLECTIONDATA_H__622FD44A_FAFF_4126_82A2_CB5D720AAA6F__INCLUDED_

#include <Utility/utility.h>
#include <PDBParser/BiomoleculeTransformations.h>

namespace PDBParser {
class BiomoleculeTransformations;
// This class contains information about the molecule including the remarks in
// the pdb file etc.
class CollectionData {
public:
  CollectionData();
  CollectionData(const CollectionData &collectionData);
  virtual ~CollectionData();
  bool print(FILE *fp);
  bool isTransformationPresent(char chainId);
  int getNumberOfTransformations(char chainId);
  bool getRotations(char chainId, double *rotation);
  bool getTranslations(char chainId, double *translation);
  int getNumberOfTransformations();
  bool getRotations(double *rotation);
  bool getTranslations(double *translation);
  // HEADER
  char m_Classification[41];
  char m_DepDate[10];
  char m_PDBID[5];
  // TITLE
  char *m_Title; // can span multiple lines.
  // COMPND
  char *m_MoleculeDescription; // can span multiple lines.
  // SOURCE
  char *m_Source; // can span multiple lines.
  // KEYWDS
  char *m_Keywords; // can span multiple lines.
  // EXPDTA
  char *m_ExperimentType; // can span multiple lines.
  // AUTHOR
  char *m_Authors; // can span multiple lines.
  // JRNL
  char *m_Citation; // can span multiple lines.
  // something for remarks. Should keep symmetry information here.
  // We need a set of <chains, transformations>
  std::vector<BiomoleculeTransformations *> m_BiomoleculeTransformations;
  // this doesnt really belong here, but this is as good as any other place!
  // chain ids is a pair, of chain id and the number of residues in it.
  // residue ids is a triplet < residue name , residue sequence number, residue
  // insertion code >
  std::vector<int> m_ChainIds;
  std::vector<int> m_ResidueIds;
};
}; // namespace PDBParser

#endif
