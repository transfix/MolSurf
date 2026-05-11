/*
  Copyright 2011 The University of Texas at Austin

	Advisor: Chandrajit
 * Bajaj <bajaj@cs.utexas.edu>

  This file is part of MolSurf.

  MolSurf is
 * free software; you can redistribute it and/or
  modify it under the terms of
 * the GNU Lesser General Public
  License version 2.1 as published by the Free
 * Software Foundation.


  MolSurf is distributed in the hope that it will be
 * useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of

 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General
 * Public License for more details.

  You should have received a copy of the
 * GNU General Public License
  along with MolSurf; if not, write to the Free
 * Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
*/
/*
Given the pqr and quad files of receptor and ligand,
and a list of docking transformations,
this class computes all relevant interface statistics
for each of the transformations.
*/

#ifndef COMP_INT_STATS_H
#define COMP_INT_STATS_H

#include <InterfaceStats/InterfaceStats.h>
#include <InterfaceStats/AtomsOnInterface.h>
#include <InterfaceStats/PointsOnInterface.h>

#include <Utility/utility.h>
#include <PDBParser/Atom.h>
#include <PDBParser/CollectionData.h>
#include <PDBParser/ElementInformation.h>
#include <PDBParser/GroupOfAtoms.h>
#include <PDBParser/ResidueInformation.h>
#include <UsefulMath/Vector.h>
#include <UsefulMath/Matrix.h>
#include <DPG/PG.h>

using namespace std;
using CCVOpenGLMath::Matrix;
using CCVOpenGLMath::Vector;

namespace PDBParser {
class GroupOfAtoms;
class Atom;
}; // namespace PDBParser

class ComputeInterfaceStats {
private:
  vector<PDBParser::Atom *>
      atomList1; // list of atoms of molecule 1, read from pqr
  vector<PDBParser::Atom *>
      atomList2; // list of atoms of molecule 2, read from pqr

  vector<DPG::WeightedPoint *>
      qPoints1; // list of quad points of molecule 1, read from quadfile
  vector<DPG::WeightedPoint *>
      qPoints2; // list of quad points of molecule 1, read from quadfile

  vector<Matrix> transformations; // list of transformations from docking output
  vector<InterfaceStats *>
      interfaces; // interface stats for each transformation
  int numInterfaces;

  vector<PDBParser::Atom *>
      interfaceAtomList1; // temporary lists used during the computation
                          // interface stats for a single transformation
  vector<PDBParser::Atom *> interfaceAtomList2;
  vector<DPG::Point *> interfacePoints;

  AtomsOnInterface *atomInterface;
  PointsOnInterface *pointInterface;

  void readAtomsFromFiles(string fileName1, string fileName2);
  void readQuadPointsFromFiles(string fileName1, string fileName2);

  void computeAllStats();

public:
  ComputeInterfaceStats(string pqrFileName1, string pqrFileName2,
                        string quadFileName1, string quadFileName2,
                        vector<Matrix> xformations);
  ~ComputeInterfaceStats();

  InterfaceStats *getInterfaceStats(int index) {
    if (index < numInterfaces)
      return interfaces[index];
    return NULL;
  }

  double getInterfaceWidth(int index) {
    if (index < numInterfaces)
      return interfaces[index]->getInterfaceWidth();
    return 0;
  }
  double getPlanarityIndex(int index) {
    if (index < numInterfaces)
      return interfaces[index]->getPlanarityIndex();
    return 0;
  }
  double getCircularityIndex(int index) {
    if (index < numInterfaces)
      return interfaces[index]->getCircularityIndex();
    return 0;
  }
  double getInterfaceArea(int index) {
    if (index < numInterfaces)
      return interfaces[index]->getInterfaceArea();
    return 0;
  }
  double getInterfaceVolume(int index) {
    if (index < numInterfaces)
      return interfaces[index]->getInterfaceVolume();
    return 0;
  }
  double getGapIndex(int index) {
    if (index < numInterfaces)
      return interfaces[index]->getGapIndex();
    return 0;
  }
  double getLocalDensityIndex(int index) {
    if (index < numInterfaces)
      return interfaces[index]->getLocalDensityIndex();
    return 0;
  }
  double getGlobalDensityIndex(int index) {
    if (index < numInterfaces)
      return interfaces[index]->getGlobalDensityIndex();
    return 0;
  }
  int getResresContact(int index, int i, int j) {
    if (index < numInterfaces)
      return interfaces[index]->getResresContact(i, j);
    return -1;
  }
  int getNumCoreResidue(int index, int stc) {
    if (index < numInterfaces)
      return interfaces[index]->getNumCoreResidue(stc);
    return -1;
  }
  int getNumRimResidue(int index, int stc) {
    if (index < numInterfaces)
      return interfaces[index]->getNumRimResidue(stc);
    return -1;
  }
  int getNumPolarResidue(int index, int stc) {
    if (index < numInterfaces)
      return interfaces[index]->getNumPolarResidue(stc);
    return -1;
  }
  int getNumNonPolarResidue(int index, int stc) {
    if (index < numInterfaces)
      return interfaces[index]->getNumNonPolarResidue(stc);
    return -1;
  }
};

#endif
