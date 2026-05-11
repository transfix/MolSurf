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
This class is instantiated with the lists of atoms for the receptor and ligand
molecules. It creates octrees for both of the molecules. Then, given a
transformation it transforms the ligand's (2nd mol) octree, and finds
neighboring leaf-pairs between the two octrees. If a leaf-pair is found, it
comapres pairwise distance between atoms in those leaves to identify the
interface atoms and residues. The identified atoms and residues are used to
compute the residue-residue contact matrix, number of polar and non-polar
residues, the local and global density indexes etc.
*/

#ifndef ATOM_INTERFACE_H
#define ATOM_INTERFACE_H

#include <UsefulMath/Vector.h>
#include <UsefulMath/Matrix.h>
#include <InterfaceStats/InterfaceStats.h>
#include <InterfaceStats/ComputeLeastSquarePlane.h>
#include <PDBParser/Atom.h>
#include <DPG/PG.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stdarg.h>
#include <map>

#if !defined(__APPLE__)
#include <malloc.h>
#endif

#ifdef _WIN32
#include <sys/types.h>
#include <sys/timeb.h>
#else
#include <sys/time.h>
#endif

#ifdef freeMem
#undef freeMem
#endif
#define freeMem(ptr)                                                           \
  {                                                                            \
    if (ptr != NULL)                                                           \
      free(ptr);                                                               \
  }

#ifdef zeroIfLess
#undef zeroIfLess
#endif
#define zeroIfLess(a, b) (((a) < (b)) ? 0 : 1)

using CCVOpenGLMath::Matrix;
using CCVOpenGLMath::Vector;

class AtomsOnInterface {
private:
  typedef struct {
    double cx, cy, cz;
    double cq, cr;
    bool leaf;
    int cPtr[8];
    int atomsStartID, atomsEndID;
  } ATOMS_OCTREE_NODE;

  // Static Atoms, static Octree

  vector<PDBParser::Atom *> staticAtoms;
  bool staticAtomsOctreeBuilt;
  vector<ATOMS_OCTREE_NODE> staticAtomsOctree;
  int staticAtomsOctreeRoot;

  // Moving Atoms, moving Octree

  vector<PDBParser::Atom *> movingAtoms;
  bool movingAtomsOctreeBuilt;
  vector<ATOMS_OCTREE_NODE> movingAtomsOctree;
  int movingAtomsOctreeRoot;

  vector<DPG::Point *> interfaceAtoms;
  ComputeLeastSquarePlane lsp;

  bool priorComputationCleared;

  double minRadius, minRadiusUsed;
  int maxLeafSize, maxLeafSizeUsed;

  double epsilon;
  double interfaceWidth;
  Matrix transMatrix;

  int curNode, maxNode;
  int curSubtreeRoot, maxSubtreeRoot;

  int *movingAtomsSubtreeRoots;
  int movingAtomsSubtreeRootsSize;

  vector<int> pairList;
  map<int, int> movingResidueIDMap;
  map<int, int>
      staticResidueIDMap; // Keeps track of previously found interface residues
                          // by their ids. Also counts the number of its atoms
                          // on the interface. This count is used to classify
                          // the residue as rim/core. If a new residue is found,
                          // then numPolar, .., numRim are updated.

  map<int, int> staticResidueIndexMap; // id to index mapping for residues
  map<int, int> movingResidueIndexMap;

  vector<bool> staticInterfaceAtoms;
  vector<bool> movingInterfaceAtoms;
  vector<int> ligandInterfaceAtomList;

  int numAtoms[2];
  int numPolar[2];
  int numNonPolar[2];
  int numCore[2];
  int numRim[2];

  bool isPolar[20];
  int numResidueAtoms[20];
  map<string, int> residueMap; // Maps residue name to an index
  int resresCont[20][20];

  double localDensityIndex;
  double globalDensityIndex;

  bool printStatus;

  void freeMemory(void);
  void setDefaults(void);
  void initStats(void);
  void copyAtoms(vector<PDBParser::Atom *> &atoms, bool stc);
  int constructAtomsOctree(int atomsStartID, int atomsEndID,
                           vector<PDBParser::Atom *> &atomsT, bool stc);
  bool buildStaticAtomsOctree(void);
  bool buildMovingAtomsOctree(void);
  bool buildOctrees(void);
  void fillMovingAtomsSubtreeRootsArray(int nodeID, int maxNodesInLevel);
  void approximateInteractions(Matrix transMat, int nodeS, int nodeM);
  void updateStats(int staticIndex, int movingIndex);
  void identifyRimCore();

public:
  AtomsOnInterface(vector<PDBParser::Atom *> &molecule1,
                   vector<PDBParser::Atom *> &molecule2);
  ~AtomsOnInterface();

  bool computeInteractions(Matrix transMat, InterfaceStats *intStat);
  vector<int> getAllPairs() { return pairList; }
  vector<int> getLigandInterfaceAtoms();
};

#endif
