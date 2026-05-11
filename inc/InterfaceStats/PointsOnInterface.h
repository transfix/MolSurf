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
This class is instantiated with the lists of quadrature points for the receptor
and ligand molecules. It creates octrees for both of the molecules. Additionally
it inserts all points of each molecule into two separate DPGs (Packing Grid).
Then, given a transformation it transforms the ligand's (2nd mol) octree, and
finds neighboring leaf-pairs between the two octrees. If a leaf-pair is found,
for each atom in a leaf, it verifies whether it has any neighbor in the other
DPG. If a neighbor is found, then the point added to the list of interface
atoms. The list of interface atoms is sent to the ComputeLeastSquarePlane class
for computation of interface statistics.
*/

#ifndef POINT_INTERFACE_H
#define POINT_INTERFACE_H

#include <UsefulMath/Vector.h>
#include <UsefulMath/Matrix.h>
#include <InterfaceStats/InterfaceStats.h>
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

class PointsOnInterface {
private:
  typedef struct {
    double cx, cy, cz;
    double cq, cr;
    bool leaf;
    int cPtr[8];
    int pointsStartID, pointsEndID;
  } POINTS_OCTREE_NODE;

  vector<DPG::WeightedPoint *> pointList1;
  vector<DPG::WeightedPoint *> pointList2;

  vector<bool> staticInterfacePoints;
  vector<bool> movingInterfacePoints;

  vector<DPG::WeightedPoint *> interfacePointList;

  bool staticOctreeBuilt;
  bool movingOctreeBuilt;

  vector<POINTS_OCTREE_NODE> staticOctree;
  vector<POINTS_OCTREE_NODE> movingOctree;

  int staticOctreeRoot;
  int movingOctreeRoot;

  double minRadius, minRadiusUsed;
  int maxLeafSize, maxLeafSizeUsed;

  double epsilon;
  double interfaceWidth;
  double area;
  Matrix transMatrix;

  int curNode, maxNode;
  int curSubtreeRoot, maxSubtreeRoot;

  DPG::PG *pgPoints1;
  DPG::PG *pgPoints2;
  int xlatePG1, xlatePG2;

  bool printStatus;

  void freeMemory(void);
  void setDefaults(void);
  void copyPoints(vector<DPG::WeightedPoint *> &pts, bool firstSet);

  void computeXlateForPG();
  void buildPGs();

  int constructOctree(int atomsStartID, int atomsEndID,
                      vector<DPG::WeightedPoint *> &atomsT, bool firstSet);
  bool buildStaticOctree(void);
  bool buildMovingOctree(void);
  bool buildOctrees(void);

  void approximateInteractions(Matrix transMat, Matrix transMatInv, int nodeS,
                               int nodeM);

  Matrix computeInverse(Matrix transMat);

public:
  PointsOnInterface(vector<DPG::WeightedPoint *> &molecule1,
                    vector<DPG::WeightedPoint *> &molecule2);
  ~PointsOnInterface();

  bool computeInteractions(Matrix transMat, InterfaceStats *intStat);
};

#endif
