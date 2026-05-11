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
#include <iostream>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <cstring>
#include <pthread.h>
#ifndef __APPLE__
#include <malloc.h>
#else
#include <malloc/malloc.h>
#endif
#include <time.h>

#ifdef _WIN32
#include <sys/types.h>
#include <sys/timeb.h>
#else
#include <sys/time.h>
#endif

#include <DPG/PG.h>

/**
Class for computing the Hausdorff distance and the average distance between two
point clouds.

We assume that the inputs are meshes in the form of RAW or RAWN file format. The
distance is computed based on the vertices; simplices are ignored.

Given two sets of points X and Y,
1. For every x \in X, we define minDist(x,Y) = min_{y \in Y} (dist(x,y))
2. We define oneWayHausdorffDist(X,Y) = max_{x \in X} minDist(x,Y)
3. We define hausdorffDist(X,Y) = max(oneWayHausdorffDist(X,Y),
oneWayHausdorffDist(Y,X))
4. We define oneWayAverageDist(X,Y) = avg_{x \in X} minDist(x,Y)
4. We define averageDist(X,Y) = avg(oneWayAverageDist(X,Y),
oneWayAverageDist(Y,X))

We use Dynamic Packing Grid (DPG) to speed up search for a nearest neighbor
**/

class HausdorffDistance {
private:
  vector<DPG::Point *> pointSetX;
  vector<DPG::Point *> pointSetY;

  int numPointsX;
  int numPointsY;

  DPG::PG *pgX;
  DPG::PG *pgY;

  double oneWayHausdorffDistXY;
  double oneWayHausdorffDistYX;

  double oneWayAverageDistXY;
  double oneWayAverageDistYX;

  double hausdorffDist;
  double averageDist;

  bool readFile(string fileName1, string fileName2, bool isRAWN);
  bool computeAllDistances();
  bool cleanup();

public:
  HausdorffDistance(string fileName1, string fileName2, bool isRAWN);
  ~HausdorffDistance();

  double getHausdorffDistance() { return hausdorffDist; }
  double getAverageDistance() { return averageDist; }
};
