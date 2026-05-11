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
#ifndef ARAND_DECIMATE_H
#define ARAND_DECIMATE_H

#include <Utility/utility.h>
#include <Decimation/point_3.h>
#include <LBIE_lib/Geoframe.h>

class arandTri {
public:
  int a;
  int b;
  int c;
};

class ArandDecimation {
public:
  ArandDecimation() {}
  void printRaw(string fname);
  void arandDecimate(Geoframe &geo, int nIt, float percentage, float min_angle,
                     float min_variation);

private:
  // data
  vector<arandTri> triangulation;
  vector<Point_3> vertices;
  vector<Vector_3> normals;
  vector<set<int>> vert2tri;
  vector<int> canDecimate; // -1 removed, 0 cannot decimate, 1 ok to decimate
  vector<bool> triExists;
  // helper functions
  void arandNeighbors(int v, set<int> &neighbors);
  void arandNeighbors(int v, vector<int> &neighbors);
  void arandRetriangulate(int v, vector<arandTri> &newTri);
  // can't move these to class arandTri because tri is just indices, not
  // verticies
  Vector_3 getNormal(arandTri &tri);
  double minAngle(arandTri &tri);
};

#endif
