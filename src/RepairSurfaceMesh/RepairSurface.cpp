/*
  Copyright 2011 The University of Texas at Austin

        Authors: Alex Rand <arand@ices.utexas.edu>
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
#include <RepairSurfaceMesh/RepairMesh.h>
// #include <GeometryFileTypes/GeometryLoader.h>

#include <cstdlib>
#include <cmath>

#include <iostream>
#include <fstream>
#include <map>
#include <vector>

using namespace std;

void usageRepairSurface() {
  cout << endl
       << "Usage: MolSurf -repairSurface <raw input file> <raw output file> "
          "[epsilon]"
       << endl
       << "    Remove duplicate vertices from a mesh." << endl;
}

bool repairSurface(int argc, char **argv) {
  if (argc != 5 && argc != 4) {
    usageRepairSurface();
    return false;
  }

  // Geometry* geometry = GeometryLoader().loadFile(argv[2]);

  double ep = 0.00000000001;
  if (argc == 5)
    ep = atof(argv[4]);

  int nV, nT;
  vector<double> x;
  vector<double> y;
  vector<double> z;
  vector<int> t1;
  vector<int> t2;
  vector<int> t3;

  // read the surface
  ifstream fin(argv[2]);
  fin >> nV >> nT;
  for (int i = 0; i < nV; i++) {
    double xx, yy, zz;
    fin >> xx >> yy >> zz;
    x.push_back(xx);
    y.push_back(yy);
    z.push_back(zz);
  }
  for (int i = 0; i < nT; i++) {
    int tt1, tt2, tt3;
    fin >> tt1 >> tt2 >> tt3;
    t1.push_back(tt1);
    t2.push_back(tt2);
    t3.push_back(tt3);
  }

  // eliminate redundant points
  int nVnew = 0;

  map<int, int> old2new;
  vector<int> redun;
  for (int i = 0; i < nV; i++) {

    bool redundant = false;
    for (int j = 0; j < i && !redundant; j++) {
      // check the point
      // FIXME: not implemented yet...

      double dx, dy, dz;
      dx = x[i] - x[j];
      dy = y[i] - y[j];
      dz = z[i] - z[j];

      double dist = sqrt(dx * dx + dy * dy + dz * dz);
      if (dist < ep) {
        redundant = true;
        redun.push_back(1);
        old2new[i] = old2new[j];
      }
    }

    if (!redundant) {
      old2new[i] = nVnew;
      nVnew++;
      redun.push_back(0);
    }
  }

  // write the surface
  ofstream fout(argv[3]);

  fout << nVnew << " " << nT << endl;
  for (int i = 0; i < nV; i++) {
    if (redun[i] == 0) {
      fout << x[i] << " " << y[i] << " " << z[i] << endl;
    }
  }

  for (int i = 0; i < nT; i++) {
    fout << old2new[t1[i]] << " " << old2new[t2[i]] << " " << old2new[t3[i]]
         << endl;
  }

  return true;
}
