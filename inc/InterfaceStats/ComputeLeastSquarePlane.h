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
Given a list of interface points, this class first computes the least
squares plane through them.

The problem is defined as,
        minimize \sum_i d_i
        where, d_i = A x_i + B y_i + C z_i + D

Taking partial derivatives we can show that
        D = A x_c + B y_c + C z_c
where, (x_c, y_c, z_c) is the centroid of the points.

If we translate the points so that the centroid lies on the origin, then
the computation can be simplified to show that the solution to the least
squares problem is the eigenvectors of the symmetric matrix M, where-

m00 = \sum_i (x'_i)^2
m11 = \sum_i (y'_i)^2
m22 = \sum_i (z'_i)^2

m01 = \sum_i (x'_i)(y'_i)
m02 = \sum_i (x'_i)(z'_i)
m12 = \sum_i (y'_i)(z'_i)

where, x' = x - x_c and so on.

The eigenvector corresponding to the lowest eigenvalue is the normal to the
least square plane and gives us the values of A, B and C. D can be derived
from these.

Also the eigenvectors corresponding to the other two eigenvalues are orthogonal
to each other and are the major and minor directions of the projected point set
on the least square plane.

The least square plane, and the interface point set is used to compute interface
statistics including the interface area, interface volume, interface width,
gapvolume index etc.

*/

#ifndef LEAST_SQ_PLANE_H
#define LEAST_SQ_PLANE_H

#include <vector>
#include <DPG/PG.h>
#include <InterfaceStats/InterfaceStats.h>

using namespace std;

class ComputeLeastSquarePlane {
private:
  double plane[4];
  double eigenValues[3];
  InterfaceStats *intStat;

  bool computePlane(vector<DPG::Point *> &pts);
  void computeStats(vector<DPG::Point *> &pts);

public:
  ComputeLeastSquarePlane();
  ~ComputeLeastSquarePlane();
  bool compute(vector<DPG::Point *> &, InterfaceStats *is);
};

#endif
