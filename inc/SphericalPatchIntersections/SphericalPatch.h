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
#ifndef CCV_SPHERICAL_PATCH_H
#define CCV_SPHERICAL_PATCH_H

#include <Utility/utility.h>

class SphericalPatch {
public:
  SphericalPatch();
  virtual ~SphericalPatch();
  double getClosestDist(double r, double c[3], double point[3], double bmin[3],
                        double bmax[3]);

protected:
  double euclideanDist3D(double p1[3], double p2[3]);
  double intersectCircle(double center[2], double radius, double point[2],
                         double intersection[2]);
  bool intersectionInCube(double intersection[3], double bmin[3],
                          double bmax[3]);
  bool intersectionInCube2D(double intersection[2], double minx, double maxx,
                            double miny, double maxy);
  double pointToCircle(int dim0, int dim1, int dim2, double sphereToCircle,
                       double r, double c[3], double point[3], double bmin[3],
                       double bmax[3], double circleDim0Coordinate);
  double pointToSphere(double r, double c[3], double point[3], double bmin[3],
                       double bmax[3]);
  double pointToLineSegment(double r, double c[3], double point[3],
                            double bmin[3], double bmax[3], double vertex[3],
                            int rayD);
};

#endif
