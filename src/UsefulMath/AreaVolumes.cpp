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
#include <UsefulMath/AreaVolumes.h>

using namespace CCVOpenGLMath;

AreaVolumes::AreaVolumes() {}

AreaVolumes::~AreaVolumes() {}

double AreaVolumes::getTriangleArea(double v1x, double v1y, double v1z,
                                    double v2x, double v2y, double v2z,
                                    double v3x, double v3y, double v3z) {
  // the area of a triangle with lengths a,b,c is given as
  // 1/4 sqrt( (a+b+c) (b+c-a) (c+a-b) (a+b-c) )
  double a = sqrt((v1x - v2x) * (v1x - v2x) + (v1y - v2y) * (v1y - v2y) +
                  (v1z - v2z) * (v1z - v2z));
  double b = sqrt((v1x - v3x) * (v1x - v3x) + (v1y - v3y) * (v1y - v3y) +
                  (v1z - v3z) * (v1z - v3z));
  double c = sqrt((v3x - v2x) * (v3x - v2x) + (v3y - v2y) * (v3y - v2y) +
                  (v3z - v2z) * (v3z - v2z));
  return 1.0 / 4.0 *
         sqrt(fabs((a + b + c) * (b + c - a) * (c + a - b) * (a + b - c)));
}
