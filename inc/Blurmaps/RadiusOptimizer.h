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
#ifndef CVC_TEXMOL_RADIUS_OPTIMIZER_H
#define CVC_TEXMOL_RADIUS_OPTIMIZER_H

#include <Utility/utility.h>

class RadiusOptimizer {
public:
  RadiusOptimizer();
  virtual ~RadiusOptimizer();
  bool init();
  bool getOptimizedRadii(double *centers, double *radii, int numberOfAtoms);

protected:
  bool getExtents(double *centers, double *radii, int numberOfAtoms,
                  double *min, double *max, double *minRadius,
                  double *maxRadius, double *avgRadius);
  bool classifyVolume(double *centers, double *radii, int numberOfAtoms,
                      double *minExt, double *maxExt, double *span,
                      unsigned char *volume, int size);
  void Func_Gradient_Sum_Gauss(double *centers, int numb_atom, double *p,
                               double C, double *func, double *grad,
                               int *classify, int curAtom, double *radii);
  void Curvature_Sum_Gauss(double *atoms, int numb_atom, double *p, double C,
                           int *classify, int mo, double *maxk);
  void Spherical_Fit(double *x, double *y, double *z, int npts, double *radius);
  void gaussinverse(double *a, int n, double eps, int *message);
  void Exchangerowcolumn(double *a, int n, int k, int ik, int jk);
  void MatrixTimesMatrix(double *a, double *v, double *x, int m, int n, int l);
  double DotProduct(double *p1, double *p2);
  int numb_vert;
  double *vert;
  bool m_Initialized;
};

#endif
