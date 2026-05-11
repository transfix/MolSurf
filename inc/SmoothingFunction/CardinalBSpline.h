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
#if !defined(                                                                  \
    AFX_CARDINALBSPLINE_H__750E6941_C414_46FE_BCCD_A668221A2CA8__INCLUDED_)
#define AFX_CARDINALBSPLINE_H__750E6941_C414_46FE_BCCD_A668221A2CA8__INCLUDED_

#include <Utility/utility.h>
#include <SmoothingFunction/SmoothingFunction.h>

class CardinalBSpline : public SmoothingFunction {
public:
  CardinalBSpline(double alpha, int m, int n, int N);
  virtual ~CardinalBSpline();
  inline double getPhiAtRealPos(double pos);
  inline double getPhiHatAtRealPos(double pos);

protected:
  bool precompute();
  void iterate(int k, double x, double *A, int j, int low_g, int high_g, int r);
  double spline(int k, double x);
};

double CardinalBSpline::getPhiAtRealPos(double pos) { return 0; }

double CardinalBSpline::getPhiHatAtRealPos(double pos) { return 0; }

#endif
