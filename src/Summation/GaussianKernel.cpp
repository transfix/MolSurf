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
#include <Summation/GaussianKernel.h>

using namespace CCVSummationModule;

GaussianKernel::GaussianKernel(double blobbiness, double error) {
  m_Blobbiness = blobbiness;
  m_Error = error;
}

GaussianKernel::~GaussianKernel() {}

bool GaussianKernel::isDecayingKernel() { return true; }

bool GaussianKernel::getLength(double *length, double radius) {
  if (!isDecayingKernel()) {
    return false;
  }
  if (!length || (radius <= 0)) {
    return false;
  }
  double l2 =
      (log(m_Error) / m_Blobbiness + 1) * radius * radius; // in one direction
  (*length) = sqrt(fabs(l2));
  if ((*length) < 1) {
    (*length) =
        1; // what to do ? otherwise may fail on very narrow kernels! SKVINAY
  }
  return true;
}

double GaussianKernel::getFunctionAt(double x, double y, double z, double cx,
                                     double cy, double cz, double radius) {
  double dist2 =
      (x - cx) * (x - cx) + (y - cy) * (y - cy) + (z - cz) * (z - cz);
  double expVal = m_Blobbiness * (dist2 / (radius * radius) - 1);
  return exp(expVal);
}
