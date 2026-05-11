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
#include <SmoothingFunction/Gaussian.h>

Gaussian::Gaussian(double alpha, int m, int n, int N)
    : SmoothingFunction(alpha, m, n, N) {
  precompute();
}

Gaussian::~Gaussian() {}

bool Gaussian::precompute() {
  if ((m_PhiLength < 1) || (m_PhiHatLength < 1)) {
    return false;
  }
  m_Phi = new double[m_PhiLength];
  m_PhiHat = new double[m_PhiHatLength];
  double b = 2 * alpha * m / (M_PI * (2.0 * alpha - 1));
  {
    int i;
    for (i = 0; i < m_PhiHatLength; i++) {
      m_PhiHat[i] =
          exp(-(M_PI / ((double)n)) * (M_PI / ((double)n)) * b * i * i);
    }
  }
  {
    double b = 2 * alpha * m / (M_PI * (2.0 * alpha - 1));
    double expFactor = -n * n / b;
    double expAmplitude = pow((M_PI * b), (-0.5));
    int i;
    for (i = 0; i < m_PhiLength; i++) {
      m_Phi[i] =
          expAmplitude * exp(i * i * expFactor / ((double)n * n * N * N));
    }
  }
  m_Initialized = true;
  return true;
}
