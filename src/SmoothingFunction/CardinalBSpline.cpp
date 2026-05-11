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
#include <SmoothingFunction/CardinalBSpline.h>

CardinalBSpline::CardinalBSpline(double alpha, int m, int n, int N)
    : SmoothingFunction(alpha, m, n, N) {
  precompute();
}

CardinalBSpline::~CardinalBSpline() {}

bool CardinalBSpline::precompute() {
  if ((m_PhiLength < 1) || (m_PhiHatLength < 1)) {
    return false;
  }
  m_Phi = new double[m_PhiLength];
  m_PhiHat = new double[m_PhiHatLength];
  {
    int i;
    double dn = (double)n;
    m_PhiHat[0] = 1 / dn;
    for (i = 1; i < m_PhiHatLength; i++) {
      m_PhiHat[i] = pow(sin(i * M_PI / dn) / (i * M_PI / dn), 2 * m) / dn;
    }
  }
  {
    int i;
    double dn = (double)n;
    for (i = 0; i < m_PhiLength; i++) {
      double pos = m + i / ((double)(N));
      m_Phi[i] = spline(2 * m, pos) / dn;
    }
  }
  m_Initialized = true;
  return true;
}

void CardinalBSpline::iterate(int k, double x, double *A, int j, int low_g,
                              int high_g, int r) {
  int index;
  int i;
  double a;
  i = high_g + r - k + 1;
  index = high_g;
  while (index >= low_g) {
    a = (x - i) / (k - j);
    A[index] = (1 - a) * A[index - 1] + a * A[index];
    index--;
    i--;
  }
}

double CardinalBSpline::spline(int k, double x) {
  double value; /* Feld fuer Rueckgabewert             */
  int r;        /* Welcher der Splines vom Grad 1 ist  */
  /* ungleich 0?                         */
  int g1, g2;           /* Laufgrenzen zu Vermeidung der Null  */
  int j, index, ug, og; /* Hilfsindices                        */
  double a;             /* Var. alpha des de Boor-Algorithmus  */
  double *A = new double[2 * m];
  /* Determine r */
  if (x <= 0 || x >= k) /* x in range ?*/
  {
    value = 0; /* No : value 0*/
  } else       /* Yes: */
  {
    r = 0; /* Calculate r*/
    while (r < x) {
      r++;
    }
    r--;
    /* Initialization */
    for (index = 0; index < k; index++) {
      A[index] = 0.0;
    }
    A[k - r - 1] = 1.0;
    g1 = r;
    g2 = k - 1 - r;
    ug = g2;
    og = g2;
    /* Calculate the values */
    if (g1 > g2) /* de - Boor - algorithm imple-      */
    {
      /* mentation, for the evaluation     */
      j = 1;          /* of einzelnen B-Splines only the   */
      while (j <= g2) /* not-Null-coefficient in de-       */
      {
        /* Boor- Schema considered           */
        og++;
        a = (x - r + k - 1 - og) / (k - j);
        A[og] = (1 - a) * A[og - 1];
        iterate(k, x, A, j, ug + 1, og - 1, r);
        a = (x - r + k - 1 - ug) / (k - j);
        A[ug] = a * A[ug];
        j++;
      }
      while (j <= g1) {
        ug++;
        og++;
        a = (x - r + k - 1 - og) / (k - j);
        A[og] = (1 - a) * A[og - 1];
        iterate(k, x, A, j, ug, og - 1, r);
        j++;
      }
      while (j < k) {
        ug++;
        iterate(k, x, A, j, ug, og, r);
        j++;
      }
    } else {
      j = 1;
      while (j <= g1) {
        og++;
        a = (x - r + k - 1 - og) / (k - j);
        A[og] = (1 - a) * A[og - 1];
        iterate(k, x, A, j, ug + 1, og - 1, r);
        a = (x - r + k - 1 - ug) / (k - j);
        A[ug] = a * A[ug];
        j++;
      }
      while (j <= g2) {
        iterate(k, x, A, j, ug + 1, og, r);
        a = (x - r + k - 1 - ug) / (k - j);
        A[ug] = a * A[ug];
        j++;
      }
      while (j < k) {
        ug++;
        iterate(k, x, A, j, ug, og, r);
        j++;
      }
    }
    value = A[k - 1];
  }
  delete[] A;
  return (value);
}
