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
#include <ASpline/quadrature.h>

Quadrature::Quadrature() { numOfPts_perTri = 0; }

Quadrature::Quadrature(const char *type, int numOfPts) {
  numOfPts_perTri = numOfPts;
  xi1 = new double[numOfPts_perTri];
  xi2 = new double[numOfPts_perTri];
  xi3 = new double[numOfPts_perTri];
  w = new double[numOfPts_perTri];
  if (strcmp(type, "Gaussian") == 0 || strcmp(type, "gaussian") == 0) {
    if (numOfPts == 1) {
      onePtGaussian();
    } else if (numOfPts == 3) {
      threePtGaussian();
    } else if (numOfPts == 6) {
      sixPtGaussian();
    } else if (numOfPts == 7) {
      sevenPtGaussian();
    } else if (numOfPts == 12) {
      twelvePtGaussian();
    } else {
      printf("Gaussian rule: 1, 3, 6, 7, or 12 pt\n");
      exit(0);
    }
  } else {
    printf("(1, 3, 6, 7, or 12 pt) Guassian or (...pt) ...\n");
    exit(0);
  }
}

// Gaussian
void Quadrature::onePtGaussian() {
  xi1[0] = 1.0 / 3.0;
  xi2[0] = 1.0 / 3.0;
  xi3[0] = 1.0 / 3.0;
  w[0] = 1.0;
}

void Quadrature::threePtGaussian() {
  xi1[0] = 2.0 / 3.0;
  xi2[0] = 1.0 / 6.0;
  xi3[0] = 1.0 / 6.0;
  xi1[1] = 1.0 / 6.0;
  xi2[1] = 2.0 / 3.0;
  xi3[1] = 1.0 / 6.0;
  xi1[2] = 1.0 / 6.0;
  xi2[2] = 1.0 / 6.0;
  xi3[2] = 2.0 / 3.0;
  w[0] = 1.0 / 3.0;
  w[1] = w[0];
  w[2] = w[0];
}

void Quadrature::sixPtGaussian() {
  xi1[0] = 0.81684757;
  xi2[0] = 0.09157621;
  xi3[0] = xi2[0];
  xi1[1] = xi2[0];
  xi2[1] = xi1[0];
  xi3[1] = xi3[0];
  xi1[2] = xi2[0];
  xi2[2] = xi2[0];
  xi3[2] = xi1[0];
  xi1[3] = 0.10810302;
  xi2[3] = 0.44594849;
  xi3[3] = xi2[3];
  xi1[4] = xi2[3];
  xi2[4] = xi1[3];
  xi3[4] = xi3[3];
  xi1[5] = xi2[3];
  xi2[5] = xi2[3];
  xi3[5] = xi1[3];
  w[0] = 0.10995174;
  w[1] = w[0];
  w[2] = w[0];
  w[3] = 0.22338159;
  w[4] = w[3];
  w[5] = w[3];
}

void Quadrature::sevenPtGaussian() {
  xi1[0] = 1.0 / 3.0;
  xi2[0] = 1.0 / 3.0;
  xi3[0] = 1.0 / 3.0;
  xi1[1] = 0.79742699;
  xi2[1] = 0.10128651;
  xi3[1] = xi2[1];
  xi1[2] = xi2[1];
  xi2[2] = xi1[1];
  xi3[2] = xi3[1];
  xi1[3] = xi2[1];
  xi2[3] = xi2[1];
  xi3[3] = xi1[1];
  xi1[4] = 0.05971587;
  xi2[4] = 0.47014206;
  xi3[4] = xi2[4];
  xi1[5] = xi2[4];
  xi2[5] = xi1[4];
  xi3[5] = xi3[4];
  xi1[6] = xi2[4];
  xi2[6] = xi2[4];
  xi3[6] = xi1[4];
  w[0] = 0.225;
  w[1] = 0.12593918;
  w[2] = w[1];
  w[3] = w[1];
  w[4] = 0.13239415;
  w[5] = w[4];
  w[6] = w[4];
}

void Quadrature::twelvePtGaussian() {
  int i;
  xi1[0] = 0.87382197;
  xi2[0] = 0.06308901;
  xi3[0] = xi2[0];
  xi1[1] = xi2[0];
  xi2[1] = xi1[0];
  xi3[1] = xi3[0];
  xi1[2] = xi2[0];
  xi2[2] = xi2[0];
  xi3[2] = xi1[0];
  xi1[3] = 0.50142651;
  xi2[3] = 0.24928675;
  xi3[3] = xi2[3];
  xi1[4] = xi2[3];
  xi2[4] = xi1[3];
  xi3[4] = xi3[3];
  xi1[5] = xi2[3];
  xi2[5] = xi2[3];
  xi3[5] = xi1[3];
  xi1[6] = 0.63650250;
  xi2[6] = 0.31035245;
  xi3[6] = 0.05314505;
  xi1[7] = xi3[6];
  xi2[7] = xi1[6];
  xi3[7] = xi2[6];
  xi1[8] = xi3[7];
  xi2[8] = xi1[7];
  xi3[8] = xi2[7];
  xi1[9] = xi1[6];
  xi2[9] = xi3[6];
  xi3[9] = xi2[6];
  xi1[10] = xi3[9];
  xi2[10] = xi1[9];
  xi3[10] = xi2[9];
  xi1[11] = xi3[10];
  xi2[11] = xi1[10];
  xi3[11] = xi2[10];
  w[0] = 0.05084491;
  w[1] = w[0];
  w[2] = w[0];
  w[3] = 0.11678628;
  w[4] = w[3];
  w[5] = w[3];
  for (i = 6; i < 12; i++) {
    w[i] = 0.08285108;
  }
}
