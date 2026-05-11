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
// Probably should be member functions of a vector object. . .

#include <AdaptiveIsocontouring/comp.h>
#include <Utility/utility.h>

// given three points, put the unit normal in computed
void computeNormal(const float pt0[3], const float pt1[3], const float pt2[3],
                   float computed[3]) {
  float a[3] = {pt0[0] - pt2[0], pt0[1] - pt2[1], pt0[2] - pt2[2]};
  // find the unit normal;
  // calculate the vector pt2-pt0, denote it a
  float b[3] = {pt1[0] - pt2[0], pt1[1] - pt2[1], pt1[2] - pt2[2]};
  float b_norm = sqrt((b[0] * b[0]) + (b[1] * b[1]) + (b[2] * b[2]));
  float a_norm = sqrt((a[0] * a[0]) + (a[1] * a[1]) + (a[2] * a[2]));
  // compute the unit norm
  for (int j = 0; j < 3; j++) {
    if (a_norm != 0) {
      a[j] /= a_norm;
    }
    if (b_norm != 0) {
      b[j] /= b_norm;
    }
  }
  computed[0] = a[1] * b[2] - a[2] * b[1];
  computed[1] = -(a[0] * b[2] - a[2] * b[0]);
  computed[2] = a[0] * b[1] - a[1] * b[0];
}

void findMidpt(const float pt1[3], const float pt2[3], float mid[3]) {
  for (int j = 0; j < 3; j++) {
    mid[j] = (pt1[j] + pt2[j]) / 2;
  }
}

void copyPoints(const float src[3], float dest[3]) {
  for (int j = 0; j < 3; j++) {
    dest[j] = src[j];
  }
}

float distance2Pts2D(const float pt1[2], const float pt2[2]) {
  return sqrt((pt1[0] - pt2[0]) * (pt1[0] - pt2[0]) +
              (pt1[1] - pt2[1]) * (pt1[1] - pt2[1]));
}

float distance2Pts(const float pt1[3], const float pt2[3]) {
  return sqrt((pt1[0] - pt2[0]) * (pt1[0] - pt2[0]) +
              (pt1[1] - pt2[1]) * (pt1[1] - pt2[1]) +
              (pt1[2] - pt2[2]) * (pt1[2] - pt2[2]));
}

float powerDistance(const float pt1[3], const float pt2[3]) {
  return (pt1[0] - pt2[0]) * (pt1[0] - pt2[0]) +
         (pt1[1] - pt2[1]) * (pt1[1] - pt2[1]) +
         (pt1[2] - pt2[2]) * (pt1[2] - pt2[2]);
}

void printVect(const float *in, int size) {
  cout << "(";
  for (int j = 0; j < size - 1; j++) {
    cout << in[j] << ",";
  }
  cout << in[size - 1] << ")";
}

void printVect(const int *in, int size) {
  cout << "(";
  for (int j = 0; j < size - 1; j++) {
    cout << in[j] << ",";
  }
  cout << in[size - 1] << ")";
}

void printVect(const double *in, int size) {
  cout << "(";
  for (int j = 0; j < size - 1; j++) {
    cout << in[j] << ",";
  }
  cout << in[size - 1] << ")";
}

void copyVect(const float *in, float *out, int size) {
  for (int j = 0; j < size; j++) {
    out[j] = in[j];
  }
}

void copyVect(const int *in, int *out, int size) {
  for (int j = 0; j < size; j++) {
    out[j] = in[j];
  }
}

void copyVect(const double *in, double *out, int size) {
  for (int j = 0; j < size; j++) {
    out[j] = in[j];
  }
}
