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
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <Decimation/point_3.h>

istream &operator>>(istream &input, Point_3 &P) {
  input >> P.xyz[0];
  input >> P.xyz[1];
  input >> P.xyz[2];
  return input;
}

ostream &operator<<(ostream &output, Point_3 P) {
  output << "(" << P.xyz[0] << ", " << P.xyz[1] << ", " << P.xyz[2] << ")";
  return output;
}

ostream &operator<<(ostream &output, Vector_3 P) {
  output << "<" << P.xyz[0] << ", " << P.xyz[1] << ", " << P.xyz[2] << ">";
  return output;
}

Point_3 &Point_3::operator=(const Point_3 &p) {
  xyz[0] = p.xyz[0];
  xyz[1] = p.xyz[1];
  xyz[2] = p.xyz[2];
  return *this;
}

int Point_3::operator==(Point_3 Q) {
  return (xyz[0] == Q.xyz[0] && xyz[1] == Q.xyz[1] && xyz[2] == Q.xyz[2]);
}

int Point_3::operator!=(Point_3 Q) {
  return (xyz[0] != Q.xyz[0] || xyz[1] != Q.xyz[1] || xyz[2] != Q.xyz[2]);
}

Vector_3 Point_3::operator-(Point_3 Q) {
  Vector_3 v;
  v.xyz[0] = xyz[0] - Q.xyz[0];
  v.xyz[1] = xyz[1] - Q.xyz[1];
  v.xyz[2] = xyz[2] - Q.xyz[2];
  return v;
}

Point_3 Point_3::operator+(Vector_3 v) {
  Point_3 P;
  P.xyz[0] = xyz[0] + v.xyz[0];
  P.xyz[1] = xyz[1] + v.xyz[1];
  P.xyz[2] = xyz[2] + v.xyz[2];
  return P;
}

Point_3 Point_3::operator-(Vector_3 v) {
  Point_3 P;
  P.xyz[0] = xyz[0] - v.xyz[0];
  P.xyz[1] = xyz[1] - v.xyz[1];
  P.xyz[2] = xyz[2] - v.xyz[2];
  return P;
}

Point_3 &Point_3::operator+=(Vector_3 v) {
  xyz[0] += v.xyz[0];
  xyz[1] += v.xyz[1];
  xyz[2] += v.xyz[2];
  return *this;
}

Point_3 &Point_3::operator-=(Vector_3 v) {
  xyz[0] -= v.xyz[0];
  xyz[1] -= v.xyz[1];
  xyz[2] -= v.xyz[2];
  return *this;
}

Point_3 operator*(int c, Point_3 Q) {
  Point_3 P;
  P.xyz[0] = c * Q.xyz[0];
  P.xyz[1] = c * Q.xyz[1];
  P.xyz[2] = c * Q.xyz[2];
  return P;
}

Point_3 operator*(double c, Point_3 Q) {
  Point_3 P;
  P.xyz[0] = c * Q.xyz[0];
  P.xyz[1] = c * Q.xyz[1];
  P.xyz[2] = c * Q.xyz[2];
  return P;
}

Point_3 operator*(Point_3 Q, int c) {
  Point_3 P;
  P.xyz[0] = c * Q.xyz[0];
  P.xyz[1] = c * Q.xyz[1];
  P.xyz[2] = c * Q.xyz[2];
  return P;
}

Point_3 operator*(Point_3 Q, double c) {
  Point_3 P;
  P.xyz[0] = c * Q.xyz[0];
  P.xyz[1] = c * Q.xyz[1];
  P.xyz[2] = c * Q.xyz[2];
  return P;
}

Point_3 operator/(Point_3 Q, int c) {
  Point_3 P;
  P.xyz[0] = Q.xyz[0] / c;
  P.xyz[1] = Q.xyz[1] / c;
  P.xyz[2] = Q.xyz[2] / c;
  return P;
}

Point_3 operator/(Point_3 Q, double c) {
  Point_3 P;
  P.xyz[0] = Q.xyz[0] / c;
  P.xyz[1] = Q.xyz[1] / c;
  P.xyz[2] = Q.xyz[2] / c;
  return P;
}

Point_3 operator+(Point_3 Q, Point_3 R) {
  Point_3 P;
  P.xyz[0] = Q.xyz[0] + R.xyz[0];
  P.xyz[1] = Q.xyz[1] + R.xyz[1];
  P.xyz[2] = Q.xyz[2] + R.xyz[2];
  return P;
}

bool operator<(const Point_3 &p1, const Point_3 &p2) {
  if (p1.xyz[0] < p2.xyz[0] && p1.xyz[1] < p2.xyz[1] && p1.xyz[2] < p2.xyz[2]) {
    return true;
  }
  return false;
}

double d(Point_3 P, Point_3 Q) {
  // Euclidean distance
  double dx = P.xyz[0] - Q.xyz[0];
  double dy = P.xyz[1] - Q.xyz[1];
  double dz = P.xyz[2] - Q.xyz[2];
  return sqrt(dx * dx + dy * dy + dz * dz);
}

double d2(Point_3 P, Point_3 Q) {
  // squared distance
  double dx = P.xyz[0] - Q.xyz[0];
  double dy = P.xyz[1] - Q.xyz[1];
  double dz = P.xyz[2] - Q.xyz[2];
  return (dx * dx + dy * dy + dz * dz);
}

void Point_3::push_back(double d) {
  if (_ptr < 3) {
    xyz[_ptr] = d;
    _ptr++;
  }
}

double Point_3::at(int index) const { return xyz[index]; }

double &Point_3::operator[](int index) { return xyz[index]; }

double Point_3::operator[](int index) const { return xyz[index]; }

Vector_3 Vector_3::operator-() {
  Vector_3 v;
  v.xyz[0] = -xyz[0];
  v.xyz[1] = -xyz[1];
  v.xyz[2] = -xyz[2];
  return v;
}

//  Scalar Ops

// Scalar multiplication
Vector_3 operator*(int c, Vector_3 w) {
  Vector_3 v;
  v.xyz[0] = c * w.xyz[0];
  v.xyz[1] = c * w.xyz[1];
  v.xyz[2] = c * w.xyz[2];
  return v;
}

Vector_3 operator*(double c, Vector_3 w) {
  Vector_3 v;
  v.xyz[0] = c * w.xyz[0];
  v.xyz[1] = c * w.xyz[1];
  v.xyz[2] = c * w.xyz[2];
  return v;
}

Vector_3 operator*(Vector_3 w, int c) {
  Vector_3 v;
  v.xyz[0] = c * w.xyz[0];
  v.xyz[1] = c * w.xyz[1];
  v.xyz[2] = c * w.xyz[2];
  return v;
}

Vector_3 operator*(Vector_3 w, double c) {
  Vector_3 v;
  v.xyz[0] = c * w.xyz[0];
  v.xyz[1] = c * w.xyz[1];
  v.xyz[2] = c * w.xyz[2];
  return v;
}

Vector_3 operator/(Vector_3 w, int c) {
  Vector_3 v;
  v.xyz[0] = w.xyz[0] / c;
  v.xyz[1] = w.xyz[1] / c;
  v.xyz[2] = w.xyz[2] / c;
  return v;
}

Vector_3 operator/(Vector_3 w, double c) {
  Vector_3 v;
  v.xyz[0] = w.xyz[0] / c;
  v.xyz[1] = w.xyz[1] / c;
  v.xyz[2] = w.xyz[2] / c;
  return v;
}

Vector_3 Vector_3::operator+(Vector_3 w) {
  Vector_3 v;
  v.xyz[0] = xyz[0] + w.xyz[0];
  v.xyz[1] = xyz[1] + w.xyz[1];
  v.xyz[2] = xyz[2] + w.xyz[2];
  return v;
}

Vector_3 Vector_3::operator-(Vector_3 w) {
  Vector_3 v;
  v.xyz[0] = xyz[0] - w.xyz[0];
  v.xyz[1] = xyz[1] - w.xyz[1];
  v.xyz[2] = xyz[2] - w.xyz[2];
  return v;
}

double Vector_3::operator*(Vector_3 w) {
  return (xyz[0] * w.xyz[0] + xyz[1] * w.xyz[1] + xyz[2] * w.xyz[2]);
}

double Vector_3::operator*(Point_3 w) {
  return (xyz[0] * w.xyz[0] + xyz[1] * w.xyz[1] + xyz[2] * w.xyz[2]);
}

Vector_3 Vector_3::operator^(Vector_3 w) {
  Vector_3 v;
  v.xyz[0] = xyz[1] * w.xyz[2] - xyz[2] * w.xyz[1];
  v.xyz[1] = xyz[2] * w.xyz[0] - xyz[0] * w.xyz[2];
  v.xyz[2] = xyz[0] * w.xyz[1] - xyz[1] * w.xyz[0];
  return v;
}

Vector_3 &Vector_3::operator*=(double c) {
  xyz[0] *= c;
  xyz[1] *= c;
  xyz[2] *= c;
  return *this;
}

Vector_3 &Vector_3::operator/=(double c) {
  xyz[0] /= c;
  xyz[1] /= c;
  xyz[2] /= c;
  return *this;
}

Vector_3 &Vector_3::operator+=(Vector_3 w) {
  xyz[0] += w.xyz[0];
  xyz[1] += w.xyz[1];
  xyz[2] += w.xyz[2];
  return *this;
}

Vector_3 &Vector_3::operator-=(Vector_3 w) {
  xyz[0] -= w.xyz[0];
  xyz[1] -= w.xyz[1];
  xyz[2] -= w.xyz[2];
  return *this;
}

Vector_3 &Vector_3::operator^=(Vector_3 w) {
  double ox = xyz[0], oy = xyz[1], oz = xyz[2];
  xyz[0] = oy * w.xyz[2] - oz * w.xyz[1];
  xyz[1] = oz * w.xyz[0] - ox * w.xyz[2];
  xyz[2] = ox * w.xyz[1] - oy * w.xyz[0];
  return *this;
}

void Vector_3::normalize() {
  double ln = sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
  if (ln == 0) {
    return; // do nothing to zero vector...
  }
  xyz[0] /= ln;
  xyz[1] /= ln;
  xyz[2] /= ln;
}
