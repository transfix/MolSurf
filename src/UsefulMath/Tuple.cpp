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
#include <UsefulMath/Tuple.h>

using CCVOpenGLMath::Tuple;

Tuple::Tuple(float x, float y, float z, float w) { set(x, y, z, w); }

Tuple::Tuple() { set(0.0, 0.0, 0.0, 0.0); }

Tuple::~Tuple() {}

Tuple::Tuple(const Tuple &copy) { set(copy); }

Tuple &Tuple::operator=(const Tuple &copy) { return set(copy); }

Tuple &Tuple::set(float x, float y, float z, float w) {
  p[0] = x;
  p[1] = y;
  p[2] = z;
  p[3] = w;
  return *this;
}

Tuple &Tuple::set(float *array) {
  p[0] = array[0];
  p[1] = array[1];
  p[2] = array[2];
  p[3] = array[3];
  return *this;
}

Tuple &Tuple::set(const Tuple &copy) {
  if (this != &copy) {
    p[0] = copy.p[0];
    p[1] = copy.p[1];
    p[2] = copy.p[2];
    p[3] = copy.p[3];
  }
  return *this;
}

float &Tuple::operator[](unsigned int i) { return p[i]; }

const float &Tuple::operator[](unsigned int i) const { return p[i]; }

void Tuple::print() const {
  printf("[%5.3lf %5.3lf %5.3lf %5.3lf ]\n", p[0], p[1], p[2], p[3]);
}
