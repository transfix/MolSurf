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
#include <AdaptiveIsocontouring/Cell.h>

using namespace AdaptiveIsocontouring;

// default constructor, this DOES NOT supply the sufficient data for a cell
Cell::Cell() {}

Cell::Cell(float data[8][3], Cell *parent_, int ci, int level_) {
  level = level_;
  parent = parent_;
  hasContour = false;
  hasChildren = false;
  child_index = ci;
  mesh_vert_index = -1;
  for (int r = 0; r < 8; r++) {
    for (int c = 0; c < 3; c++) {
      corners[r][c] = data[r][c];
    }
    // taking care of setting func_vals to NULL in here as well
    func_vals[r] = NULL;
    children[r] = NULL;
  }
  for (int j = 0; j < 12; j++) {
    edges_tagged[j] = false;
  }
}

Cell::~Cell() {
  for (int j = 0; j < 8; j++) {
    if (children[j] != NULL) {
      delete children[j];
    }
    if (func_vals[j] != NULL) {
      delete func_vals[j];
    }
  }
}

void Cell::setParticlesForCell(Cell *c, vector<float *> p) {
  for (int pc = 0; pc < p.size(); pc++) {
    if ((p[pc][0] >= c->min_X() && p[pc][0] <= c->max_X() &&
         p[pc][1] >= c->min_Y() && p[pc][1] <= c->max_Y() &&
         p[pc][2] >= c->min_Z() && p[pc][2] <= c->max_Z())) {
      c->particles.push_back(&p[pc][0]);
    }
  }
}

// finds the min max of the points
void Cell::findMinMax(float *min, float *max, vector<float *> &in) {
  if (in.size() > 0) {
    for (int j = 0; j < 3; j++) {
      min[j] = max[j] = in[0][j];
    }
  } else {
    return;
  }
  for (int j = 0; j < in.size(); j++) {
    if (in[j][0] > max[0]) {
      max[0] = in[j][0];
    }
    if (in[j][0] < min[0]) {
      min[0] = in[j][0];
    }
    if (in[j][1] > max[1]) {
      max[1] = in[j][1];
    }
    if (in[j][1] < min[1]) {
      min[1] = in[j][1];
    }
    if (in[j][2] > max[2]) {
      max[2] = in[j][2];
    }
    if (in[j][2] < min[2]) {
      min[2] = in[j][2];
    }
  }
}

void Cell::findMinMax(float *min, float *max, vector<double *> &in) {
  if (in.size() > 0) {
    for (int j = 0; j < 3; j++) {
      min[j] = max[j] = (float)in[0][j];
    }
  } else {
    return;
  }
  for (int j = 0; j < in.size(); j++) {
    if (in[j][0] > max[0]) {
      max[0] = (float)in[j][0];
    }
    if (in[j][0] < min[0]) {
      min[0] = (float)in[j][0];
    }
    if (in[j][1] > max[1]) {
      max[1] = (float)in[j][1];
    }
    if (in[j][1] < min[1]) {
      min[1] = (float)in[j][1];
    }
    if (in[j][2] > max[2]) {
      max[2] = (float)in[j][2];
    }
    if (in[j][2] < min[2]) {
      min[2] = (float)in[j][2];
    }
  }
}

float Cell::max_X() { return corners[7][0]; }

float Cell::min_X() { return corners[0][0]; }

float Cell::max_Y() { return corners[7][1]; }

float Cell::min_Y() { return corners[0][1]; }

float Cell::max_Z() { return corners[7][2]; }

float Cell::min_Z() { return corners[0][2]; }
