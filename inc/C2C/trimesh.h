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
#ifndef TRI_MESH_H
#define TRI_MESH_H

#include <Utility/utility.h>

#define EPS 0.00001

typedef int Triangle[3];
typedef float Position[3];
typedef struct _Position {
  float cord[3];
  int idx;
} VPosition;

typedef struct _EdgeIndex {
  int ix, iy, iz;
  int dir;
  int idx;

  _EdgeIndex(int i, int j, int k, int n) {
    ix = i;
    iy = j;
    iz = k;
    dir = n;
  }

} EdgeIndex;

struct LtEdge {
  bool operator()(const EdgeIndex &e1, const EdgeIndex &e2) const {
    if (e1.iz < e2.iz) {
      return true;
    } else if (e1.iz == e2.iz) {
      if (e1.iy < e2.iy) {
        return true;
      } else if (e1.iy == e2.iy) {
        if (e1.ix < e2.ix) {
          return true;
        } else if (e1.ix == e2.ix) {
          if (e1.dir < e2.dir) {
            return true;
          }
        }
      }
    }
    return false;
  }
};

struct EqPos {
  bool operator()(const VPosition &p1, const VPosition &p2) const {
    return (fabs(p1.cord[0] - p2.cord[0]) < EPS &&
            fabs(p1.cord[1] - p2.cord[1]) < EPS &&
            fabs(p1.cord[2] - p2.cord[2]) < EPS);
  }

  bool operator()(const Position &p1, const Position &p2) const {
    return (fabs(p1[0] - p2[0]) < EPS && fabs(p1[1] - p2[1]) < EPS &&
            fabs(p1[2] - p2[2]) < EPS);
  }
};

struct DiffPos {

  bool operator()(const VPosition &p1, const VPosition &p2) const {
    return (fabs(p1.cord[0] - p2.cord[0]) >= EPS ||
            fabs(p1.cord[1] - p2.cord[1]) >= EPS ||
            fabs(p1.cord[2] - p2.cord[2]) >= EPS);
  }

  bool operator()(const Position &p1, const Position &p2) const {
    return (fabs(p1[0] - p2[0]) >= EPS || fabs(p1[1] - p2[1]) >= EPS ||
            fabs(p1[2] - p2[2]) >= EPS);
  }
};

struct LtPos {
  bool operator()(const VPosition &p1, const VPosition &p2) const {
    if (p1.cord[2] < p2.cord[2] - EPS) {
      return true;
    } else if (fabs(p1.cord[2] - p2.cord[2]) < EPS) {
      if (p1.cord[1] < p2.cord[1] - EPS) {
        return true;
      } else if (fabs(p1.cord[1] - p2.cord[1]) < EPS) {
        if (p1.cord[0] < p2.cord[0] - EPS) {
          return true;
        }
      }
    }
    return false;
  }

  bool operator()(const Position &p1, const Position &p2) const {
    return ((p1[2] < p2[2]) ||
            ((fabs(p1[2] - p2[2]) < EPS) && (p1[1] < p2[1])) ||
            ((fabs(p1[2] - p2[2]) < EPS) && (fabs(p1[1] - p2[1]) < EPS) &&
             (p1[0] < p2[0])));
  }
};

#endif
