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
#include <LBIE_lib/contour3d.h>

#if !defined(__APPLE__)
#include <malloc.h>
#else
#include <stdlib.h>
#endif

#include <stdio.h>

Contour3d::Contour3d() {
  nvert = 0;
  ntri = 0;
  vsize = 500;
  tsize = 1000;
  vert = (float(*)[3])malloc(sizeof(float[3]) * vsize);
  tri = (unsigned int(*)[3])malloc(sizeof(unsigned int[3]) * tsize);
}

Contour3d::~Contour3d() {
  free(vert);
  free(tri);
}

int Contour3d::AddVert(float v_pos[3], float norm[3]) {
  int i;
  if (nvert + 1 > vsize) {
    vsize <<= 1;
    vert = (float(*)[3])realloc(vert, sizeof(float[3]) * vsize);
  }
  for (i = 0; i < 3; i++) {
    vert[nvert][i] = v_pos[i];
  }
  return nvert++;
}

int Contour3d::AddTri(unsigned int v1, unsigned int v2, unsigned int v3) {
  if (ntri + 1 > tsize) {
    tsize <<= 1;
    tri = (unsigned int(*)[3])realloc(tri, sizeof(unsigned int[3]) * tsize);
  }
  tri[ntri][0] = v1;
  tri[ntri][1] = v2;
  tri[ntri][2] = v3;
  return ntri++;
}

void Contour3d::Clear() {
  nvert = 0;
  ntri = 0;
}

void Contour3d::display() {}

int Contour3d::center_vtx(int v1, int v2, int v3) {
  float center_vtx[3], norm[3];
  center_vtx[0] = (vert[v1][0] + vert[v2][0] + vert[v3][0]) / 3.;
  center_vtx[1] = (vert[v1][1] + vert[v2][1] + vert[v3][1]) / 3.;
  center_vtx[2] = (vert[v1][2] + vert[v2][2] + vert[v3][2]) / 3.;
  return AddVert(center_vtx, norm);
}
