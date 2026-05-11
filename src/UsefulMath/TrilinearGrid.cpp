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
#include <UsefulMath/TrilinearGrid.h>

using CCVOpenGLMath::TrilinearGrid;

TrilinearGrid::TrilinearGrid() {}

TrilinearGrid::~TrilinearGrid() {}

void TrilinearGrid::cell2xyz(int cell_index, int &x, int &y, int &z,
                             unsigned int *gdim) {
  uint temp0 = gdim[0] - 1;
  uint temp1 = gdim[1] - 1;

  x = cell_index % temp0;
  y = (cell_index / temp0) % temp1;
  z = cell_index / (temp0 * temp1);
}

int TrilinearGrid::xyz2cell(int x, int y, int z, unsigned int *gdim) {
  return x + y * (gdim[0] - 1) + z * (gdim[0] - 1) * (gdim[1] - 1);
}

void TrilinearGrid::vtx2xyz(int vtx_idx, int &x, int &y, int &z,
                            unsigned int *gdim) {
  x = vtx_idx % (gdim[0]);
  y = (vtx_idx / (gdim[0])) % (gdim[1]);
  z = vtx_idx / ((gdim[0]) * (gdim[1]));
}

void TrilinearGrid::getCellVertices(int cell_index, int *cellVertexArray,
                                    unsigned int *gdim) {
  int x, y, z;
  cell2xyz(cell_index, x, y, z, gdim);
  cellVertexArray[0] = TrilinearGrid::xyz2vtx(x, y, z, gdim);
  cellVertexArray[1] = TrilinearGrid::xyz2vtx(x + 1, y, z, gdim);
  cellVertexArray[2] = TrilinearGrid::xyz2vtx(x, y + 1, z, gdim);
  cellVertexArray[3] = TrilinearGrid::xyz2vtx(x + 1, y + 1, z, gdim);
  cellVertexArray[4] = TrilinearGrid::xyz2vtx(x, y, z + 1, gdim);
  cellVertexArray[5] = TrilinearGrid::xyz2vtx(x + 1, y, z + 1, gdim);
  cellVertexArray[6] = TrilinearGrid::xyz2vtx(x, y + 1, z + 1, gdim);
  cellVertexArray[7] = TrilinearGrid::xyz2vtx(x + 1, y + 1, z + 1, gdim);
}

int TrilinearGrid::xyz2vtx(int x, int y, int z, unsigned int *dim) {
  return x + y * dim[0] + z * dim[0] * dim[1];
}

int TrilinearGrid::getNeighbor(int i, int j, int k, int *neighborArray,
                               unsigned int *dim) {
  int neighborIndex = 0;
  if (i <= 0) {
    neighborArray[neighborIndex++] = xyz2vtx(i + 1, j, k, dim);
  } else if (i >= (int)(dim[0]) - 1) {
    neighborArray[neighborIndex++] = xyz2vtx(i - 1, j, k, dim);
  } else {
    neighborArray[neighborIndex++] = xyz2vtx(i - 1, j, k, dim);
    neighborArray[neighborIndex++] = xyz2vtx(i + 1, j, k, dim);
  }
  if (j <= 0) {
    neighborArray[neighborIndex++] = xyz2vtx(i, j + 1, k, dim);
  } else if (j >= (int)(dim[1]) - 1) {
    neighborArray[neighborIndex++] = xyz2vtx(i, j - 1, k, dim);
  } else {
    neighborArray[neighborIndex++] = xyz2vtx(i, j - 1, k, dim);
    neighborArray[neighborIndex++] = xyz2vtx(i, j + 1, k, dim);
  }
  if (k <= 0) {
    neighborArray[neighborIndex++] = xyz2vtx(i, j, k + 1, dim);
  } else if (k >= (int)(dim[2]) - 1) {
    neighborArray[neighborIndex++] = xyz2vtx(i, j, k - 1, dim);
  } else {
    neighborArray[neighborIndex++] = xyz2vtx(i, j, k - 1, dim);
    neighborArray[neighborIndex++] = xyz2vtx(i, j, k + 1, dim);
  }
  return neighborIndex;
}
