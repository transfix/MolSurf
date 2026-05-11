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
#if !defined(                                                                  \
    AFX_TRILINEARGRID_H__E1E51C41_6207_4C9F_AC66_FA5385A0E475__INCLUDED_)
#define AFX_TRILINEARGRID_H__E1E51C41_6207_4C9F_AC66_FA5385A0E475__INCLUDED_

#include <Utility/utility.h>

namespace CCVOpenGLMath {
// Used to go from vertex to grid indices and obtain neighbors.
// All functions in this class are static members
class TrilinearGrid {
public:
  // There is no need to create these objects as all member functions are static
  // in nature.
  TrilinearGrid();
  virtual ~TrilinearGrid();
  /* 	Given a cell index, we return the indices of the vertex at the first
     (lowest) corner. The last parameter contains the dimensions. \param
     cell_index is the index of a cell in a uniform 3d grid \param x, y, z are
     passed by reference and are returned as the indices of the lowest corner
          \param gdim is of length 3, containing the dimensions of the grid.
  */
  static void cell2xyz(int cell_index, int &x, int &y, int &z,
                       unsigned int *gdim);
  static int xyz2cell(int x, int y, int z, unsigned int *gdim);
  /*	A 1D vertex index is converted to a 3D set of indices
          \param vtx_idx is the index in 1D of a vertex
          \param x, y, z are passed by reference and are returned as the 3D
     indices \param gdim is of length 3, containing the dimensions of the grid.
  */
  static void vtx2xyz(int vtx_idx, int &x, int &y, int &z, unsigned int *gdim);
  static void getCellVertices(int cell_index, int *cellVertexArray,
                              unsigned int *gdim);
  static int xyz2vtx(int x, int y, int z, unsigned int *dim);
  static int getNeighbor(int i, int j, int k, int *neighborArray,
                         unsigned int *dim);
};
}; // namespace CCVOpenGLMath
#endif
