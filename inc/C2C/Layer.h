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
// Layer: class of one layer of cells

#ifndef C2C_LAYER_H
#define C2C_LAYER_H

#include <Utility/utility.h>
#include <C2C/arithlib/bitbuffer.h>
#include <C2C/arithlib/encode.h>
#include <C2C/cubes.h>

#ifdef WIN32
#include <C2C/arithlib/unitypes.h>
#else
#include <unistd.h>
#endif

class Layer {
private:
  // dimension of cell layer , use number of vertices
  // in each dimension
  int dim[2];
  int nc;        // # of cells in the layer
  u_char *codes; // code of cells

public:
  Layer(int d1, int d2, u_char *cells = NULL);
  Layer(const Layer &lay);
  ~Layer();

  u_char *getCells(void) const {
    u_char *c_codes =
        (u_char *)malloc(sizeof(u_char) * (dim[0] - 1) * (dim[1] - 1));
    memcpy(c_codes, codes, (dim[0] - 1) * (dim[1] - 1));
    return c_codes;
  }
  int *getDimen(void) const { return (int *)dim; }
  int getNC(void) const { return nc; }
  // diffBits: return difference of two cell layer bitmap
  BIT *diffBits(Layer &lay);
  void writeOut(FILE *fp);
};

#endif
