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
#ifndef C2C_SINGLECONSTEP_H
#define C2C_SINGLECONSTEP_H

#include <Utility/utility.h>
#include <C2C/Layer.h>
#include <C2C/Slice.h>

template <class T> class SingleConStep {
private:
  int dim[3];        // dimension of the mesh
  int c_slice;       // current number of slice
  int c_layer;       // current number of layer
  Slice<T> **slices; // slices of the mesh
  Layer **layers;    // cell layers of the mesh

public:
  SingleConStep(int dx, int dy, int dz);
  ~SingleConStep();

  void addSlice(Slice<T> &sl);
  void addLayer(Layer &lay);
  void info(void);

  // write out a singleconstep without differing from another one
  void writeSingleCS(float val, FILE *ofp);
  // write the difference between two singleconsteps
  void writeDiffCS(const SingleConStep &prev, float val, FILE *ofp);
};

template <class T> SingleConStep<T>::SingleConStep(int dx, int dy, int dz) {
  int i;
  dim[0] = dx;
  dim[1] = dy;
  dim[2] = dz;
  c_slice = c_layer = 0;
  slices = (Slice<T> **)malloc(sizeof(Slice<T> *) * dim[2]);
  for (i = 0; i < dim[2]; i++) {
    slices[i] = NULL;
  }
  layers = (Layer **)malloc(sizeof(Layer *) * (dim[2] - 1));
  for (i = 0; i < dim[2] - 1; i++) {
    layers[i] = NULL;
  }
}

template <class T> SingleConStep<T>::~SingleConStep() {
  if (slices != NULL) {
    for (int i = 0; i < dim[2]; i++) {
      delete slices[i];
    }
    free(slices);
  }
  if (layers != NULL) {
    for (int j = 0; j < dim[2] - 1; j++) {
      delete layers[j];
    }
    free(layers);
  }
}

template <class T> void SingleConStep<T>::addSlice(Slice<T> &sl) {
  slices[c_slice++] = new Slice<T>(sl);
  // printf("# of vertices = %d, %d\n", sl.nv, slices[c_slice-1]->nv);
}

template <class T> void SingleConStep<T>::addLayer(Layer &lay) {
  layers[c_layer++] = new Layer(lay);
}

template <class T> void SingleConStep<T>::info() {
  fprintf(stderr, "%d slices & %d layers\n", c_slice, c_layer);
}

template <class T> void SingleConStep<T>::writeSingleCS(float val, FILE *ofp) {
  int i;
  for (i = 0; i < dim[2] - 1; i++) {
    slices[i]->writeOut(val); // slice inherits out fp from contour object
    layers[i]->writeOut(ofp);
  }
  slices[dim[2] - 1]->writeOut(val); // the last slice
}

// another way of doing temporal differential coding
template <class T>
void SingleConStep<T>::writeDiffCS(const SingleConStep &prev, float val,
                                   FILE *ofp) {
  BIT *sl_diff, *lay_diff;
  BitBuffer *vbuf = new BitBuffer();
  BitBuffer *cbuf = new BitBuffer();
  int count_c = 0, count_v = 0, count_vnew = 0;
  T *fvals = new T[dim[0] * dim[1] * dim[2]]; // function value array
  T *fvals_new = new T[dim[0] * dim[1] * dim[2]];
  for (int i = 0; i < dim[2]; i++) {
    sl_diff = slices[i]->diffBits(*(prev.slices[i]));
    vbuf->put_bits(dim[0] * dim[1], sl_diff);
    for (int j = 0; j < dim[0] * dim[1]; j++) {
      if (slices[i]->verts[j].isUsed()) {
        if (prev.slices[i]->verts[j].isUsed())
          fvals[count_v++] = (slices[i]->verts[j].getValue() >=
                              prev.slices[i]->verts[j].getValue())
                                 ? slices[i]->verts[j].getValue() -
                                       prev.slices[i]->verts[j].getValue()
                                 : prev.slices[i]->verts[j].getValue() -
                                       slices[i]->verts[j].getValue();
        else {
          fvals_new[count_vnew++] = slices[i]->verts[j].getValue();
        }
      }
    }
    delete[] sl_diff;
    if (i != dim[2] - 1) // not the last slice
    {
      lay_diff = layers[i]->diffBits(*(prev.layers[i]));
      cbuf->put_bits((dim[0] - 1) * (dim[1] - 1), lay_diff);
      delete[] lay_diff;
    }
  }
#ifdef ZP_CODEC
  vbuf->zp_encode();
#else
  vbuf->arith_encode();
#endif
  vbuf->writeFile(ofp);
  delete vbuf;
  BitBuffer *valbuf = encode_vals(fvals, count_v, val);
  valbuf->writeFile(ofp);
  delete valbuf;
  valbuf = encode_vals(fvals_new, count_vnew, val);
  valbuf->writeFile(ofp);
  delete valbuf;
  delete[] fvals;
#ifdef ZP_CODEC
  cbuf->zp_encode();
#else
  cbuf->arith_encode();
#endif
  cbuf->writeFile(ofp);
  delete cbuf;
}
#endif
