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
#ifndef DECODER_H
#define DECODER_H

#include <Utility/utility.h>
#include <C2C/arithlib/bitbuffer.h>
#include <C2C/arithlib/encode.h>
#include <C2C/C2CBuffer.h>
#include <C2C/ContourGeom.h>
#include <C2C/cubes.h>
#include <C2C/vtkMarchingCubesCases.h>

#ifdef WIN32
#include <C2C/arithlib/unitypes.h>
#else
#include <unistd.h>
#endif

#define VERBOSE

template <class T> class Decoder {
public:
  Decoder(C2CBuffer<T> *buf, bool color = false);
  ~Decoder();

  // Construct contour from compressed data
  ContourGeom *constructCon();

  void trackContour(ContourGeom *con, Cell *, int);

protected:
  int InterpEdge(ContourGeom *con, T *val, u_char *color, T grad[3][8], int i,
                 int j, int edge);
  void interpRect3Dpts_x(int i1, int j1, int k1, T *data, u_char *color,
                         T grad[3][8], int d1, int d2, float pt[3],
                         float norm[3], float clr[3]);
  void interpRect3Dpts_y(int i1, int j1, int k1, T *data, u_char *color,
                         T grad[3][8], int d1, int d2, float pt[3],
                         float norm[3], float clr[3]);
  void interpRect3Dpts_z(int i1, int j1, int k1, T *data, u_char *color,
                         T grad[3][8], int d1, int d2, float pt[3],
                         float norm[3], float clr[3]);

  void getCellValues(int i, int j, T *data, u_char *color);

  void getConfig(int nc, Cell *cells);

  T getValue(int i, int j, int v); // get the vth value of cell (i, j)

  void markNeigbor(int, int, int, int *, u_char);

  u_char getCode(int ix, int iy);

  bool isDegenerate(ContourGeom *con, int v1, int v2, int v3);

  // swap the first and the next slice
  void swapOrder() {
    T *td = first;
    first = next;
    next = td;
    u_char *tb = bf;
    bf = bn;
    bn = tb;
    u_char *tc = colf;
    colf = coln;
    coln = tc;
  }

private:
  int count; // current slice
  float value;
  int dim[3];
  float orig[3], span[3];
  T *first, *next; // current two slices
  u_char *colf, *coln;
  u_char *bf, *bn;
  C2CBuffer<T> *c2c_in;
};

template <class T>
Decoder<T>::Decoder(C2CBuffer<T> *buf, bool color) : c2c_in(buf) {
  value = c2c_in->getIsovalue();
  c2c_in->getDimension(dim, orig, span);
#ifdef _DEBUG
  cout << "dim: " << dim[0] << " " << dim[1] << " " << dim[2] << endl;
  cout << "orig: " << orig[0] << " " << orig[1] << " " << orig[2] << endl;
  cout << "span: " << span[0] << " " << span[1] << " " << span[2] << endl;
  cout << "isovalue = " << value << endl;
#endif
  count = 0;
  first = new T[dim[0] * dim[1]];
  next = new T[dim[0] * dim[1]];
  bf = new u_char[dim[0] * dim[1]];
  bn = new u_char[dim[0] * dim[1]];
  if (color) {
    colf = new u_char[3 * dim[0] * dim[1]];
    coln = new u_char[3 * dim[0] * dim[1]];
  } else {
    colf = 0;
    coln = 0;
  }
}

template <class T> Decoder<T>::~Decoder() {
  delete[] first;
  delete[] next;
  delete[] bf;
  delete[] bn;
  if (colf && coln) {
    delete[] colf;
    delete[] coln;
  }
}

template <class T> void Decoder<T>::getConfig(int nc, Cell *cells) {
  int i, j, k;
  for (k = 0; k < nc; k++) {
    i = cells[k].ix;
    j = cells[k].iy;
    // printf("cell (%d  %d) config= %d\n", cells[k].ix, cells[k].iy,
    // cells[k].code);
    cells[k].code = getCode(i, j);
    // printf("cell (%d  %d) config= %d\n", cells[k].ix, cells[k].iy,
    // cells[k].code);
  }
}

template <class T>
void Decoder<T>::markNeigbor(int ix, int iy, int vo, int *sign, u_char used) {
  u_char mark;
  for (int i = 0; i < 3; i++) {
    if (sign[vertneigbors[vo][i]] == 0) {
      mark = 1 & (used >> vertneigbors[vo][i]);
      if (mark) // the vertex is marked
      {
        sign[vertneigbors[vo][i]] =
            (getValue(ix, iy, vertneigbors[vo][i]) < value) ? -1 : 1;
      } else {
        assert(getValue(ix, iy, vertneigbors[vo][i]) == 0.0);
        sign[vertneigbors[vo][i]] = sign[vo];
      }
      markNeigbor(ix, iy, vertneigbors[vo][i], sign, used);
    }
  }
  return;
}

template <class T> T Decoder<T>::getValue(int ix, int iy, int v) {
  int offset = ix + iy * dim[0];
  switch (v) {
  case 0:
    return first[offset];
  case 1:
    return first[offset + 1];
  case 2:
    return next[offset + 1];
  case 3:
    return next[offset];
  case 4:
    return first[offset + dim[0]];
  case 5:
    return first[offset + dim[0] + 1];
  case 6:
    return next[offset + dim[0] + 1];
  case 7:
    return next[offset + dim[0]];
  }
  // should not come here
  assert(-1);
  return (0);
}

template <class T> u_char Decoder<T>::getCode(int ix, int iy) {
  u_char used = 0;
  int sign[8], vo;
  int offset = ix + iy * dim[0];
  if (bf[offset] == 1) {
    used |= 0x01;
    vo = 0;
  }
  if (bf[offset + 1] == 1) {
    used |= 0x02;
    vo = 1;
  }
  if (bn[offset + 1] == 1) {
    used |= 0x04;
    vo = 2;
  }
  if (bn[offset] == 1) {
    used |= 0x08;
    vo = 3;
  }
  if (bf[offset + dim[0]] == 1) {
    used |= 0x10;
    vo = 4;
  }
  if (bf[offset + dim[0] + 1] == 1) {
    used |= 0x20;
    vo = 5;
  }
  if (bn[offset + dim[0] + 1] == 1) {
    used |= 0x40;
    vo = 6;
  }
  if (bn[offset + dim[0]] == 1) {
    used |= 0x80;
    vo = 7;
  }
  // next to get the sign of very vertex
  memset(sign, 0, sizeof(int) * 8);
  assert(used != 0);
  if (getValue(ix, iy, vo) < value) {
    sign[vo] = -1;
  } else {
    sign[vo] = 1;
  }
  markNeigbor(ix, iy, vo, sign, used);
  u_char code = 0;
  for (int i = 0; i < 8; i++) {
    assert(sign[i] != 0);
    if (sign[i] == -1) {
      code |= 1 << i;
    }
  }
  return code;
}

template <class T> ContourGeom *Decoder<T>::constructCon() {
  int nc, nv, relev = 0;
  ContourGeom *con = new ContourGeom();
  Cell *cells;
  nv = c2c_in->getSlice(first, colf, bf);
  relev += nv;
  count++;
  nc = c2c_in->getLayer(cells);
  nv = c2c_in->getSlice(next, coln, bn);
  relev += nv;
  count++;
  trackContour(con, cells, nc);
  for (int i = 1; i < dim[2] - 1; i++) {
    swapOrder();
    delete[] cells;
    nc = c2c_in->getLayer(cells);
    nv = c2c_in->getSlice(next, coln, bn);
    count++;
    relev += nv;
    trackContour(con, cells, nc);
  }
  delete[] cells;
#ifdef _DEBUG
  printf("total number of relevant vertice: %d\n", relev);
#endif
  return con;
}

template <class T>
void Decoder<T>::trackContour(ContourGeom *con, Cell *cells, int nc) {
  int code, e, t, i;
  int edge;
  int edge_v[12];
  T val[8];
  u_char color[24];
  T grad[3][8];
  u_int v1, v2, v3;
  getConfig(nc, cells);
  for (i = 0; i < nc; i++) {
    code = cells[i].code;
    getCellValues(cells[i].ix, cells[i].iy, val, color);
    for (e = 0; e < cubeedges[code][0]; e++) {
      edge = cubeedges[code][1 + e];
      edge_v[edge] =
          InterpEdge(con, val, color, grad, cells[i].ix, cells[i].iy, edge);
    }
    /*for(t=0; triCases[code].edges[t] != -1; ) {
      v1 = edge_v[triCases[code].edges[t++]];
      v2 = edge_v[triCases[code].edges[t++]];
      v3 = edge_v[triCases[code].edges[t++]];
      if(v1 != v2 && v1 != v3 && v2 != v3) {
    con->AddTri(v1, v2, v3);
      }
      }*/
    for (t = 0; t < cubes[code][0]; t++) {
      v1 = edge_v[cubes[code][t * 3 + 1]];
      v2 = edge_v[cubes[code][t * 3 + 2]];
      v3 = edge_v[cubes[code][t * 3 + 3]];
#ifdef UNIQUE
      if (v1 != v2 && v1 != v3 && v2 != v3) {
        con->AddTri(v1, v2, v3);
      }
#else
      if (!isDegenerate(con, v1, v2, v3)) {
        con->AddTri(v1, v2, v3);
      }
#endif
    }
  }
}

template <class T>
bool Decoder<T>::isDegenerate(ContourGeom *con, int v1, int v2, int v3) {
  /*
  if ((con->vert[v1][0] == con->vert[v2][0]) &&
          (con->vert[v1][1] == con->vert[v2][1]) &&
          (con->vert[v1][2] == con->vert[v2][2]))	 return true;
  else if ((con->vert[v1][0] == con->vert[v3][0]) &&
                   (con->vert[v1][1] == con->vert[v3][1]) &&
                   (con->vert[v1][2] == con->vert[v3][2])) return true;
  else if ((con->vert[v2][0] == con->vert[v3][0]) &&
                   (con->vert[v2][1] == con->vert[v3][1]) &&
                   (con->vert[v2][2] == con->vert[v3][2])) return true;
  */
  if (v1 == v2 || v1 == v3 || v2 == v3) {
    return true;
  }
  return false;
}

template <class T>
int Decoder<T>::InterpEdge(ContourGeom *con, T *val, u_char *color,
                           T grad[3][8], int i, int j, int edge) {
  float pt[3];
  float norm[3];
  float clr[3];
  int k = count - 2;
  EdgeInfo *ei = &edgeinfo[edge];
  int v;
  switch (ei->dir) {
  case 0:
    interpRect3Dpts_x(i + ei->di, j + ei->dj, k + ei->dk, val, color, grad,
                      ei->d1, ei->d2, pt, norm, clr);
    break;
  case 1:
    interpRect3Dpts_y(i + ei->di, j + ei->dj, k + ei->dk, val, color, grad,
                      ei->d1, ei->d2, pt, norm, clr);
    break;
  case 2:
    // assert(ei->dk == 0);
    interpRect3Dpts_z(i + ei->di, j + ei->dj, k + ei->dk, val, color, grad,
                      ei->d1, ei->d2, pt, norm, clr);
    break;
  }
#ifdef UNIQUE
  EdgeIndex eindex(i + ei->di, j + ei->dj, k + ei->dk, ei->dir);
  v = con->AddVertUnique(pt, norm, clr, eindex);
  // v = con->AddVertUnique(pt, norm, 0);
#else
  v = con->AddVert(pt, norm, clr);
#endif
  return v;
}

template <class T>
void Decoder<T>::interpRect3Dpts_x(int i1, int j1, int k1, T *data,
                                   u_char *color, T grad[3][8], int d1, int d2,
                                   float pt[3], float norm[3], float clr[3]) {
  double ival;
  if (data[d2] == data[d1]) {
    printf("x1 = %f, x2 = %f, isoval = %f\n", (float)data[d1], (float)data[d2],
           value);
  }
#ifdef RECONSTRUCT
  if (data[d2] == data[d1]) {
    ival = 0.5;
  } else {
    ival = (value - (float)data[d1]) / ((float)data[d2] - (float)data[d1]);
  }
#else
  ival = (value - (float)data[d1]) / ((float)data[d2] - (float)data[d1]);
#endif
  // ival = 0.5;
  pt[0] = orig[0] + span[0] * (i1 + ival);
  pt[1] = orig[1] + span[1] * j1;
  pt[2] = orig[2] + span[2] * k1;
  clr[0] =
      (color[d1 * 3 + 0] * (1.0 - ival) + color[d2 * 3 + 0] * ival) / 255.0f;
  clr[1] =
      (color[d1 * 3 + 1] * (1.0 - ival) + color[d2 * 3 + 1] * ival) / 255.0f;
  clr[2] =
      (color[d1 * 3 + 2] * (1.0 - ival) + color[d2 * 3 + 2] * ival) / 255.0f;
#ifdef _DEBUG
  // printf("interpRect3Dpts_x(%d,%d,%d)\n", i1, j1, k1);
#endif
}

template <class T>
void Decoder<T>::interpRect3Dpts_y(int i1, int j1, int k1, T *data,
                                   u_char *color, T grad[3][8], int d1, int d2,
                                   float pt[3], float norm[3], float clr[3]) {
  double ival;
#ifdef RECONSTRUCT
  if (data[d2] == data[d1]) {
    ival = 0.5;
  } else {
    ival = (value - (float)data[d1]) / ((float)data[d2] - (float)data[d1]);
  }
#else
  ival = (value - (float)data[d1]) / ((float)data[d2] - (float)data[d1]);
#endif
  // printf("value = %f: %f %f\n", value, data[d1], data[d2]);
  // ival = 0.5;
  pt[0] = orig[0] + span[0] * i1;
  pt[1] = orig[1] + span[1] * (j1 + ival);
  pt[2] = orig[2] + span[2] * k1;
  clr[0] =
      (color[d1 * 3 + 0] * (1.0 - ival) + color[d2 * 3 + 0] * ival) / 255.0f;
  clr[1] =
      (color[d1 * 3 + 1] * (1.0 - ival) + color[d2 * 3 + 1] * ival) / 255.0f;
  clr[2] =
      (color[d1 * 3 + 2] * (1.0 - ival) + color[d2 * 3 + 2] * ival) / 255.0f;
#ifdef _DEBUG
  // printf("interpRect3Dpts_y(%d,%d,%d)\n", i1, j1, k1);
#endif
}

template <class T>
void Decoder<T>::interpRect3Dpts_z(int i1, int j1, int k1, T *data,
                                   u_char *color, T grad[3][8], int d1, int d2,
                                   float pt[3], float norm[3], float clr[3]) {
  double ival;
  if (data[d2] == data[d1]) {
    printf("x1 = %f, x2 = %f, isoval = %f\n", (float)data[d1], (float)data[d2],
           value);
  }
#ifdef RECONSTRUCT
  if (data[d2] == data[d1]) {
    ival = 0.5;
  } else {
    ival = (value - (float)data[d1]) / ((float)data[d2] - (float)data[d1]);
  }
#else
  ival = (value - (float)data[d1]) / ((float)data[d2] - (float)data[d1]);
#endif
  pt[0] = orig[0] + span[0] * i1;
  pt[1] = orig[1] + span[1] * j1;
  pt[2] = orig[2] + span[2] * (k1 + ival);
  clr[0] =
      (color[d1 * 3 + 0] * (1.0 - ival) + color[d2 * 3 + 0] * ival) / 255.0f;
  clr[1] =
      (color[d1 * 3 + 1] * (1.0 - ival) + color[d2 * 3 + 1] * ival) / 255.0f;
  clr[2] =
      (color[d1 * 3 + 2] * (1.0 - ival) + color[d2 * 3 + 2] * ival) / 255.0f;
#ifdef _DEBUG
  // printf("interpRect3Dpts_z(%d,%d,%d)\n", i1, j1, k1);
#endif
}

template <class T>
void Decoder<T>::getCellValues(int i, int j, T data[8], u_char color[24]) {
  int offset = i + j * dim[0];
  data[0] = first[offset];
  data[1] = first[offset + 1];
  data[2] = next[offset + 1];
  data[3] = next[offset];
  data[4] = first[offset + dim[0]];
  data[5] = first[offset + dim[0] + 1];
  data[6] = next[offset + dim[0] + 1];
  data[7] = next[offset + dim[0]];
  if (colf && coln) {
    color[0] = colf[(offset) * 3 + 0];
    color[1] = colf[(offset) * 3 + 1];
    color[2] = colf[(offset) * 3 + 2];
    color[3] = colf[(offset + 1) * 3 + 0];
    color[4] = colf[(offset + 1) * 3 + 1];
    color[5] = colf[(offset + 1) * 3 + 2];
    color[6] = coln[(offset + 1) * 3 + 0];
    color[7] = coln[(offset + 1) * 3 + 1];
    color[8] = coln[(offset + 1) * 3 + 2];
    color[9] = coln[(offset) * 3 + 0];
    color[10] = coln[(offset) * 3 + 1];
    color[11] = coln[(offset) * 3 + 2];
    color[12] = colf[(offset + dim[0]) * 3 + 0];
    color[13] = colf[(offset + dim[0]) * 3 + 1];
    color[14] = colf[(offset + dim[0]) * 3 + 2];
    color[15] = colf[(offset + dim[0] + 1) * 3 + 0];
    color[16] = colf[(offset + dim[0] + 1) * 3 + 1];
    color[17] = colf[(offset + dim[0] + 1) * 3 + 2];
    color[18] = coln[(offset + dim[0] + 1) * 3 + 0];
    color[19] = coln[(offset + dim[0] + 1) * 3 + 1];
    color[20] = coln[(offset + dim[0] + 1) * 3 + 2];
    color[21] = coln[(offset + dim[0]) * 3 + 0];
    color[22] = coln[(offset + dim[0]) * 3 + 1];
    color[23] = coln[(offset + dim[0]) * 3 + 2];
  } else {
    memset(color, 0, 24);
  }
}

#endif
