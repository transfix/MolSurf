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
#include <Blurmaps/AreaVolume.h>

AreaVolume::AreaVolume() {}

AreaVolume::~AreaVolume() {}

double AreaVolume::sqr(double a) { return a * a; }

void AreaVolume::cross_product(vtx *triangle, double *normal) {
  vtx v1, v2;
  double abs_val;
  v1.x = triangle[1].x - triangle[0].x;
  v1.y = triangle[1].y - triangle[0].y;
  v1.z = triangle[1].z - triangle[0].z;
  v2.x = triangle[2].x - triangle[0].x;
  v2.y = triangle[2].y - triangle[0].y;
  v2.z = triangle[2].z - triangle[0].z;
  normal[0] = v1.y * v2.z - v1.z * v2.y;
  normal[1] = v1.z * v2.x - v1.x * v2.z;
  normal[2] = v1.x * v2.y - v1.y * v2.x;
  abs_val = sqrt(normal[0] * normal[0] + normal[1] * normal[1] +
                 normal[2] * normal[2]);
  if (abs_val < 0.0000001) {
    normal[0] = 0;
    normal[1] = 0;
    normal[2] = 0;
    normal[3] = 0;
  } else {
    normal[0] /= abs_val;
    normal[1] /= abs_val;
    normal[2] /= abs_val;
    normal[3] = -(normal[0] * triangle[0].x + normal[1] * triangle[0].y +
                  normal[2] * triangle[0].z);
  }
}

double AreaVolume::tri_area(vtx x0, vtx x1, vtx x2) {
  double area;
  vtx vec1, vec2, cp;
  vec1.x = x1.x - x0.x;
  vec1.y = x1.y - x0.y;
  vec1.z = x1.z - x0.z;
  vec2.x = x2.x - x0.x;
  vec2.y = x2.y - x0.y;
  vec2.z = x2.z - x0.z;
  cp.x = vec1.y * vec2.z - vec1.z * vec2.y;
  cp.y = vec1.z * vec2.x - vec1.x * vec2.z;
  cp.z = vec1.x * vec2.y - vec1.y * vec2.x;
  area = 0.5 * fabs(sqrt(sqr(cp.x) + sqr(cp.y) + sqr(cp.z)));
  return area;
}

double AreaVolume::tri_area_p(vtx x0, vtx x1, vtx x2) {
  double area;
  vtx vec1, vec2, cp;
  vec1.x = x1.x - x0.x;
  vec1.y = x1.y - x0.y;
  vec1.z = x1.z - x0.z;
  vec2.x = x2.x - x0.x;
  vec2.y = x2.y - x0.y;
  vec2.z = x2.z - x0.z;
  cp.x = vec1.y * vec2.z - vec1.z * vec2.y;
  cp.y = vec1.z * vec2.x - vec1.x * vec2.z;
  cp.z = vec1.x * vec2.y - vec1.y * vec2.x;
  area = 0.5 * fabs(sqrt(sqr(cp.x) + sqr(cp.y) + sqr(cp.z)));
  return area;
}

double AreaVolume::tetra_vol(vtx x0, vtx x1, vtx x2, vtx x3) {
  double eq[4], vol;
  double triangle_area;
  double d;
  vtx tri[3];
  tri[0] = x0;
  tri[1] = x1;
  tri[2] = x2;
  triangle_area = tri_area(x0, x1, x2);
  cross_product(tri, eq);
  d = fabs((double)(eq[0] * x3.x + eq[1] * x3.y + eq[2] * x3.z + eq[3]));
  vol = d * triangle_area / 3.;
  return vol;
}

double AreaVolume::getRatio(double a, double b, float isoval) {
  if (b == a) {
    return 0.5;
  }
  return (isoval - a) / (b - a);
}

void AreaVolume::getPts(vtx &inter_v, vtx x1, vtx x2, float ival) {
  inter_v.x = x1.x + (x2.x - x1.x) * ival;
  inter_v.y = x1.y + (x2.y - x1.y) * ival;
  inter_v.z = x1.z + (x2.z - x1.z) * ival;
}

void AreaVolume::interval_vol(vtx x[4], double v[4], float isovalue,
                              double *Area, double *Volume) {
  int code;
  int plus_no = 0;
  int in_idx = 0;
  int out_idx = 0;
  int in_vidx[4];
  int out_vidx[4];
  vtx inter_vtx[4];
  int i;
  double ival;
  code = 0;
  vtx _t;
  vtx x11, x12, x13, x21, x22, x23;
  double t;
  double V1, V2, V3, V, A;
  if (v[3] < v[2]) {
    _t = x[3];
    x[3] = x[2];
    x[2] = _t;
    t = v[3];
    v[3] = v[2];
    v[2] = t;
  }
  if (v[2] < v[1]) {
    _t = x[2];
    x[2] = x[1];
    x[1] = _t;
    t = v[2];
    v[2] = v[1];
    v[1] = t;
  }
  if (v[1] < v[0]) {
    _t = x[1];
    x[1] = x[0];
    x[0] = _t;
    t = v[1];
    v[1] = v[0];
    v[0] = t;
  }
  if (v[3] < v[2]) {
    _t = x[3];
    x[3] = x[2];
    x[2] = _t;
    t = v[3];
    v[3] = v[2];
    v[2] = t;
  }
  if (v[2] < v[1]) {
    _t = x[2];
    x[2] = x[1];
    x[1] = _t;
    t = v[2];
    v[2] = v[1];
    v[1] = t;
  }
  if (v[3] < v[2]) {
    _t = x[3];
    x[3] = x[2];
    x[2] = _t;
    t = v[3];
    v[3] = v[2];
    v[2] = t;
  }
  if (v[0] >= isovalue) {
    code += 0x01;
    in_vidx[in_idx++] = 0;
  } else {
    out_vidx[out_idx++] = 0;
  }
  if (v[1] >= isovalue) {
    code += 0x02;
    in_vidx[in_idx++] = 1;
  } else {
    out_vidx[out_idx++] = 1;
  }
  if (v[2] >= isovalue) {
    code += 0x04;
    in_vidx[in_idx++] = 2;
  } else {
    out_vidx[out_idx++] = 2;
  }
  if (v[3] >= isovalue) {
    code += 0x08;
    in_vidx[in_idx++] = 3;
  } else {
    out_vidx[out_idx++] = 3;
  }
  switch (in_idx) {
  case 0:
    V = 0;
    A = 0;
    break;
  case 1:
    for (i = 0; i < 3; i++) {
      ival = getRatio(v[in_vidx[0]], v[out_vidx[i]], isovalue);
      getPts(inter_vtx[i], x[in_vidx[0]], x[out_vidx[i]], (float)ival);
    }
    V = tetra_vol(inter_vtx[0], inter_vtx[1], inter_vtx[2], x[in_vidx[0]]);
    // V=tetra_vol(x[0],x[1],x[2],x[3]);
    // V=1/6.
    A = tri_area(inter_vtx[0], inter_vtx[1], inter_vtx[2]);
    break;
  case 2:
    for (i = 0; i < 2; i++) {
      ival = getRatio(v[in_vidx[0]], v[out_vidx[i]], isovalue);
      getPts(inter_vtx[i], x[in_vidx[0]], x[out_vidx[i]], (float)ival);
    }
    for (i = 0; i < 2; i++) {
      ival = getRatio(v[in_vidx[1]], v[out_vidx[i]], isovalue);
      getPts(inter_vtx[2 + i], x[in_vidx[1]], x[out_vidx[i]], (float)ival);
    }
    x11 = x[in_vidx[0]];
    x12 = inter_vtx[0];
    x13 = inter_vtx[1];
    x21 = x[in_vidx[1]];
    x22 = inter_vtx[2];
    x23 = inter_vtx[3];
    V1 = tetra_vol(x11, x12, x13, x21);
    V2 = tetra_vol(x21, x12, x13, x22);
    V3 = tetra_vol(x21, x13, x22, x23);
    V = V1 + V2 + V3;
    // V=tetra_vol(x[0],x[1],x[2],x[3]);
    // V=0;
    A = tri_area(x12, x22, x13) + tri_area(x23, x22, x13);
    break;
  case 3:
    for (i = 0; i < 3; i++) {
      ival = getRatio(v[in_vidx[i]], v[out_vidx[0]], isovalue);
      // getPts(inter_vtx[i],x[out_vidx[0]],x[in_vidx[i]],ival);
      getPts(inter_vtx[i], x[in_vidx[i]], x[out_vidx[0]], (float)ival);
    }
    V1 = tetra_vol(inter_vtx[0], inter_vtx[1], inter_vtx[2], x[out_vidx[0]]);
    V2 = tetra_vol(x[0], x[1], x[2], x[3]);
    V = V2 - V1;
    // V=tetra_vol(x[0],x[1],x[2],x[3]);
    // V=0;
    A = tri_area(inter_vtx[0], inter_vtx[1], inter_vtx[2]);
    break;
  case 4:
    V = tetra_vol(x[0], x[1], x[2], x[3]);
    A = 0;
    break;
  }
  *Volume = V;
  *Area = A;
}

void AreaVolume::cell2index(int c, int &x, int &y, int &z) {
  x = c % (dim[0] - 1);
  y = (c / (dim[0] - 1)) % (dim[1] - 1);
  z = c / ((dim[0] - 1) * (dim[1] - 1));
}

int AreaVolume::idx2vtx(int i, int j, int k) {
  if (((i < 0) || (i >= dim[0])) || ((j < 0) || (j >= dim[1])) ||
      ((k < 0) || (k >= dim[2]))) {
    return -1;
  }
  return i + (j * dim[0]) + (k * dim[0] * dim[1]);
}

void AreaVolume::getCellVerts(int c, int i, int v[4], double x0[3],
                              double x1[3], double x2[3], double x3[3]) {
  int x, y, z;
  cell2index(c, x, y, z);
  switch (i) {
  case 0:
    v[0] = idx2vtx(x, y, z + 1);
    v[1] = idx2vtx(x + 1, y, z + 1);
    v[2] = idx2vtx(x + 1, y + 1, z + 1);
    v[3] = idx2vtx(x + 1, y + 1, z);
    x0[0] = x;
    x0[1] = y;
    x0[2] = z + 1;
    x1[0] = x + 1;
    x1[1] = y;
    x1[2] = z + 1;
    x2[0] = x + 1;
    x2[1] = y + 1;
    x2[2] = z + 1;
    x3[0] = x + 1;
    x3[1] = y + 1;
    x3[2] = z;
    break;
  case 1:
    v[0] = idx2vtx(x, y, z + 1);
    v[1] = idx2vtx(x + 1, y, z);
    v[2] = idx2vtx(x + 1, y, z + 1);
    v[3] = idx2vtx(x + 1, y + 1, z);
    x0[0] = x;
    x0[1] = y;
    x0[2] = z + 1;
    x1[0] = x + 1;
    x1[1] = y;
    x1[2] = z;
    x2[0] = x + 1;
    x2[1] = y;
    x2[2] = z + 1;
    x3[0] = x + 1;
    x3[1] = y + 1;
    x3[2] = z;
    break;
  case 2:
    v[0] = idx2vtx(x, y, z + 1);
    v[2] = idx2vtx(x, y + 1, z + 1);
    v[1] = idx2vtx(x + 1, y + 1, z + 1);
    v[3] = idx2vtx(x + 1, y + 1, z);
    x0[0] = x;
    x0[1] = y;
    x0[2] = z + 1;
    x1[0] = x;
    x1[1] = y + 1;
    x1[2] = z + 1;
    x2[0] = x + 1;
    x2[1] = y + 1;
    x2[2] = z + 1;
    x3[0] = x + 1;
    x3[1] = y + 1;
    x3[2] = z;
    break;
  case 3:
    v[0] = idx2vtx(x, y, z + 1);
    v[1] = idx2vtx(x, y + 1, z + 1);
    v[2] = idx2vtx(x, y + 1, z);
    v[3] = idx2vtx(x + 1, y + 1, z);
    x0[0] = x;
    x0[1] = y;
    x0[2] = z + 1;
    x1[0] = x;
    x1[1] = y + 1;
    x1[2] = z + 1;
    x2[0] = x;
    x2[1] = y + 1;
    x2[2] = z;
    x3[0] = x + 1;
    x3[1] = y + 1;
    x3[2] = z;
    break;
  case 4:
    v[0] = idx2vtx(x, y, z + 1);
    v[1] = idx2vtx(x, y, z);
    v[2] = idx2vtx(x + 1, y, z);
    v[3] = idx2vtx(x + 1, y + 1, z);
    x0[0] = x;
    x0[1] = y;
    x0[2] = z + 1;
    x1[0] = x;
    x1[1] = y;
    x1[2] = z;
    x2[0] = x + 1;
    x2[1] = y;
    x2[2] = z;
    x3[0] = x + 1;
    x3[1] = y + 1;
    x3[2] = z;
    break;
  case 5:
    v[0] = idx2vtx(x, y, z + 1);
    v[2] = idx2vtx(x, y, z);
    v[1] = idx2vtx(x, y + 1, z);
    v[3] = idx2vtx(x + 1, y + 1, z);
    x0[0] = x;
    x0[1] = y;
    x0[2] = z + 1;
    x1[0] = x;
    x1[1] = y;
    x1[2] = z;
    x2[0] = x;
    x2[1] = y + 1;
    x2[2] = z;
    x3[0] = x + 1;
    x3[1] = y + 1;
    x3[2] = z;
    break;
  }
}

void AreaVolume::Transform(double x[3]) {
  int i;
  for (i = 0; i < 3; i++) {
    x[i] = (double)(orig[i] + x[i] * span[i]);
  }
}

int check_total_inside;

bool AreaVolume::getTetCellAreaAndVolume(float *vol, double isovalue, int c,
                                         double *sum_A, double *sum_V,
                                         double *app_vol_lower,
                                         double *app_vol_upper,
                                         double *app_boundary) {
  int i;
  for (i = 0; i < 6; i++) {
    int v[4];
    vtx x[4];
    double val[4];
    double x0[3], x1[3], x2[3], x3[3];
    getCellVerts(c, i, v, x0, x1, x2, x3);
    Transform(x0);
    Transform(x1);
    Transform(x2);
    Transform(x3);
    x[0].x = x0[0];
    x[0].y = x0[1];
    x[0].z = x0[2];
    x[1].x = x1[0];
    x[1].y = x1[1];
    x[1].z = x1[2];
    x[2].x = x2[0];
    x[2].y = x2[1];
    x[2].z = x2[2];
    x[3].x = x3[0];
    x[3].y = x3[1];
    x[3].z = x3[2];
    val[0] = vol[v[0]];
    val[1] = vol[v[1]];
    val[2] = vol[v[2]];
    val[3] = vol[v[3]];
    double V, A;
    interval_vol(x, val, (float)isovalue, &A, &V);
    // assert(V>=0 && V <=1/6);
    if (V < -0.0001 || V > (1. / 3.) * 1.7 * 1.7 * 1.7) {
      // fprintf(stderr,"abnormal volume[%d] : %lf\n",count_volume++,V);
      V = 0;
    }
    if (A < -0.0001 || A > 1.7 * 1.7) {
      // fprintf(stderr,"abnormal area[%d] : %lf\n",count_area++,A);
      A = 0;
    }
    // if (V>0.00001) printf("V : %lf\n",V);
    *sum_A += A;
    *sum_V += V;
    if (val[0] > isovalue && val[1] > isovalue && val[2] > isovalue &&
        val[3] > isovalue) {
      *app_vol_lower += 1 / 6.;
    } else if (val[0] > isovalue || val[1] > isovalue || val[2] > isovalue ||
               val[3] > isovalue) {
      *app_vol_upper += 1 / 6.;
    } else {
      *app_boundary += 1 / 6.;
    }
  }
  return true;
}

void AreaVolume::getVolume(float *vol, float isovalue, int *_dim, float *_orig,
                           float *_span, double *Area, double *Volume) {
  if (!vol || !_dim || !_orig || !_span || !Area || !Volume) {
    return;
  }
  dim[0] = _dim[0];
  dim[1] = _dim[1];
  dim[2] = _dim[2];
  orig[0] = _orig[0];
  orig[1] = _orig[1];
  orig[2] = _orig[2];
  span[0] = _span[0];
  span[1] = _span[1];
  span[2] = _span[2];
  nverts = dim[0] * dim[1] * dim[2];
  ncells = (dim[0] - 1) * (dim[1] - 1) * (dim[2] - 1);
  int c;
  int count_volume, count_area;
  double sum_V = 0;
  double sum_A = 0;
  double app_vol_lower = 0;
  double app_vol_upper = 0;
  double app_boundary = 0;
  count_volume = 0;
  count_area = 0;
  double singleCellVolume = span[0] * span[1] * span[2];
  double volBelow = 0;
  check_total_inside = 0;
  for (c = 0; c < ncells; c++) {
    getTetCellAreaAndVolume(vol, isovalue, c, &sum_A, &sum_V, &app_vol_lower,
                            &app_vol_upper, &app_boundary);
    if (c % 10000 == 0) {
      printf("%lf\n", sum_V);
    }
  }
  printf("total inside = %d\n", check_total_inside);
  *Volume = sum_V;
  *Area = volBelow;
}
