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
#if !defined(AFX_AREAVOLUME_H__8548386B_9854_46C5_9988_D549CC2B1646__INCLUDED_)
#define AFX_AREAVOLUME_H__8548386B_9854_46C5_9988_D549CC2B1646__INCLUDED_

#include <Utility/utility.h>
#include <UsefulMath/Matrix.h>
#include <UsefulMath/TrilinearGrid.h>

typedef struct _vtx {
  double x;
  double y;
  double z;
} vtx;

class AreaVolume {
public:
  AreaVolume();
  virtual ~AreaVolume();
  void getVolume(float *vol, float isovalue, int *_dim, float *_orig,
                 float *_span, double *Area, double *Volume);

protected:
  bool getTetCellAreaAndVolume(float *vol, double isovalue, int c,
                               double *sum_A, double *sum_V,
                               double *app_vol_lower, double *app_vol_upper,
                               double *app_boundary);
  bool getCellVolume(float *vol, double isovalue, int c, double *volAbove,
                     double *volBelow, double singleCellVolume);
  double sqr(double a);
  void cross_product(vtx *triangle, double *normal);
  double tri_area(vtx x0, vtx x1, vtx x2);
  double tri_area_p(vtx x0, vtx x1, vtx x2);
  double tetra_vol(vtx x0, vtx x1, vtx x2, vtx x3);
  double getRatio(double a, double b, float isoval);
  void getPts(vtx &inter_v, vtx x1, vtx x2, float ival);
  void interval_vol(vtx x[4], double v[4], float isovalue, double *Area,
                    double *Volume);
  void cell2index(int c, int &x, int &y, int &z);
  int idx2vtx(int i, int j, int k);
  void getCellVerts(int c, int i, int v[4], double x0[3], double x1[3],
                    double x2[3], double x3[3]);
  void Transform(double x[3]);
  int dim[3];
  float orig[3];
  float span[3];
  int nverts, ncells;
};

#endif
