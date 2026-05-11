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
#ifndef CONTOUR_H
#define CONTOUR_H

#include <Utility/utility.h>
#include <AdaptiveIsocontouring/comp.h>
#include <AdaptiveIsocontouring/Geoframe.h>

namespace AdaptiveIsocontouring {

class Contour {
public:
  Contour();
  ~Contour();
  void initialize(int fine_level, int coarse_level, float *vol, float *min,
                  float *max, LittleGeoframe &geo);
  void initialize(int fine_level, int coarse_level, float *vol, float *min,
                  float *max, LittleGeoframe *geo);
  void initialize(int fine_level, int coarse_level, double *vol,
                  LittleGeoframe &geo);
  void initialize_redbox(int resol, float *b_min, float *b_max, float *sum_in);
  void set_isovalue_outer(float iso);
  void set_isovalue_inner(float iso);
  void update_box(float *new_box, float *old_box);
  void update_func_val(float *func_val,
                       int grid_size); // gives the grid size in one dimension.
  void polygonize();
  void polygonize_redbox();
  void put_bishoulder();
  void put_bishoulder_redbox();

private:
  void update_quad_render();
  int fine_level;
  int coarse_level;
  int num_outer_quads;
  LittleGeoframe *geo;
  float iso_val_outer;
  float iso_val_inner;
  float *vd_coarse;
  float *vd_fine;
  float *bi_vert_in;
  float *bi_vert_out;
  float dim[3];
  float in_min[3];
  float in_max[3];
  // variables for the red_box
  int box_resol;
  float box_min[3];
  float box_max[3];
  float box_dim[3];
  int box_rel_min[3]; // gives the relative index in terms of box to coarse
  int box_rel_max[3];
};

} // namespace AdaptiveIsocontouring

#endif
