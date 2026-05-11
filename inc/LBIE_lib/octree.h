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
#ifndef __OCTREE_H__
#define __OCTREE_H__

#include <LBIE_lib/contour3d.h>
#include <LBIE_lib/e_face.h>
#include <LBIE_lib/Geoframe.h>
#include <SimpleVolumeData/SimpleVolumeData.h>
#include <stdio.h>

typedef struct _octcell {
  char refine_flag;
} octcell;

typedef struct _MinMax {
  float min;
  float max;
} MinMax;

#define NUM_PARENT_CELL 6
#define FLOAT_MINIMUM -10000000
#define FLOAT_MAXIMUM 10000000

typedef struct {
  int dir;
  int di, dj, dk;
  int d1, d2;
} EdgeInfo;

static int po[8][6][3] = {
    {{0, -1, -1}, {-1, 0, -1}, {0, 0, -1}, {-1, -1, 0}, {0, -1, 0}, {-1, 0, 0}},
    {{0, -1, -1}, {0, 0, -1}, {1, 0, -1}, {0, -1, 0}, {1, -1, 0}, {1, 0, 0}},
    {{0, 1, -1}, {-1, 0, -1}, {0, 0, -1}, {-1, 1, 0}, {0, 1, 0}, {-1, 0, 0}},
    {{0, 1, -1}, {0, 0, -1}, {1, 0, -1}, {0, 1, 0}, {1, 1, 0}, {1, 0, 0}},
    {{0, -1, 1}, {-1, 0, 1}, {0, 0, 1}, {-1, -1, 0}, {0, -1, 0}, {-1, 0, 0}},
    {{0, -1, 1}, {0, 0, 1}, {1, 0, 1}, {0, -1, 0}, {1, -1, 0}, {1, 0, 0}},
    {{0, 1, 1}, {-1, 0, 1}, {0, 0, 1}, {-1, 1, 0}, {0, 1, 0}, {-1, 0, 0}},
    {{0, 1, 1}, {0, 0, 1}, {1, 0, 1}, {0, 1, 0}, {1, 1, 0}, {1, 0, 0}}};

static EdgeInfo edgeinfo[12] = {
    {0, 0, 0, 0, 0, 1}, {2, 1, 0, 0, 1, 2}, {0, 0, 0, 1, 3, 2},
    {2, 0, 0, 0, 0, 3}, {0, 0, 1, 0, 4, 5}, {2, 1, 1, 0, 5, 6},
    {0, 0, 1, 1, 7, 6}, {2, 0, 1, 0, 4, 7}, {1, 0, 0, 0, 0, 4},
    {1, 1, 0, 0, 1, 5}, {1, 0, 0, 1, 3, 7}, {1, 1, 0, 1, 2, 6}};

static EdgeInfo edge_dir[6][4] = {{{1, 1, 0, 0, 0, 4},
                                   {0, 1, 1, 0, 4, 1},
                                   {1, 1, 1, 0, 4, 2},
                                   {0, 0, 1, 0, 3, 4}},
                                  {{2, 1, 0, 0, 0, 4},
                                   {0, 1, 0, 1, 4, 1},
                                   {2, 1, 0, 1, 4, 2},
                                   {0, 0, 0, 1, 3, 4}},
                                  {{1, 2, 0, 1, 0, 4},
                                   {2, 2, 1, 0, 1, 4},
                                   {1, 2, 1, 1, 4, 2},
                                   {2, 2, 1, 1, 4, 3}},
                                  {{2, 1, 2, 0, 0, 4},
                                   {0, 1, 2, 1, 4, 1},
                                   {2, 1, 2, 1, 4, 2},
                                   {0, 0, 2, 1, 3, 4}},
                                  {{1, 0, 0, 1, 0, 4},
                                   {2, 0, 1, 0, 1, 4},
                                   {1, 0, 1, 1, 4, 2},
                                   {2, 0, 1, 1, 4, 3}},
                                  {{1, 1, 0, 2, 0, 4},
                                   {0, 1, 1, 2, 4, 1},
                                   {1, 1, 1, 2, 4, 2},
                                   {0, 0, 1, 2, 3, 4}}};

static int level_id[] = {0,
                         1,
                         1 + 8,
                         1 + 8 + 64,
                         1 + 8 + 64 + 512,
                         1 + 8 + 64 + 512 + 4096,
                         1 + 8 + 64 + 512 + 4096 + 32768,
                         1 + 8 + 64 + 512 + 4096 + 32768 + 262144,
                         1 + 8 + 64 + 512 + 4096 + 32768 + 262144 + 2097152,
                         1 + 8 + 64 + 512 + 4096 + 32768 + 262144 + 2097152 +
                             16777216};

static int interp_2vtx[18][5] = {
    {2, 0, 1},       {2, 0, 2},       {4, 0, 1, 2, 3}, {2, 1, 3},
    {2, 2, 3},

    {2, 0, 4},       {4, 0, 1, 4, 5}, {2, 1, 5},       {4, 0, 2, 4, 6},
    {4, 1, 3, 5, 7}, {2, 2, 6},       {4, 2, 3, 6, 7}, {2, 3, 7},

    {2, 4, 5},       {2, 4, 6},       {4, 4, 5, 6, 7}, {2, 5, 7},
    {2, 6, 7}};

static int cube_eid[12][2] = {{0, 1}, {1, 2}, {2, 3}, {0, 3}, {4, 5}, {5, 6},
                              {6, 7}, {4, 7}, {0, 4}, {1, 5}, {3, 7}, {2, 6}};

class Octree {

private:
  //    Contour3d curcon;

  FILE *vol_fp;
  float iso_val, iso_val_in;
  int leaf_num;
  octcell *oct_array;
  int octcell_num;
  int cell_num;
  int oct_depth;
  int level_res[10];
  int *cut_array;
  MESHTYPE meshType;
  int flag_normal;

  int in_out;
  int flag_extend;

  edge_face e_face[12][12];

  double **qef_array;
  double **qef_array_in;
  int *vtx_idx_arr;
  int *vtx_idx_arr_in;
  int *grid_idx_arr;
  int *vtx_idx_arr_refine;

  float *orig_vol;
  char *ebit;
  char *vbit;
  float *BSplineCoeff;

  MinMax *minmax;

  SimpleVolumeData *volumeData;

  int prop_flag;
  int interior_flag;
  char prop_fname[1000];

  int nverts, ncells;
  int dim[3];
  int olddim[3];
  float orig[3];
  float span[3];

  int get_neighbor_bit(int id, int level);
  int is_refined(int x, int y, int z, int level);
  int is_refined2(int x, int y, int z, int level);
  int get_depth(int res);
  int get_octcell_num(int depth);
  int get_level(int oc_id);
  void idx2vtx(int oc_id, int level, int *vtx, int *interpol_vtx);
  void idx2vtx(int oc_id, int level, int *vtx);
  void interpRect3Dpts_x(int i1, int j1, int k1, float f1, float f2, float val,
                         float *pt, float norm[3], int level);
  void interpRect3Dpts_y(int i1, int j1, int k1, float f1, float d2, float val,
                         float *pt, float norm[3], int level);
  void interpRect3Dpts_z(int i1, int j1, int k1, float f1, float d2, float val,
                         float *pt, float norm[3], int level);
  void e_face_initialization();

  float compute_error(int oc_idx, int level, float &min, float &max);
  void construct_octree(char *);

  int child(int oc_id, int level, int i);

  void octcell2xyz(int oc_id, int &x, int &y, int &z, int level);
  int xyz2octcell(int x, int y, int z, int level);
  int xyz2vtx(int x, int y, int z);
  void getCellValues(int oc_id, int level, float *val);
  int is_intersect(int in_idx, float isovalue, float *in_value, int &iv_idx,
                   int i, int j, int k, int level, int faceidx,
                   Geoframe &geofrm);

  int is_intersect(float *val, int e_id);
  int is_intersect_interval(float *val, int e_id);

  void getVertGrad(int i, int j, int k, float g[3]);
  float getValue(int i, int j, int k);

  void add_tetra_face(int oc_id, int level, Geoframe &geofrm);
  void march_one_face(int face_id, int oc_id, int level, unsigned int *my_vtx,
                      unsigned int center, Geoframe &geofrm);
  void march_each_face(int oc_id, int level, int face_id, unsigned int center,
                       Geoframe &geofrm);
  void march_each_edge(int oc_id, int level, int edge_id, int *my_edge,
                       Geoframe &geofrm);
  void add_tetra_face_interval(int oc_id, int level, Geoframe &geofrm);

  void add_tetra_cube(int oc_id, int level, Geoframe &geofrm);
  void add_tetra_cube_adaptive(int oc_id, int level, Geoframe &geofrm);
  void add_tetra_hexa(int x, int y, int z, int *vtx, Geoframe &geofrm);

  void add_hexa(Geoframe &geofrm, unsigned int *vtx);

  void add_cube(int oc_id, int level, Geoframe &geofrm);
  int march_hexa_face(int oc_id, int level, int face_id, int *vtx,
                      Geoframe &geofrm);
  void add_one_vertex(int x, int y, int z, int cell_size,
                      unsigned int &my_vertex, Geoframe &geofrm);
  void add_middle_vertex(int x, int y, int z, float dx, float dy, float dz,
                         int cell_size, unsigned int &my_vertex,
                         Geoframe &geofrm);

  void march_edge(int x, int y, int z, int cell_size, int e_id, int num,
                  int *temp, int *index, int &m_id, unsigned int *my_array,
                  Geoframe &geofrm);
  void get_index_array(int level, int &num, int *index);

  void face_0(int x, int y, int z, int cell_size, int face_id, int v0, int v1,
              int v2, int v3, int v4, Geoframe &geofrm);
  void face_1(int v0, int v1, int v2, int v3, int v4, unsigned int *array,
              int leng, Geoframe &geofrm);
  void face_2_0(int x, int y, int z, int face_id, int v0, int v1, int v2,
                int v3, int v4, unsigned int *array0, unsigned int *array1,
                int m_id0, int m_id1, Geoframe &geofrm);
  void face_2_1(int v0, int v1, int v2, int v3, int v4, unsigned int *array0,
                unsigned int *array1, int m_id0, int m_id1, Geoframe &geofrm);
  void face_3(int x, int y, int z, int face_id, int cell_size, int v0, int v1,
              int v2, int v3, int v4, unsigned int *array0,
              unsigned int *array1, unsigned int *array2, int m_id0, int m_id1,
              int m_id2, Geoframe &geofrm);
  void face_4(int x, int y, int z, int face_id, int v0, int v1, int v2, int v3,
              int v4, unsigned int face_center, unsigned int *array0,
              unsigned int *array1, unsigned int *array2, unsigned int *array3,
              int m_id0, int m_id1, int m_id2, int m_id3, Geoframe &geofrm);

  void permute_1(int &i0, int &i1, int &i2, int &i3, int s0, int s1, int s2,
                 int s3);
  void permute_2(int &i0, int &i1, int &i2, int &i3, int &s0, int &s1, int &s2,
                 int &s3);
  void permute_3(int &i0, int &i1, int &i2, int &i3, int s0, int s1, int s2,
                 int s3);

  void get_middle_array_1(int face_id, int *my_edge0, int *my_edge1,
                          int *my_edge2, int *my_edge3, unsigned int *my_array,
                          int &m_id, int x, int y, int z, int level,
                          Geoframe &geofrm);
  void get_middle_array_2(int face_id, int *my_edge0, int *my_edge1,
                          int *my_edge2, int *my_edge3, unsigned int *my_array0,
                          unsigned int *my_array1, int &m_id0, int &m_id1,
                          int x, int y, int z, int level, Geoframe &geofrm);
  void get_middle_array_3(int face_id, int *my_edge0, int *my_edge1,
                          int *my_edge2, int *my_edge3, unsigned int *my_array0,
                          unsigned int *my_array1, unsigned int *my_array2,
                          int &m_id0, int &m_id1, int &m_id2, int x, int y,
                          int z, int level, Geoframe &geofrm);
  void get_middle_array_4(int face_id, int *my_edge0, int *my_edge1,
                          int *my_edge2, int *my_edge3, unsigned int *my_array0,
                          unsigned int *my_array1, unsigned int *my_array2,
                          unsigned int *my_array3, int &m_id0, int &m_id1,
                          int &m_id2, int &m_id3, unsigned int &face_center,
                          int x, int y, int z, int level, Geoframe &geofrm);

  int is_skipcell(int oc_id);
  int is_skipcell_in(int oc_id);
  int is_skipcell_interval(int oc_id);
  void get_solution(int oc_id, float *pos);
  void get_vtx(int x, int y, int z, int level, float *pos);
  void get_VtxNorm(float *vtx, float *norm);
  void get_min_vertex(int e_id, int intersect_id, int x, int y, int z, int &tx,
                      int &ty, int &tz);
  void get_vertex_index(int v_id, int x, int y, int z, int &tx, int &ty,
                        int &tz);

  int cell_comp(int oc_id, int level, float pt[12][3], float norm[12][3]);
  int cell_comp_in(int oc_id, int level, float pt[12][3], float norm[12][3]);
  void put_qef(int oc_id, double *sigma_ni_2, double *sigma_ni_2_pi,
               double *sigma_ni_2_pi_2, double *solution, double qef);
  void put_qef_in(int oc_id, double *sigma_ni_2, double *sigma_ni_2_pi,
                  double *sigma_ni_2_pi_2, double *solution, double qef);
  void get_qef(int child_idx, double *temp_sigma_ni_2,
               double *temp_sigma_ni_2_pi, double *temp_sigma_ni_2_pi_2);
  void get_qef_in(int child_idx, double *temp_sigma_ni_2,
                  double *temp_sigma_ni_2_pi, double *temp_sigma_ni_2_pi_2);
  float get_err(int oc_id);
  float get_err_grad(int oc_id);
  float get_err_grad_test(int oc_id);

  void eflag_clear();
  int is_eflag_on(int x, int y, int z, int level, int e);
  void eflag_on(int x, int y, int z, int level, int e);

  void vflag_clear();
  int is_vflag_on(int x, int y, int z, int level, int v);
  void vflag_on(int x, int y, int z, int level, int v);
  int is_mf_vflag_on(int x, int y, int z, int level, int f);
  void mf_vflag_on(int x, int y, int z, int level, int f);

  int is_min_edge(int oc_id, int e_id, unsigned int *vtx, int &vtx_num,
                  int intersect_id, Geoframe &geofrm);
  int is_min_edge_2(int oc_id, int e_id, int *vtx, int &vtx_num, int *con_id,
                    int intersect_id, Geoframe &geofrm);
  int is_min_vertex(int oc_id, int v_id, unsigned int *vtx, Geoframe &geofrm);

  int min_vtx(int x, int y, int z, int level, Geoframe &geofrm);
  int min_vtx_tetra(int x, int y, int z, int e_id_0, int e_id, int level,
                    int &con_ind, int &con_id_1, Geoframe &geofrm);
  int min_vtx_hexa(int x, int y, int z, int level, Geoframe &geofrm);

  void clear(double *a, double *b, double *c);
  void clear(double *a);

public:
  Octree();
  // construction of octree from the given volume

  ~Octree();

  float vol_min, vol_max;

  float spans[3];
  float minext[3], maxext[3];

  // isosurface simplification operations
  void set_isovalue(float val) { iso_val = val; }
  void set_isovalue_in(float val) { iso_val_in = val; }
  void setMeshType(MESHTYPE meshType_) { meshType = meshType_; }
  void setNormalType(int n) { flag_normal = n; }

  // dual contouring operations
  void collapse();
  void collapse_interval();
  void compute_qef();
  void compute_qef_interval();
  void traverse_qef(float err_tol);
  void traverse_qef_interval(float err_tol, float err_tol_in);
  void mesh_extract(Geoframe &geofrm, float err_tol);
  void polygonize(Geoframe &geofrm);
  void polygonize_interval(Geoframe &geofrm);

  void polygonize_quad(Geoframe &geofrm, float err_tol);
  void tetrahedralize(Geoframe &geofrm);
  void hexahedralize(Geoframe &geofrm, float err_tol);
  void hexa_adaptive_1(Geoframe &geofrm, int *oc_id, int *edge_id,
                       float err_tol, unsigned int *vtx);
  void hexa_adaptive_1_top(Geoframe &geofrm, unsigned int *vtx,
                           unsigned int *vtx_new);
  void hexa_adaptive_2(Geoframe &geofrm, int *oc_id, int *edge_id,
                       float err_tol, unsigned int *vtx);
  void add_hexa_adaptive_2(Geoframe &geofrm, unsigned int *vtx_new);
  void add_hexa_adaptive_2_1(Geoframe &geofrm, unsigned int *vtx,
                             unsigned int *vtx_temp);
  void add_hexa_adaptive_2_2(Geoframe &geofrm, unsigned int *vtx_new);
  void add_hexa_adaptive_2_4(Geoframe &geofrm, unsigned int *vtx_new);
  void tetra_to_4_hexa(Geoframe &geofrm);
  void tetrahedralize_interval(Geoframe &geofrm);
  void quality_improve(Geoframe &geofrm, int improvemethod);
  void quality_improve(Geoframe &geofrm, int improvemethod, float delta_t);
  void func_val(Geoframe &geofrm);
  void find_oc_id(int x, int y, int z, int level, int j, int intersect_id,
                  int *oc_id);
  void find_oc_id_hexa(int x, int y, int z, int level, int j, int *oc_id);
  void find_edge_id_hexa(int x, int y, int z, int cell_size, int j,
                         int *edge_id);
  void find_vtx_new(Geoframe &geofrm, int x, int y, int z, int level, int j,
                    int intersect_id, unsigned int *vtx_new);
  void quad_adaptive(Geoframe &geofrm, int *oc_id, float err_tol,
                     unsigned int *vtx, int flag_method);
  void quad_adaptive_method1(Geoframe &geofrm, int *oc_id, float err_tol,
                             unsigned int *vtx);
  void quad_adaptive_method2(Geoframe &geofrm, int *oc_id, float err_tol,
                             unsigned int *vtx);
  void quad_adaptive_method3(Geoframe &geofrm, int *oc_id, float err_tol,
                             unsigned int *vtx, int flag_method);
  void quad_adaptive_method4(Geoframe &geofrm, int *oc_id, float err_tol,
                             unsigned int *vtx);
  void quad_adaptive_method5(Geoframe &geofrm, int *oc_id, float err_tol,
                             unsigned int *vtx);
  void get_vtx_new(Geoframe &geofrm, int oc_id, unsigned int vtx);
  void check_topo();

  void setInEx() { interior_flag = !interior_flag; }
  void prop_init(const char *rawiv_fname);
  void Octree_loadVolume(const char *dx_fname);
  void initMem();
  void pdbVolume(const char *pdb_fname);

  void setVolume(SimpleVolumeData *volData);
  void loadVolume();
  void loadData(SimpleVolumeData *volumeData, float *vol);
  bool isInterior() { return interior_flag; }
  void flipInterior() { interior_flag = !interior_flag; }

  // quality improve
  void smoothing_joeliu_volume(Geoframe &geofrm, int sign);
  void geometric_flow(Geoframe &geofrm);
  void geometric_flow(Geoframe &geofrm, float delta_t);
  void geometric_flow_tri(Geoframe &geofrm);
  void geometric_flow_tri(Geoframe &geofrm, float delta_t);
  void geometric_flow_tet(Geoframe &geofrm);
  void geometric_flow_quad(Geoframe &geofrm);
  void geometric_flow_hex(Geoframe &geofrm);
  void crossproduct(float v0[3], float v1[3], float v2[3], float *normal);
  float area_tri(float v0[3], float v1[3], float v2[3]);
  float area_quad(float v0[3], float v1[3], float v2[3], float v3[3]);
  float volume_tet(float v0[3], float v1[3], float v2[3], float v3[3]);
  float volume_hex(float v0[3], float v1[3], float v2[3], float v3[3],
                   float v4[3], float v5[3], float v6[3], float v7[3]);
  void edge_contraction(Geoframe &geofrm);
  void optimization(Geoframe &geofrm);
  void testHexa(float v0[3], float v1[3], float v2[3], float v3[3],
                float &jacob, float &condition_no, float &oddy);
  void getGrad(float v0[3], float v1[3], float v2[3], float v3[3],
               float &condition_no, float *grad);

  void assign_refine_sign_quad(Geoframe &geofrm, float err_tol);
  void assign_refine_sign_hexa(Geoframe &geofrm, float err_tol);
};

#endif //__OCTREE_H__
