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
#ifndef __Geoframe_H__
#define __Geoframe_H__

#include <assert.h>

#if !defined(__APPLE__)
#include <malloc.h>
#endif

#include <math.h>
#include <stdio.h>
#include <string>
#include <Utility/utility.h>

enum MESHTYPE { TRI, TRI2, TET, TET2, QUAD, HEX, NURBS, INVALIDMESH };
enum FILETYPE { RAW, RAWC, RAWN, RAWNC, INVALIDFILE };

class Geoframe {
public:
  Geoframe(const char *rawiv_fname);
  Geoframe();
  ~Geoframe();
  void LoadFile(const char *raw_fname);

  double testTetrahedron1(float v0[3], float v1[3], float v2[3], float v3[3]);
  float get_aspect_ratio(uint v0, uint v1, uint v2);
  float getRadius(float a[3][3], float b[3], float v0[3]);
  int AddQuad(uint *v, int num);
  int AddTri(uint v1, uint v2, uint v3);
  int AddVert(float v_pos[3], float norm[3]);
  int center_vtx(int v1, int v2, int v3);
  int CheckBound_edge(uint index_1, uint index_2);
  int getNumHexas(void);
  int getNumQuads(void);
  int getNumTris(void);
  int getNumVerts(void);
  void setNumVerts(int n);
  void setNumTris(int n);
  int TestNum(uint *v);
  int testRHS(float v0[3], float v1[3], float v2[3], float v3[3]);
  int testTetrahedron(float v0[3], float v1[3], float v2[3], float v3[3]);
  uint AddRefine_edgevtx_hex(uint index_1, uint index_2);
  uint AddRefine_edgevtx(uint index_1, uint index_2);
  uint AddRefine_facevtx(uint index_1, uint index_2, uint index_3,
                         uint index_4);
  void Add_2_Tetra(uint *v, uint my_vertex);
  void Add_2_Tri(uint *v);
  void AddBound_edge(uint index_1, uint index_2);
  void AddBound(int index, int sign);
  void AddPyramid(uint *v, uint my_vertex, int num);
  void AddQuad_adaptive_2_1(uint *v, uint *v_new, int num);
  void AddQuad_adaptive_2_3(uint *v, uint *v_new, int num);
  void AddQuad_adaptive_3_1(uint *v, uint *v_new, int num);
  void AddQuad_adaptive_3_2a(uint *v, uint *v_new, int num);
  void AddQuad_adaptive_3_2b(uint *v, uint *v_new, int num);
  void AddQuad_adaptive_3_3(uint *v, uint *v_new, int num);
  void AddQuad_adaptive_4_2b(uint *v, uint *v_new, int num);
  void AddQuad_adaptive_4(uint *v, uint *v_new, int num);
  void AddQuad_adaptive(uint *v, uint *v_new, int num);
  void AddQuad_hexa(uint *v, int num);
  void AddQuad_indirect(uint *v);
  void AddTetra(uint v1, uint v2, uint v3, uint v4);
  void AddVert_adaptive_2_1(uint *v, uint *v_new);
  void AddVert_adaptive_2_3(uint *v, uint *v_new);
  void AddVert_adaptive_3_1(uint *v, uint *v_new);
  void AddVert_adaptive_3_2a(uint *v, uint *v_new);
  void AddVert_adaptive_3_2b(uint *v, uint *v_new);
  void AddVert_adaptive_3_3(uint *v, uint *v_new);
  void AddVert_adaptive_4_2b(uint *v, uint *v_new);
  void AddVert_adaptive_4(uint *v, uint *v_new);
  void AddVert_adaptive(uint *v, uint *v_new);
  void AddVert_hexa_adaptive_1_center(uint *v, uint *v_new);
  void AddVert_hexa_adaptive_1_top(uint *v, uint *v_new);
  void AddVert_hexa_adaptive_2_1(uint *v, int *edge_id, uint *v_new);
  void AddVert_hexa_adaptive_2_2(uint *v, int *edge_id, uint *v_new);
  void AddVert_hexa_adaptive_2_4(uint *v, int *edge_id, uint *v_new);
  void AddVert_hexa_adaptive_2(uint *v, int *edge_id, uint *v_new);
  void AddVtxNew(int index, int sign);
  void calculateAspectRatio();
  void calculateExtents();
  void calculatenormals();
  void calculateTriangleNormal(float *norm, uint c);
  void Clear();
  void cross(float *dest, const float *v1, const float *v2);
  void edge_contraction_tetra(uint *v, int num);
  void edge_contraction_tri(uint *v);
  void edge_contraction(uint *v, int num);
  void Extend_Tetra(uint *v);
  void Extend_Tri(uint v0, uint v1, uint v2);
  void saveHexa(const char *filename);
  void saveQuad(const char *filename);
  void saveTetra(const char *filename);
  void saveTriangle(const char *filename);
  void setMin(float min0, float min1, float min2);
  void setSpan(float span0, float span1, float span2);
  void updateBySpan();
  void write_raw(const char *raw_fname, MESHTYPE meshtype);
  MESHTYPE getMeshType() { return meshType; }
  FILETYPE getFileType() { return fileType; }
  float avg_aspect, max_aspect, min_aspect;
  float (*curvatures)[2];
  float max_x, min_x;
  float max_y, min_y;
  float max_z, min_z;
  float span[3];
  int numQuads;
  int numTris;
  int numVerts;
  int (*refine_edge)[18];
  int tsize, vsize, qsize;
  int *vtx_idx_arr_extend;
  uint (*bound_edge)[18];
  uint *bound_tri;
  uint (*refine_edgevtx)[18];
  uint *vtxnew_sign;
  MESHTYPE meshType;
  FILETYPE fileType;
  void readVertex(int i, bool bound, FILE *fp);

public: // these are already ingrained in the code
  double biggestDim;
  double centerx, centery, centerz;
  float (*funcs)[1];
  float (*normals)[3];
  float (*verts)[3];
  float (*colors)[3];
  int numHexas;
  uint *bound_sign;
  uint (*quads)[4];
  uint (*triangles)[3];

private:
  void initializeGeoframe();
};

#endif
