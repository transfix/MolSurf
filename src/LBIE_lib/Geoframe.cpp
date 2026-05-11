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
#include <LBIE_lib/Geoframe.h>
#include <Utility/utility.h>

using namespace CVCUtility;

void Geoframe::setNumVerts(int n) { numVerts = n; }

void Geoframe::setNumTris(int n) { numTris = n; }

// Code common to all constructors
void Geoframe::initializeGeoframe() {
  min_x = min_y = min_z = 0;
  span[0] = span[1] = span[2] = 1;
  numVerts = 0;
  numTris = 0;
  numQuads = 0;
  numHexas = 0;
  vsize = 3000000;
  tsize = 3000000;
  qsize = 3000000;
  verts = (float(*)[3])mallocSafely(sizeof(float[3]) * vsize);
  funcs = (float(*)[1])mallocSafely(sizeof(float[1]) * vsize);
  normals = (float(*)[3])mallocSafely(sizeof(float[3]) * tsize);
  curvatures = (float(*)[2])mallocSafely(sizeof(float[2]) * tsize);
  triangles = (uint(*)[3])mallocSafely(sizeof(uint[3]) * tsize);
  quads = (uint(*)[4])mallocSafely(sizeof(uint[4]) * qsize);
  bound_sign = (uint(*))mallocSafely(sizeof(uint) * vsize);
  bound_tri = (uint(*))mallocSafely(sizeof(uint) * tsize);
  vtxnew_sign = (uint(*))mallocSafely(sizeof(uint) * vsize);
  bound_edge = (uint(*)[18])mallocSafely(sizeof(uint[18]) * vsize);
  refine_edge = (int(*)[18])mallocSafely(sizeof(int[18]) * vsize);
  refine_edgevtx = (uint(*)[18])mallocSafely(sizeof(uint[18]) * vsize);
  colors = (float(*)[3])mallocSafely(sizeof(float[3]) * vsize);
  vtx_idx_arr_extend = (int *)mallocSafely(sizeof(int) * vsize);
  for (int k = 0; k < 3000000; k++) {
    vtx_idx_arr_extend[k] = -1;
  }
}

Geoframe::Geoframe() { initializeGeoframe(); }

Geoframe::~Geoframe() {
  free(triangles);
  free(quads);
  free(verts);
  free(funcs);
  free(normals);
  free(curvatures);
  free(bound_sign);
  free(bound_tri);
  free(vtx_idx_arr_extend);
  free(vtxnew_sign);
  free(bound_edge);
  free(refine_edge);
  free(refine_edgevtx);
}

// Cross product
void Geoframe::cross(float *dest, const float *v1, const float *v2) {
  dest[0] = v1[1] * v2[2] - v1[2] * v2[1];
  dest[1] = v1[2] * v2[0] - v1[0] * v2[2];
  dest[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

void Geoframe::calculateAspectRatio() {
  avg_aspect = 0;
  max_aspect = 0;
  printf("calculate aspect ratio\n");
  int vert1, vert2, vert3;
  float tmp_aspect = 0;
  for (int c = 0; c < numTris; c++) {
    vert1 = triangles[c][0];
    vert2 = triangles[c][1];
    vert3 = triangles[c][2];
    tmp_aspect = get_aspect_ratio(vert1, vert2, vert3);
    if (tmp_aspect > max_aspect) {
      max_aspect = tmp_aspect;
    }
    avg_aspect += tmp_aspect;
  }
  avg_aspect /= (float)numTris;
  printf("max aspect ratio:%f avg aspect ratio:%f\n", max_aspect, avg_aspect);
}

void Geoframe::calculateTriangleNormal(float *norm, uint c) {
  float v1[3], v2[3];
  int vert = triangles[c][0];
  v1[0] = v2[0] = -verts[vert][0];
  v1[1] = v2[1] = -verts[vert][1];
  v1[2] = v2[2] = -verts[vert][2];
  vert = triangles[c][1];
  v1[0] += verts[vert][0];
  v1[1] += verts[vert][1];
  v1[2] += verts[vert][2];
  vert = triangles[c][2];
  v2[0] += verts[vert][0];
  v2[1] += verts[vert][1];
  v2[2] += verts[vert][2];
  cross(norm, v1, v2);
}

void Geoframe::calculatenormals() {
  float normal[3];
  // for each triangle
  for (int c = 0; c < numTris; c++) {
    calculateTriangleNormal(normal, c);
    normals[c][0] = -normal[0];
    normals[c][1] = -normal[1];
    normals[c][2] = -normal[2];
  }
  // normalize the vectors
  for (int vert = 0; vert < numTris; vert++) {
    float len = (float)sqrt(normals[vert][0] * normals[vert][0] +
                            normals[vert][1] * normals[vert][1] +
                            normals[vert][2] * normals[vert][2]);
    if (len == 0.0) {
      normals[vert][0] = 1.0;
      normals[vert][1] = 0.0;
      normals[vert][2] = 0.0;
      cout << "error" << endl;
      continue;
    }
    normals[vert][0] /= len;
    normals[vert][1] /= len;
    normals[vert][2] /= len;
  }
}

void Geoframe::calculateExtents() {
  max_x = 0, min_x = 0;
  max_y = 0, min_y = 0;
  max_z = 0, min_z = 0;
  float value;
  for (int c = 0; c < numVerts; c++) {
    if (c == 0) {
      max_x = min_x = verts[c][0];
      max_y = min_y = verts[c][1];
      max_z = min_z = verts[c][2];
    } else {
      value = verts[c][0];
      max_x = (value > max_x ? value : max_x);
      min_x = (value < min_x ? value : min_x);
      value = verts[c][1];
      max_y = (value > max_y ? value : max_y);
      min_y = (value < min_y ? value : min_y);
      value = verts[c][2];
      max_z = (value > max_z ? value : max_z);
      min_z = (value < min_z ? value : min_z);
    }
  }
  biggestDim = (max_y - min_y > max_x - min_x ? max_y - min_y : max_x - min_x);
  biggestDim = (max_z - min_z > biggestDim ? max_z - min_z : biggestDim);
  centerx = (max_x + min_x) / 2.0;
  centery = (max_y + min_y) / 2.0;
  centerz = (max_z + min_z) / 2.0;
}

// Helper function that reads a single vertex in from the file.
void Geoframe::readVertex(int i, bool bound, FILE *file) {
  int size;
  if (bound) {
    switch (fileType) {
    case RAWNC:
      size = fscanf(file, "%f %f %f %f %f %f %f %f %f %d\n", &verts[i][0],
                    &verts[i][1], &verts[i][2], &normals[i][0], &normals[i][1],
                    &normals[i][2], &colors[i][0], &colors[i][1], &colors[i][2],
                    &bound_sign[i]);
      if (size != 10)
        error("Wrong number of arguments read!");
      break;
    case RAWN:
      size = fscanf(file, "%f %f %f %f %f %f %d\n", &verts[i][0], &verts[i][1],
                    &verts[i][2], &normals[i][0], &normals[i][1],
                    &normals[i][2], &bound_sign[i]);
      colors[i][0] = 1.0;
      colors[i][1] = 0.7;
      colors[i][2] = 0.0;
      if (size != 7)
        error("Wrong number of arguments read!");
      break;
    case RAWC:
      size = fscanf(file, "%f %f %f %f %f %f %d\n",

                    &verts[i][0], &verts[i][1], &verts[i][2], &colors[i][0],
                    &colors[i][1], &colors[i][2], &bound_sign[i]);
      normals[i][0] = 0;
      normals[i][1] = 0;
      normals[i][2] = 0;
      if (size != 7)
        error("Wrong number of arguments read!");
      break;
    case RAW:
      size = fscanf(file, "%f %f %f %d\n", &verts[i][0], &verts[i][1],
                    &verts[i][2], &bound_sign[i]);
      normals[i][0] = 0;
      normals[i][1] = 0;
      normals[i][2] = 0;
      colors[i][0] = 1.0;
      colors[i][1] = 0.7;
      colors[i][2] = 0.0;
      if (size != 4)
        error("Wrong number of arguments read!");
      break;
    default:
      error("Did not recognize the file type. Exiting.");
    }
  } else {
    bound_sign[i] = 1;
    switch (fileType) {
    case RAWNC:
      size =
          fscanf(file, "%f %f %f %f %f %f %f %f %f\n", &verts[i][0],
                 &verts[i][1], &verts[i][2], &normals[i][0], &normals[i][1],
                 &normals[i][2], &colors[i][0], &colors[i][1], &colors[i][2]);
      if (size != 9)
        error("Wrong number of arguments read!");
      break;
    case RAWN:
      size =
          fscanf(file, "%f %f %f %f %f %f\n", &verts[i][0], &verts[i][1],
                 &verts[i][2], &normals[i][0], &normals[i][1], &normals[i][2]);
      colors[i][0] = 1.0;
      colors[i][1] = 0.7;
      colors[i][2] = 0.0;
      if (size != 6)
        error("Wrong number of arguments read!");
      break;
    case RAWC:
      size = fscanf(file, "%f %f %f %f %f %f\n", &verts[i][0], &verts[i][1],
                    &verts[i][2], &colors[i][0], &colors[i][1], &colors[i][2]);
      normals[i][0] = 0;
      normals[i][1] = 0;
      normals[i][2] = 0;
      if (size != 6)
        error("Wrong number of arguments read!");
      break;
    case RAW:
      size =
          fscanf(file, "%f %f %f\n", &verts[i][0], &verts[i][1], &verts[i][2]);
      normals[i][0] = 0;
      normals[i][1] = 0;
      normals[i][2] = 0;
      colors[i][0] = 1.0;
      colors[i][1] = 0.7;
      colors[i][2] = 0.0;
      if (size != 3)
        error("Wrong number of arguments read!");
      break;
    }
  }
}

Geoframe::Geoframe(const char *rawivFileName) {

  initializeGeoframe();

  LoadFile(rawivFileName);
}

void Geoframe::LoadFile(const char *rawivFileName) {
  // Open file
  FILE *file = fopen(rawivFileName, "r");
  if (file == NULL) {
    cout << "ERROR: Could not open file: " << rawivFileName << endl;
    exit(-1);
  }
  // Get mesh type
  meshType = INVALIDMESH;
  if (strstr(rawivFileName, "_tri"))
    meshType = TRI;
  if (strstr(rawivFileName, "_tet"))
    meshType = TET;
  if (strstr(rawivFileName, "_tri2"))
    meshType = TRI2;
  if (strstr(rawivFileName, "_tet2"))
    meshType = TET2;
  if (strstr(rawivFileName, "_quad"))
    meshType = QUAD;
  if (strstr(rawivFileName, "_hex"))
    meshType = HEX;
  if (strstr(rawivFileName, "_nurbs"))
    meshType = NURBS;
  if (meshType == INVALIDMESH) {
    cout << "WARNING: Filename does not have a tri/tet/quad/hex label." << endl;
    cout << "WARNING: assuming tri mesh type." << endl;
    meshType = TRI;
  }
  // Get file type
  fileType = INVALIDFILE;
  // arand: added else statements below... previously every file was being
  // treated as a RAW file
  if (strstr(rawivFileName, ".rawnc"))
    fileType = RAWNC;
  else if (strstr(rawivFileName, ".rawn"))
    fileType = RAWN;
  else if (strstr(rawivFileName, ".rawc"))
    fileType = RAWC;
  else if (strstr(rawivFileName, ".raw"))
    fileType = RAW;
  if (fileType == INVALIDFILE) {
    cout << "ERROR: Incorrect file extension: " << rawivFileName << endl;
    exit(-1);
  }

  float r = 0;
  float g = 0;
  float b = 0;
  int **neighbor = 0;
  // Parse depending on mesh type
  switch (meshType) {
  case TRI2: // Craig: I think geoframe doesn't care...
  case TRI:  // read triangular mesh
    if (fscanf(file, "%d %d\n", &numVerts, &numTris) != 2)
      error("Couldn't read number of vertices and triangles!");
    verts = (float(*)[3])reallocSafely(verts, sizeof(float[3]) * numVerts);
    normals = (float(*)[3])reallocSafely(
        normals, sizeof(float[3]) * minimum(numVerts, numTris));
    triangles = (uint(*)[3])reallocSafely(triangles, sizeof(uint[3]) * numTris);
    colors = (float(*)[3])reallocSafely(colors, sizeof(float[3]) * numVerts);
    // Craig: changed these from realloc to reallocSafely, added a color init
    bound_sign = (uint(*))reallocSafely(bound_sign, sizeof(uint) * vsize);
    bound_tri = (uint(*))reallocSafely(bound_tri, sizeof(uint) * numTris);
    for (int i = 0; i < numVerts; i++) {
      readVertex(i, false, file);
    }
    for (int i = 0; i < numTris; i++) {
      // bound_tri[i] = 1;
      if (fscanf(file, "%d %d %d\n", &triangles[i][0], &triangles[i][1],
                 &triangles[i][2]) != 3)
        error("Error reading triangle!");
    }
    break;
  case TET2: // Craig: I think geoframe doesn't care...
  case TET:  // read tetra mesh
    if (fscanf(file, "%d %d\n", &numVerts, &numTris) != 2)
      error("Couldn't read number of vertices and triangles!");
    numTris = numTris * 4;
    verts = (float(*)[3])reallocSafely(verts, sizeof(float[3]) * numVerts);
    colors = (float(*)[3])reallocSafely(colors, sizeof(float[3]) * numVerts);
    normals = (float(*)[3])reallocSafely(normals, sizeof(float[3]) * numVerts);
    triangles = (uint(*)[3])reallocSafely(triangles, sizeof(uint[3]) * numTris);
    bound_sign = (uint(*))reallocSafely(bound_sign, sizeof(uint) * vsize);
    bound_tri = (uint(*))reallocSafely(bound_tri, sizeof(uint) * numTris);
    for (int i = 0; i < numVerts; i++) {
      readVertex(i, true, file);

      verts[i][0] = verts[i][0] + 32;
      verts[i][1] = verts[i][1] + 32;
      verts[i][2] = verts[i][2] + 32;
    }
    neighbor = new int *[numVerts];
    for (int i = 0; i < numVerts; i++) {
      neighbor[i] = 0;
    }
    for (int i = 0; i < numVerts; i++) {
      neighbor[i] = new int[50];
    }
    for (int i = 0; i < numVerts; i++) {
      for (int j = 0; j < 50; j++) {
        neighbor[i][j] = -1;
      }
    }
    for (int i = 0; i < numTris / 4; i++) {
      int t0, t1, t2, t3;
      if (fscanf(file, "%d %d %d %d\n", &t0, &t1, &t2, &t3) != 4)
        error("Error reading quad!");
      triangles[4 * i + 0][0] = t0;
      triangles[4 * i + 0][1] = t2;
      triangles[4 * i + 0][2] = t1;
      triangles[4 * i + 1][0] = t1;
      triangles[4 * i + 1][1] = t2;
      triangles[4 * i + 1][2] = t3;
      triangles[4 * i + 2][0] = t0;
      triangles[4 * i + 2][1] = t3;
      triangles[4 * i + 2][2] = t2;
      triangles[4 * i + 3][0] = t0;
      triangles[4 * i + 3][1] = t1;
      triangles[4 * i + 3][2] = t3;
      // find which face is on the boundary
      for (int j = 0; j < 50; j++) {
        if (neighbor[t0][j] == -1) {
          neighbor[t0][j] = i;
          break;
        }
      }
      for (int j = 0; j < 50; j++) {
        if (neighbor[t1][j] == -1) {
          neighbor[t1][j] = i;
          break;
        }
      }
      for (int j = 0; j < 50; j++) {
        if (neighbor[t2][j] == -1) {
          neighbor[t2][j] = i;
          break;
        }
      }
      for (int j = 0; j < 50; j++) {
        if (neighbor[t3][j] == -1) {
          neighbor[t3][j] = i;
          break;
        }
      }
    }
    for (int i = 0; i < numTris / 4; i++) {
      for (int j = 0; j < 4; j++) {
        int vv0 = triangles[4 * i + j][0];
        int vv1 = triangles[4 * i + j][1];
        int vv2 = triangles[4 * i + j][2];
        bound_tri[4 * i + j] = 0;
        for (int k = 0; k < 50; k++) {
          int itet = neighbor[vv0][k];
          if (itet == -1) {
            break;
          }
          int t0, t1, t2, t3;
          t0 = triangles[4 * itet][0];
          t1 = triangles[4 * itet][1];
          t2 = triangles[4 * itet][2];
          t3 = triangles[4 * itet + 1][2];
          int bool_0 = 0;
          if (t0 == vv0 || t1 == vv0 || t2 == vv0 || t3 == vv0) {
            bool_0++;
          }
          if (t0 == vv1 || t1 == vv1 || t2 == vv1 || t3 == vv1) {
            bool_0++;
          }
          if (t0 == vv2 || t1 == vv2 || t2 == vv2 || t3 == vv2) {
            bool_0++;
          }
          if (bool_0 == 3) {
            bound_tri[4 * i + j]++;
          }
        }
      }
    }
    for (int i = 0; i < numVerts; i++) {
      delete[] neighbor[i];
    }
    delete[] neighbor;
    break;
  case QUAD: // read quad mesh
    if (fscanf(file, "%d %d\n", &numVerts, &numQuads) != 2)
      error("Couldn't read number of vertices and quads!");
    verts = (float(*)[3])reallocSafely(verts, sizeof(float[3]) * numVerts);
    normals = (float(*)[3])reallocSafely(normals, sizeof(float[3]) * numVerts);
    quads = (uint(*)[4])reallocSafely(quads, sizeof(uint[4]) * numQuads);
    bound_sign = (uint(*))reallocSafely(bound_sign, sizeof(uint) * numVerts);
    normals[0][0] = 9999;
    for (int i = 0; i < numVerts; i++) {
      bound_sign[i] = 1;
      switch (fileType) { readVertex(i, false, file); }
    }
    for (int i = 0; i < numQuads; i++) {
      if (fscanf(file, "%d %d %d %d\n", &quads[i][0], &quads[i][1],
                 &quads[i][2], &quads[i][3]) != 4)
        error("Error reading quad!");
    }
    break;
  case HEX: // read hex mesh
    if (fscanf(file, "%d %d\n", &numVerts, &numHexas) != 2)
      error("Couldn't read number of vertices and hexagons!");
    numQuads = numHexas * 6;
    verts = (float(*)[3])reallocSafely(verts, sizeof(float[3]) * numVerts);
    normals = (float(*)[3])reallocSafely(normals, sizeof(float[3]) * numVerts);
    quads = (uint(*)[4])reallocSafely(quads, sizeof(uint[4]) * numHexas * 6);
    colors = (float(*)[3])reallocSafely(colors, sizeof(float[3]) * numVerts);
    bound_sign = (uint(*))reallocSafely(bound_sign, sizeof(uint) * numVerts);
    for (int i = 0; i < numVerts; i++) {
      switch (fileType) { readVertex(i, true, file); }
      if (bound_sign[i] > 1) {
        r = 1;
        g = 0.59;
        b = 0.59;
        bound_sign[i] = 1;
      }
    }
    neighbor = new int *[numVerts];
    for (int i = 0; i < numVerts; i++) {
      neighbor[i] = 0;
    }
    for (int i = 0; i < numVerts; i++) {
      neighbor[i] = new int[50];
    }
    for (int i = 0; i < numVerts; i++) {
      for (int j = 0; j < 50; j++) {
        neighbor[i][j] = -1;
      }
    }
    for (int i = 0; i < numHexas; i++) {
      int t0, t1, t2, t3, t4, t5, t6, t7;
      if (fscanf(file, "%d %d %d %d %d %d %d %d\n", &t0, &t1, &t2, &t3, &t4,
                 &t5, &t6, &t7) != 8)
        error("Error reading hexagon!");
      quads[6 * i + 0][0] = t0;
      quads[6 * i + 0][1] = t3;
      quads[6 * i + 0][2] = t2;
      quads[6 * i + 0][3] = t1;
      quads[6 * i + 1][0] = t4;
      quads[6 * i + 1][1] = t5;
      quads[6 * i + 1][2] = t6;
      quads[6 * i + 1][3] = t7;
      quads[6 * i + 2][0] = t0;
      quads[6 * i + 2][1] = t4;
      quads[6 * i + 2][2] = t7;
      quads[6 * i + 2][3] = t3;
      quads[6 * i + 3][0] = t1;
      quads[6 * i + 3][1] = t2;
      quads[6 * i + 3][2] = t6;
      quads[6 * i + 3][3] = t5;
      quads[6 * i + 4][0] = t0;
      quads[6 * i + 4][1] = t1;
      quads[6 * i + 4][2] = t5;
      quads[6 * i + 4][3] = t4;
      quads[6 * i + 5][0] = t2;
      quads[6 * i + 5][1] = t3;
      quads[6 * i + 5][2] = t7;
      quads[6 * i + 5][3] = t6;
      for (int j = 0; j < 50; j++) {
        if (neighbor[t0][j] == -1) {
          neighbor[t0][j] = i;
          break;
        }
      }
      for (int j = 0; j < 50; j++) {
        if (neighbor[t1][j] == -1) {
          neighbor[t1][j] = i;
          break;
        }
      }
      for (int j = 0; j < 50; j++) {
        if (neighbor[t2][j] == -1) {
          neighbor[t2][j] = i;
          break;
        }
      }
      for (int j = 0; j < 50; j++) {
        if (neighbor[t3][j] == -1) {
          neighbor[t3][j] = i;
          break;
        }
      }
      for (int j = 0; j < 50; j++) {
        if (neighbor[t4][j] == -1) {
          neighbor[t4][j] = i;
          break;
        }
      }
      for (int j = 0; j < 50; j++) {
        if (neighbor[t5][j] == -1) {
          neighbor[t5][j] = i;
          break;
        }
      }
      for (int j = 0; j < 50; j++) {
        if (neighbor[t6][j] == -1) {
          neighbor[t6][j] = i;
          break;
        }
      }
      for (int j = 0; j < 50; j++) {
        if (neighbor[t7][j] == -1) {
          neighbor[t7][j] = i;
          break;
        }
      }
    }
    // find which face is on the boundary
    for (int i = 0; i < numHexas; i++) {
      for (int j = 0; j < 6; j++) {
        int vv0 = quads[6 * i + j][0];
        int vv1 = quads[6 * i + j][1];
        int vv2 = quads[6 * i + j][2];
        int vv3 = quads[6 * i + j][3];
        for (int k = 0; k < 50; k++) {
          int ihex = neighbor[vv0][k];
          if (ihex == -1) {
            break;
          }
          int t0 = quads[6 * ihex + 0][0];
          int t1 = quads[6 * ihex + 0][1];
          int t2 = quads[6 * ihex + 0][2];
          int t3 = quads[6 * ihex + 0][3];
          int t4 = quads[6 * ihex + 1][0];
          int t5 = quads[6 * ihex + 1][1];
          int t6 = quads[6 * ihex + 1][2];
          int t7 = quads[6 * ihex + 1][3];
          int bool_0 = 0;
          if (t0 == vv0 || t1 == vv0 || t2 == vv0 || t3 == vv0 || t4 == vv0 ||
              t5 == vv0 || t6 == vv0 || t7 == vv0) {
            bool_0++;
          }
          if (t0 == vv1 || t1 == vv1 || t2 == vv1 || t3 == vv1 || t4 == vv1 ||
              t5 == vv1 || t6 == vv1 || t7 == vv1) {
            bool_0++;
          }
          if (t0 == vv2 || t1 == vv2 || t2 == vv2 || t3 == vv2 || t4 == vv2 ||
              t5 == vv2 || t6 == vv2 || t7 == vv2) {
            bool_0++;
          }
          if (t0 == vv3 || t1 == vv3 || t2 == vv3 || t3 == vv3 || t4 == vv3 ||
              t5 == vv3 || t6 == vv3 || t7 == vv3) {
            bool_0++;
          }
        }
      }
    }
    for (int i = 0; i < numVerts; i++) {
      delete[] neighbor[i];
    }
    delete[] neighbor;
    break;
  case NURBS: // read nurbs mesh
    int nv1, npatch, *patch;
    npatch = 0;
    if (fscanf(file, "%d %d %d %d %d\n", &numVerts, &numHexas, &nv1, &numTris,
               &npatch) != 5)
      error("Couldn't parse NURBS file!");
    numVerts = numVerts + nv1;
    numQuads = numHexas * 6;
    verts = (float(*)[3])reallocSafely(verts, sizeof(float[3]) * numVerts);
    normals = (float(*)[3])reallocSafely(normals, sizeof(float[3]) * numVerts);
    colors = (float(*)[3])reallocSafely(colors, sizeof(float[3]) * numVerts);
    quads = (uint(*)[4])reallocSafely(quads, sizeof(uint[4]) * numHexas * 6);
    triangles = (uint(*)[3])reallocSafely(triangles, sizeof(uint[3]) * numTris);
    bound_sign = (uint(*))reallocSafely(bound_sign, sizeof(uint) * numVerts);
    bound_tri = (uint(*))reallocSafely(bound_tri, sizeof(uint) * numTris);
    patch = (int(*))reallocSafely(patch, sizeof(int) * npatch);
    for (int i = 0; i < numVerts; i++) {
      switch (fileType) { readVertex(i, true, file); }
    }
    for (int i = 0; i < npatch; i++) {
      int temp;
      if (fscanf(file, "%d\n", &temp) != 1)
        error("Error reading patch!");
      patch[i] = temp;
    }
    neighbor = new int *[numVerts];
    for (int i = 0; i < numVerts; i++) {
      neighbor[i] = 0;
    }
    for (int i = 0; i < numVerts; i++) {
      neighbor[i] = new int[50];
    }
    for (int i = 0; i < numVerts; i++) {
      for (int j = 0; j < 50; j++) {
        neighbor[i][j] = -1;
      }
    }

    for (int i = 0; i < numHexas; i++) {
      int t0, t1, t2, t3, t4, t5, t6, t7;
      if (fscanf(file, "%d %d %d %d %d %d %d %d\n", &t0, &t1, &t2, &t3, &t4,
                 &t5, &t6, &t7) != 8)
        error("Error reading hexagon!");
      quads[6 * i][0] = t0;
      quads[6 * i][1] = t3;
      quads[6 * i][2] = t2;
      quads[6 * i][3] = t1;
      quads[6 * i + 1][0] = t4;
      quads[6 * i + 1][1] = t5;
      quads[6 * i + 1][2] = t6;
      quads[6 * i + 1][3] = t7;
      quads[6 * i + 2][0] = t0;
      quads[6 * i + 2][1] = t4;
      quads[6 * i + 2][2] = t7;
      quads[6 * i + 2][3] = t3;
      quads[6 * i + 3][0] = t1;
      quads[6 * i + 3][1] = t2;
      quads[6 * i + 3][2] = t6;
      quads[6 * i + 3][3] = t5;
      quads[6 * i + 4][0] = t0;
      quads[6 * i + 4][1] = t1;
      quads[6 * i + 4][2] = t5;
      quads[6 * i + 4][3] = t4;
      quads[6 * i + 5][0] = t2;
      quads[6 * i + 5][1] = t3;
      quads[6 * i + 5][2] = t7;
      quads[6 * i + 5][3] = t6;
      for (int j = 0; j < 50; j++) {
        if (neighbor[t0][j] == -1) {
          neighbor[t0][j] = i;
          break;
        }
      }
      for (int j = 0; j < 50; j++) {
        if (neighbor[t1][j] == -1) {
          neighbor[t1][j] = i;
          break;
        }
      }
      for (int j = 0; j < 50; j++) {
        if (neighbor[t2][j] == -1) {
          neighbor[t2][j] = i;
          break;
        }
      }
      for (int j = 0; j < 50; j++) {
        if (neighbor[t3][j] == -1) {
          neighbor[t3][j] = i;
          break;
        }
      }
      for (int j = 0; j < 50; j++) {
        if (neighbor[t4][j] == -1) {
          neighbor[t4][j] = i;
          break;
        }
      }
      for (int j = 0; j < 50; j++) {
        if (neighbor[t5][j] == -1) {
          neighbor[t5][j] = i;
          break;
        }
      }
      for (int j = 0; j < 50; j++) {
        if (neighbor[t6][j] == -1) {
          neighbor[t6][j] = i;
          break;
        }
      }
      for (int j = 0; j < 50; j++) {
        if (neighbor[t7][j] == -1) {
          neighbor[t7][j] = i;
          break;
        }
      }
    }
    for (int i = 0; i < numHexas; i++) {
      for (int j = 0; j < 6; j++) {
        int vv0 = quads[6 * i + j][0];
        int vv1 = quads[6 * i + j][1];
        int vv2 = quads[6 * i + j][2];
        int vv3 = quads[6 * i + j][3];
        for (int k = 0; k < 50; k++) {
          int ihex = neighbor[vv0][k];
          if (ihex == -1) {
            break;
          }
          int t0 = quads[6 * ihex][0];
          int t1 = quads[6 * ihex][1];
          int t2 = quads[6 * ihex][2];
          int t3 = quads[6 * ihex][3];
          int t4 = quads[6 * ihex + 1][0];
          int t5 = quads[6 * ihex + 1][1];
          int t6 = quads[6 * ihex + 1][2];
          int t7 = quads[6 * ihex + 1][3];
          int bool_0 = 0;
          if (t0 == vv0 || t1 == vv0 || t2 == vv0 || t3 == vv0 || t4 == vv0 ||
              t5 == vv0 || t6 == vv0 || t7 == vv0) {
            bool_0++;
          }
          if (t0 == vv1 || t1 == vv1 || t2 == vv1 || t3 == vv1 || t4 == vv1 ||
              t5 == vv1 || t6 == vv1 || t7 == vv1) {
            bool_0++;
          }
          if (t0 == vv2 || t1 == vv2 || t2 == vv2 || t3 == vv2 || t4 == vv2 ||
              t5 == vv2 || t6 == vv2 || t7 == vv2) {
            bool_0++;
          }
          if (t0 == vv3 || t1 == vv3 || t2 == vv3 || t3 == vv3 || t4 == vv3 ||
              t5 == vv3 || t6 == vv3 || t7 == vv3) {
            bool_0++;
          }
        }
      }
    }
    for (int i = 0; i < numVerts; i++) {
      delete[] neighbor[i];
    }
    delete[] neighbor;
    for (int i = 0; i < numVerts; i++) {
      if (fscanf(file, "%f %f %f\n", &verts[i][0], &verts[i][1],
                 &verts[i][2]) != 3)
        error("Error reading vertex!");
      if (fileType != RAWC && fileType != RAWNC) {
        r = 1.0;
        g = 0.7;
        b = 0.0;
      }
      verts[i + numVerts][0] = verts[i][0];
      verts[i + numVerts][1] = verts[i][1];
      verts[i + numVerts][2] = verts[i][2];
      bound_sign[i + numVerts] = 1;
    }
    for (int i = 0; i < numTris; i++) {
      int t0, t1, t2;
      if (fscanf(file, "%d %d %d\n", &t0, &t1, &t2) != 3)
        error("Error reading triangle!");
      triangles[i][0] = t0 + numVerts;
      triangles[i][1] = t1 + numVerts;
      triangles[i][2] = t2 + numVerts;
      bound_tri[i] = 1;
    }
    break;
  }
  fclose(file);
}

void Geoframe::write_raw(const char *filename, MESHTYPE meshtype) {
  updateBySpan();
  if (meshtype == TRI || meshtype == TRI2) {
    saveTriangle(filename);
  } else if (meshtype == TET || meshtype == TET2) {
    saveTetra(filename);
  } else if (meshtype == HEX) {
    saveHexa(filename);
  } else if (meshtype == QUAD) {
    saveQuad(filename);
  }
}

void Geoframe::updateBySpan() {
  // arand: commented this unneeded output...
  // printf("%f %f %f\n",min_x,min_y,min_z);
  for (int i = 0; i < numVerts; i++) {
    verts[i][0] = min_x + span[0] * verts[i][0];
    verts[i][1] = min_y + span[1] * verts[i][1];
    verts[i][2] = min_z + span[2] * verts[i][2];
  }
}

void Geoframe::saveTriangle(const char *filename) {

  FILETYPE outFileType = INVALIDFILE;
  // arand: added else statements below... previously every file was being
  // treated as a RAW file
  if (strstr(filename, ".rawnc"))
    outFileType = RAWNC;
  else if (strstr(filename, ".rawn"))
    outFileType = RAWN;
  else if (strstr(filename, ".rawc"))
    outFileType = RAWC;
  else if (strstr(filename, ".raw"))
    outFileType = RAW;
  if (outFileType == INVALIDFILE) {
    cout << "ERROR: Incorrect file extension: " << filename << endl;
    exit(-1);
  }

  FILE *fp = fopen(filename, "w");
  int v0, v1, v2;
  float r, a, b, c, p, sum_area, center;
  sum_area = 0;
  center = 32;
  printf("number of vertices: %d\n", numVerts);
  printf("number of triangles: %d\n", numTris);
  fprintf(fp, "%d %d\n", numVerts, numTris);
  for (int i = 0; i < numVerts; i++) {
    r = (verts[i][0] - center) * (verts[i][0] - center) +
        (verts[i][1] - center) * (verts[i][1] - center) +
        (verts[i][2] - center) * (verts[i][2] - center);
    fprintf(fp, "%f %f %f", verts[i][0], verts[i][1], verts[i][2]);
    if (outFileType == RAWN || outFileType == RAWNC) {
      // print normals
      fprintf(fp, " %f %f %f", normals[i][0], normals[i][1], normals[i][2]);
    }
    if (outFileType == RAWC || outFileType == RAWNC) {
      // print normals
      fprintf(fp, " %f %f %f", colors[i][0], colors[i][1], colors[i][2]);
    }
    fprintf(fp, "\n");
  }
  for (int i = 0; i < numTris; i++) {
    v0 = triangles[i][0];
    v1 = triangles[i][1];
    v2 = triangles[i][2];
    r = (verts[v0][0] - center) * (verts[v0][0] - center) +
        (verts[v0][1] - center) * (verts[v0][1] - center) +
        (verts[v0][2] - center) * (verts[v0][2] - center);
    if (sqrt(r) < 17.0) {
      a = sqrt((verts[v1][0] - verts[v0][0]) * (verts[v1][0] - verts[v0][0]) +
               (verts[v1][1] - verts[v0][1]) * (verts[v1][1] - verts[v0][1]) +
               (verts[v1][2] - verts[v0][2]) * (verts[v1][2] - verts[v0][2]));
      b = sqrt((verts[v2][0] - verts[v1][0]) * (verts[v2][0] - verts[v1][0]) +
               (verts[v2][1] - verts[v1][1]) * (verts[v2][1] - verts[v1][1]) +
               (verts[v2][2] - verts[v1][2]) * (verts[v2][2] - verts[v1][2]));
      c = sqrt((verts[v0][0] - verts[v2][0]) * (verts[v0][0] - verts[v2][0]) +
               (verts[v0][1] - verts[v2][1]) * (verts[v0][1] - verts[v2][1]) +
               (verts[v0][2] - verts[v2][2]) * (verts[v0][2] - verts[v2][2]));
      p = (a + b + c) / 2.0;
      sum_area += sqrt(p * (p - a) * (p - b) * (p - c));
    }
    fprintf(fp, "%d %d %d\n", triangles[i][0], triangles[i][1],
            triangles[i][2]);
  }
  printf("finished\n");
  fclose(fp);
}

void Geoframe::saveTetra(const char *filename) {
  FILE *fp = fopen(filename, "w");
  float center = 64;
  fprintf(fp, "%d %d\n", numVerts, numTris / 4);
  for (int i = 0; i < numVerts; i++) {
    fprintf(fp, "%f %f %f\n", (verts[i][0] - center), (verts[i][1] - center),
            (verts[i][2] - center));
  }
  for (int i = 0; i < numTris / 4; i++) {
    fprintf(fp, "%d %d %d %d\n", triangles[4 * i][0], triangles[4 * i][1],
            triangles[4 * i][2], triangles[4 * i + 1][2]);
  }
  fclose(fp);
}

void Geoframe::saveHexa(const char *filename) {
  FILE *fp = fopen(filename, "w");
  fprintf(fp, "%d %d\n", numVerts, numQuads / 6);
  for (int i = 0; i < numVerts; i++) {
    fprintf(fp, "%f %f %f %d\n", verts[i][0], verts[i][1], verts[i][2],
            bound_sign[i]);
  }
  for (int i = 0; i < numQuads / 6; i++) {
    fprintf(fp, "%d %d %d %d %d %d %d %d\n", quads[6 * i][0], quads[6 * i][1],
            quads[6 * i][2], quads[6 * i][3], quads[6 * i + 1][1],
            quads[6 * i + 1][0], quads[6 * i + 1][3], quads[6 * i + 1][2]);
  }
  fclose(fp);
}

void Geoframe::saveQuad(const char *filename) {
  // Output 4-node quad element
  FILE *fp = fopen(filename, "w");
  fprintf(fp, "%d %d\n", numVerts, numQuads);
  for (int i = 0; i < numVerts; i++) {
    fprintf(fp, "%f %f %f\n", verts[i][0], verts[i][1], verts[i][2]);
  }
  for (int i = 0; i < numQuads; i++) {
    fprintf(fp, "%d %d %d %d\n", quads[i][0], quads[i][1], quads[i][2],
            quads[i][3]);
  }
  fclose(fp);
}

void Geoframe::Clear() {
  numTris = 0;
  numVerts = 0;
  numQuads = 0;
  numHexas = 0;
}

int Geoframe::getNumTris(void) { return numTris; }

int Geoframe::getNumQuads(void) { return numQuads; }

int Geoframe::getNumHexas(void) { return numHexas; }

int Geoframe::getNumVerts(void) { return numVerts; }

int Geoframe::TestNum(uint *v) {
  float dis_0 = 0.0;
  float dis_1 = 0.0;
  float dis_2 = 0.0;
  float dis_3 = 0.0;
  for (int i = 0; i < 3; i++) {
    dis_0 +=
        (verts[v[1]][i] - verts[v[0]][i]) * (verts[v[1]][i] - verts[v[0]][i]);
    dis_1 +=
        (verts[v[2]][i] - verts[v[1]][i]) * (verts[v[2]][i] - verts[v[1]][i]);
    dis_2 +=
        (verts[v[3]][i] - verts[v[2]][i]) * (verts[v[3]][i] - verts[v[2]][i]);
    dis_3 +=
        (verts[v[0]][i] - verts[v[3]][i]) * (verts[v[0]][i] - verts[v[3]][i]);
  }
  dis_0 = (float)sqrt(dis_0);
  dis_1 = (float)sqrt(dis_1);
  dis_2 = (float)sqrt(dis_2);
  dis_3 = (float)sqrt(dis_3);
  if (dis_0 == 0.0 || dis_1 == 0.0) {
    numTris++;
    return 3;
  } else if (dis_2 == 0.0 || dis_3 == 0.0) {
    numTris++;
    return 3;
  } else {
    return 4;
  }
}

int Geoframe::AddQuad(uint *v, int num) {
  assert(num == 3 || num == 4);
  num = 4;
  if (numQuads >= qsize) {
    qsize <<= 1;
    quads = (uint(*)[4])reallocSafely(quads, sizeof(uint[4]) * qsize);
  }
  if (num == 4) {
    quads[numQuads][0] = v[0];
    quads[numQuads][1] = v[1];
    quads[numQuads][2] = v[2];
    quads[numQuads][3] = v[3];
    return numQuads++;
  } else // (num==3)
  {
    triangles[numTris][0] = v[0];
    triangles[numTris][1] = v[1];
    triangles[numTris][2] = v[2];
    return numTris++;
  }
}

void Geoframe::AddQuad_indirect(uint *v) {
  float pv[3], pv0[3], pv1[3], pv2[3], pv3[3];
  float nv[3], nv0[3], nv1[3], nv2[3], nv3[3];
  uint v_new[5], v_quad[4];
  for (int i = 0; i < 3; i++) {
    pv0[i] = (verts[v[0]][i] + verts[v[1]][i]) / 2;
    pv1[i] = (verts[v[1]][i] + verts[v[2]][i]) / 2;
    pv2[i] = (verts[v[2]][i] + verts[v[3]][i]) / 2;
    pv3[i] = (verts[v[3]][i] + verts[v[0]][i]) / 2;
    nv0[i] = (normals[v[0]][i] + normals[v[1]][i]) / 2;
    nv1[i] = (normals[v[1]][i] + normals[v[2]][i]) / 2;
    nv2[i] = (normals[v[2]][i] + normals[v[3]][i]) / 2;
    nv3[i] = (normals[v[3]][i] + normals[v[0]][i]) / 2;
  }
  v_new[0] = AddVert(pv0, nv0);
  v_new[1] = AddVert(pv1, nv1);
  v_new[2] = AddVert(pv2, nv2);
  v_new[3] = AddVert(pv3, nv3);
  AddBound(v_new[0], 1);
  AddBound(v_new[1], 1);
  AddBound(v_new[2], 1);
  AddBound(v_new[3], 1);
  if (v[0] == v[1]) {
    for (int i = 0; i < 3; i++) {
      pv[i] = (verts[v[0]][i] + 2 * verts[v_new[2]][i]) / 3;
      nv[i] = (normals[v[0]][i] + 2 * normals[v_new[2]][i]) / 3;
    }
    v_new[4] = AddVert(pv, nv);
    AddBound(v_new[4], 1);
    v_quad[0] = v[0];
    v_quad[1] = v_new[1];
    v_quad[2] = v_new[4];
    v_quad[3] = v_new[3];
    AddQuad(v_quad, 4);
    v_quad[0] = v[2];
    v_quad[1] = v_new[2];
    v_quad[2] = v_new[4];
    v_quad[3] = v_new[1];
    AddQuad(v_quad, 4);
    v_quad[0] = v[3];
    v_quad[1] = v_new[3];
    v_quad[2] = v_new[4];
    v_quad[3] = v_new[2];
    AddQuad(v_quad, 4);
  } else if (v[1] == v[2]) {
    for (int i = 0; i < 3; i++) {
      pv[i] = (verts[v[1]][i] + 2 * verts[v_new[3]][i]) / 3;
      nv[i] = (normals[v[1]][i] + 2 * normals[v_new[3]][i]) / 3;
    }
    v_new[4] = AddVert(pv, nv);
    AddBound(v_new[4], 1);
    v_quad[0] = v[0];
    v_quad[1] = v_new[0];
    v_quad[2] = v_new[4];
    v_quad[3] = v_new[3];
    AddQuad(v_quad, 4);
    v_quad[0] = v[1];
    v_quad[1] = v_new[2];
    v_quad[2] = v_new[4];
    v_quad[3] = v_new[0];
    AddQuad(v_quad, 4);
    v_quad[0] = v[3];
    v_quad[1] = v_new[3];
    v_quad[2] = v_new[4];
    v_quad[3] = v_new[2];
    AddQuad(v_quad, 4);
  } else if (v[2] == v[3]) {
    for (int i = 0; i < 3; i++) {
      pv[i] = (verts[v[2]][i] + 2 * verts[v_new[0]][i]) / 3;
      nv[i] = (normals[v[2]][i] + 2 * normals[v_new[0]][i]) / 3;
    }
    v_new[4] = AddVert(pv, nv);
    AddBound(v_new[4], 1);
    v_quad[0] = v[0];
    v_quad[1] = v_new[0];
    v_quad[2] = v_new[4];
    v_quad[3] = v_new[3];
    AddQuad(v_quad, 4);
    v_quad[0] = v[1];
    v_quad[1] = v_new[1];
    v_quad[2] = v_new[4];
    v_quad[3] = v_new[0];
    AddQuad(v_quad, 4);
    v_quad[0] = v[2];
    v_quad[1] = v_new[3];
    v_quad[2] = v_new[4];
    v_quad[3] = v_new[1];
    AddQuad(v_quad, 4);
  } else if (v[3] == v[0]) {
    for (int i = 0; i < 3; i++) {
      pv[i] = (verts[v[0]][i] + 2 * verts[v_new[1]][i]) / 3;
      nv[i] = (normals[v[0]][i] + 2 * normals[v_new[1]][i]) / 3;
    }
    v_new[4] = AddVert(pv, nv);
    AddBound(v_new[4], 1);
    v_quad[0] = v[0];
    v_quad[1] = v_new[0];
    v_quad[2] = v_new[4];
    v_quad[3] = v_new[2];
    AddQuad(v_quad, 4);
    v_quad[0] = v[1];
    v_quad[1] = v_new[1];
    v_quad[2] = v_new[4];
    v_quad[3] = v_new[0];
    AddQuad(v_quad, 4);
    v_quad[0] = v[2];
    v_quad[1] = v_new[2];
    v_quad[2] = v_new[4];
    v_quad[3] = v_new[1];
    AddQuad(v_quad, 4);
  } else {
    for (int i = 0; i < 3; i++) {
      pv[i] =
          (verts[v[0]][i] + verts[v[1]][i] + verts[v[2]][i] + verts[v[3]][i]) /
          4;
      nv[i] = (normals[v[0]][i] + normals[v[1]][i] + normals[v[2]][i] +
               normals[v[3]][i]) /
              4;
    }
    v_new[4] = AddVert(pv, nv);
    AddBound(v_new[4], 1);
    v_quad[0] = v[0];
    v_quad[1] = v_new[0];
    v_quad[2] = v_new[4];
    v_quad[3] = v_new[3];
    AddQuad(v_quad, 4);
    v_quad[0] = v[1];
    v_quad[1] = v_new[1];
    v_quad[2] = v_new[4];
    v_quad[3] = v_new[0];
    AddQuad(v_quad, 4);
    v_quad[0] = v[2];
    v_quad[1] = v_new[2];
    v_quad[2] = v_new[4];
    v_quad[3] = v_new[1];
    AddQuad(v_quad, 4);
    v_quad[0] = v[3];
    v_quad[1] = v_new[3];
    v_quad[2] = v_new[4];
    v_quad[3] = v_new[2];
    AddQuad(v_quad, 4);
  }
}

uint Geoframe::AddRefine_edgevtx(uint index_1, uint index_2) {
  uint v_new;
  float pv[3], nv[3];
  for (int i = 0; i < 18; i++) {
    if (refine_edge[index_1][i] == (int)index_2) {
      return refine_edgevtx[index_1][i];
    }
    if (refine_edge[index_1][i] == -1) {
      refine_edge[index_1][i] = index_2;
      for (int j = 0; j < 3; j++) {
        pv[j] = (2 * verts[index_1][j] + verts[index_2][j]) / 3;
        nv[j] = (2 * normals[index_1][j] + normals[index_2][j]) / 3;
      }
      v_new = AddVert(pv, nv);
      AddBound(v_new, 1);
      refine_edgevtx[index_1][i] = v_new;
      return v_new;
    }
  }
  return 0;
}

uint Geoframe::AddRefine_edgevtx_hex(uint index_1, uint index_2) {
  uint v_new;
  float pv[3], nv[3];
  for (int i = 0; i < 18; i++) {
    if (refine_edge[index_1][i] == (int)index_2) {
      return refine_edgevtx[index_1][i];
    }
    if (refine_edge[index_1][i] == -1) {
      refine_edge[index_1][i] = index_2;
      for (int j = 0; j < 3; j++) {
        pv[j] = (2 * verts[index_1][j] + verts[index_2][j]) / 3;
        nv[j] = (2 * normals[index_1][j] + normals[index_2][j]) / 3;
      }
      v_new = AddVert(pv, nv);
      refine_edgevtx[index_1][i] = v_new;
      return v_new;
    }
  }
  return 0;
}

uint Geoframe::AddRefine_facevtx(uint index_1, uint index_2, uint index_3,
                                 uint index_4) {
  uint v_new;
  float pv[3], nv[3];
  for (int i = 0; i < 18; i++) {
    if (refine_edge[index_1][i] == (int)index_2) {
      return refine_edgevtx[index_1][i];
    }
    if (refine_edge[index_1][i] == -1) {
      refine_edge[index_1][i] = index_2;
      for (int j = 0; j < 3; j++) {
        pv[j] = (4 * verts[index_1][j] + 2 * verts[index_3][j] +
                 2 * verts[index_4][j] + verts[index_2][j]) /
                9;
        nv[j] = (4 * normals[index_1][j] + 2 * normals[index_3][j] +
                 2 * normals[index_4][j] + normals[index_2][j]) /
                9;
      }
      v_new = AddVert(pv, nv);
      refine_edgevtx[index_1][i] = v_new;
      return v_new;
    }
  }
  return 0;
}

void Geoframe::AddVert_adaptive(uint *v, uint *v_new) {
  float pv[3], pv0[3], pv1[3], pv2[3], pv3[3];
  float nv[3], nv0[3], nv1[3], nv2[3], nv3[3];
  for (int i = 0; i < 3; i++) {
    pv[i] =
        (verts[v[0]][i] + verts[v[1]][i] + verts[v[2]][i] + verts[v[3]][i]) / 4;
    pv0[i] = (verts[v[0]][i] + 2 * pv[i]) / 3;
    pv1[i] = (verts[v[1]][i] + 2 * pv[i]) / 3;
    pv2[i] = (verts[v[2]][i] + 2 * pv[i]) / 3;
    pv3[i] = (verts[v[3]][i] + 2 * pv[i]) / 3;
    nv[i] = (normals[v[0]][i] + normals[v[1]][i] + normals[v[2]][i] +
             normals[v[3]][i]) /
            4;
    nv0[i] = (normals[v[0]][i] + 2 * nv[i]) / 3;
    nv1[i] = (normals[v[1]][i] + 2 * nv[i]) / 3;
    nv2[i] = (normals[v[2]][i] + 2 * nv[i]) / 3;
    nv3[i] = (normals[v[3]][i] + 2 * nv[i]) / 3;
  }
  v_new[0] = AddVert(pv0, nv0);
  v_new[1] = AddVert(pv1, nv1);
  v_new[2] = AddVert(pv2, nv2);
  v_new[3] = AddVert(pv3, nv3);
  AddBound(v_new[0], 1);
  AddBound(v_new[1], 1);
  AddBound(v_new[2], 1);
  AddBound(v_new[3], 1);
}

void Geoframe::AddQuad_adaptive(uint *v, uint *v_new, int num) {
  uint vv[4], v_quad[4];
  for (int i = 0; i < 4; i++) {
    vv[i] = v_new[i];
  }
  v_quad[0] = v[0];
  v_quad[1] = v[1];
  v_quad[2] = vv[1];
  v_quad[3] = vv[0];
  AddQuad(v_quad, num);
  v_quad[0] = v[1];
  v_quad[1] = v[2];
  v_quad[2] = vv[2];
  v_quad[3] = vv[1];
  AddQuad(v_quad, num);
  v_quad[0] = v[2];
  v_quad[1] = v[3];
  v_quad[2] = vv[3];
  v_quad[3] = vv[2];
  AddQuad(v_quad, num);
  v_quad[0] = v[3];
  v_quad[1] = v[0];
  v_quad[2] = vv[0];
  v_quad[3] = vv[3];
  AddQuad(v_quad, num);
  v_quad[0] = vv[0];
  v_quad[1] = vv[1];
  v_quad[2] = vv[2];
  v_quad[3] = vv[3];
  AddQuad(v_quad, num);
}

void Geoframe::AddVert_adaptive_2_1(uint *v, uint *v_new) {
  float pv0[3], pv1[3], pv2[3], pv3[3];
  float nv0[3], nv1[3], nv2[3], nv3[3];
  for (int i = 0; i < 3; i++) {
    pv0[i] = (2 * verts[v[0]][i] + verts[v[1]][i]) / 3;
    pv1[i] = (2 * verts[v[1]][i] + verts[v[0]][i]) / 3;
    pv2[i] = (pv1[i] + 2 * (2 * verts[v[2]][i] + verts[v[3]][i]) / 3) / 3;
    pv3[i] = (2 * pv0[i] + (2 * verts[v[3]][i] + verts[v[2]][i]) / 3) / 3;
    nv0[i] = (2 * normals[v[0]][i] + normals[v[1]][i]) / 3;
    nv1[i] = (2 * normals[v[1]][i] + normals[v[0]][i]) / 3;
    nv2[i] = (nv1[i] + 2 * (2 * normals[v[2]][i] + normals[v[3]][i]) / 3) / 3;
    nv3[i] = (2 * nv0[i] + (2 * normals[v[3]][i] + normals[v[2]][i]) / 3) / 3;
  }
  v_new[2] = AddVert(pv2, nv2);
  v_new[3] = AddVert(pv3, nv3);
  AddBound(v_new[2], 1);
  AddBound(v_new[3], 1);
  v_new[0] = AddRefine_edgevtx(v[0], v[1]);
  v_new[1] = AddRefine_edgevtx(v[1], v[0]);
  v_new[4] = AddRefine_edgevtx(v[0], v[3]);
  v_new[5] = AddRefine_edgevtx(v[3], v[0]);
}

void Geoframe::AddQuad_adaptive_2_1(uint *v, uint *v_new, int num) {
  uint vv[6], v_quad[4];
  for (int i = 0; i < 6; i++) {
    vv[i] = v_new[i];
  }
  v_quad[0] = vv[0];
  v_quad[1] = vv[1];
  v_quad[2] = vv[2];
  v_quad[3] = vv[3];
  AddQuad(v_quad, num);
  v_quad[0] = v[1];
  v_quad[1] = v[2];
  v_quad[2] = vv[2];
  v_quad[3] = vv[1];
  AddQuad(v_quad, num);
  v_quad[0] = v[2];
  v_quad[1] = v[3];
  v_quad[2] = vv[5];
  v_quad[3] = vv[2];
  AddQuad(v_quad, num);
  v_quad[0] = v[0];
  v_quad[1] = vv[0];
  v_quad[2] = vv[3];
  v_quad[3] = vv[4];
  AddQuad(v_quad, num);
  v_quad[0] = vv[2];
  v_quad[1] = vv[5];
  v_quad[2] = vv[4];
  v_quad[3] = vv[3];
  AddQuad(v_quad, num);
}

void Geoframe::AddVert_adaptive_2_3(uint *v, uint *v_new) {
  float pv0[3], pv1[3], pv6[3], pv7[3];
  float nv0[3], nv1[3], nv6[3], nv7[3];
  for (int i = 0; i < 3; i++) {
    pv0[i] = (2 * verts[v[0]][i] + verts[v[1]][i]) / 3;
    pv1[i] = (2 * verts[v[1]][i] + verts[v[0]][i]) / 3;
    pv6[i] = (5 * pv1[i] + (2 * verts[v[2]][i] + verts[v[3]][i]) / 3) / 6;
    pv7[i] = (5 * pv0[i] + (2 * verts[v[3]][i] + verts[v[2]][i]) / 3) / 6;
    nv0[i] = (2 * normals[v[0]][i] + normals[v[1]][i]) / 3;
    nv1[i] = (2 * normals[v[1]][i] + normals[v[0]][i]) / 3;
    nv6[i] = (5 * nv1[i] + (2 * normals[v[2]][i] + normals[v[3]][i]) / 3) / 6;
    nv7[i] = (5 * nv0[i] + (2 * normals[v[3]][i] + normals[v[2]][i]) / 3) / 6;
  }
  v_new[6] = AddVert(pv6, nv6);
  v_new[7] = AddVert(pv7, nv7);
  AddBound(v_new[6], 1);
  AddBound(v_new[7], 1);
  v_new[0] = AddRefine_edgevtx(v[0], v[1]);
  v_new[1] = AddRefine_edgevtx(v[1], v[0]);
  v_new[2] = AddRefine_edgevtx(v[1], v[2]);
  v_new[3] = AddRefine_edgevtx(v[2], v[1]);
  v_new[4] = AddRefine_edgevtx(v[0], v[3]);
  v_new[5] = AddRefine_edgevtx(v[3], v[0]);
}

void Geoframe::AddQuad_adaptive_2_3(uint *v, uint *v_new, int num) {
  uint vv[8], v_quad[4];
  for (int i = 0; i < 8; i++) {
    vv[i] = v_new[i];
  }
  v_quad[0] = vv[0];
  v_quad[1] = vv[1];
  v_quad[2] = vv[6];
  v_quad[3] = vv[7];
  AddQuad(v_quad, num);
  v_quad[0] = v[1];
  v_quad[1] = vv[2];
  v_quad[2] = vv[6];
  v_quad[3] = vv[1];
  AddQuad(v_quad, num);
  v_quad[0] = v[0];
  v_quad[1] = vv[0];
  v_quad[2] = vv[7];
  v_quad[3] = vv[4];
  AddQuad(v_quad, num);
  v_quad[0] = vv[6];
  v_quad[1] = vv[2];
  v_quad[2] = vv[4];
  v_quad[3] = vv[7];
  AddQuad(v_quad, num);
  v_quad[0] = vv[2];
  v_quad[1] = vv[3];
  v_quad[2] = vv[5];
  v_quad[3] = vv[4];
  AddQuad(v_quad, num);
  v_quad[0] = v[2];
  v_quad[1] = v[3];
  v_quad[2] = vv[5];
  v_quad[3] = vv[3];
  AddQuad(v_quad, num);
}

void Geoframe::AddVert_adaptive_4(uint *v, uint *v_new) {
  float pv0[3], pv1[3], pv2[3], pv3[3], pv4[3], pv5[3], pv6[3], pv7[3], pv8[3],
      pv9[3], pv10[3], pv11[3];
  float nv0[3], nv1[3], nv2[3], nv3[3], nv4[3], nv5[3], nv6[3], nv7[3], nv8[3],
      nv9[3], nv10[3], nv11[3];
  for (int i = 0; i < 3; i++) {
    pv0[i] = (2 * verts[v[0]][i] + verts[v[1]][i]) / 3;
    pv1[i] = (2 * verts[v[1]][i] + verts[v[0]][i]) / 3;
    pv2[i] = (2 * verts[v[1]][i] + verts[v[2]][i]) / 3;
    pv3[i] = (2 * verts[v[2]][i] + verts[v[1]][i]) / 3;
    pv4[i] = (2 * verts[v[2]][i] + verts[v[3]][i]) / 3;
    pv5[i] = (2 * verts[v[3]][i] + verts[v[2]][i]) / 3;
    pv6[i] = (2 * verts[v[0]][i] + verts[v[3]][i]) / 3;
    pv7[i] = (2 * verts[v[3]][i] + verts[v[0]][i]) / 3;
    pv8[i] = (2 * pv0[i] + pv5[i]) / 3;
    pv9[i] = (2 * pv5[i] + pv0[i]) / 3;
    pv10[i] = (2 * pv1[i] + pv4[i]) / 3;
    pv11[i] = (2 * pv4[i] + pv1[i]) / 3;
    nv0[i] = (2 * normals[v[0]][i] + normals[v[1]][i]) / 3;
    nv1[i] = (2 * normals[v[1]][i] + normals[v[0]][i]) / 3;
    nv2[i] = (2 * normals[v[1]][i] + normals[v[2]][i]) / 3;
    nv3[i] = (2 * normals[v[2]][i] + normals[v[1]][i]) / 3;
    nv4[i] = (2 * normals[v[2]][i] + normals[v[3]][i]) / 3;
    nv5[i] = (2 * normals[v[3]][i] + normals[v[2]][i]) / 3;
    nv6[i] = (2 * normals[v[0]][i] + normals[v[3]][i]) / 3;
    nv7[i] = (2 * normals[v[3]][i] + normals[v[0]][i]) / 3;
    nv8[i] = (2 * nv0[i] + nv5[i]) / 3;
    nv9[i] = (2 * nv5[i] + nv0[i]) / 3;
    nv10[i] = (2 * nv1[i] + nv4[i]) / 3;
    nv11[i] = (2 * nv4[i] + nv1[i]) / 3;
  }
  v_new[8] = AddVert(pv8, nv8);
  v_new[9] = AddVert(pv9, nv9);
  v_new[10] = AddVert(pv10, nv10);
  v_new[11] = AddVert(pv11, nv11);
  AddBound(v_new[8], 1);
  AddBound(v_new[9], 1);
  AddBound(v_new[10], 1);
  AddBound(v_new[11], 1);
  v_new[0] = AddRefine_edgevtx(v[0], v[1]);
  v_new[1] = AddRefine_edgevtx(v[1], v[0]);
  v_new[2] = AddRefine_edgevtx(v[1], v[2]);
  v_new[3] = AddRefine_edgevtx(v[2], v[1]);
  v_new[4] = AddRefine_edgevtx(v[2], v[3]);
  v_new[5] = AddRefine_edgevtx(v[3], v[2]);
  v_new[6] = AddRefine_edgevtx(v[0], v[3]);
  v_new[7] = AddRefine_edgevtx(v[3], v[0]);
}

void Geoframe::AddQuad_adaptive_4(uint *v, uint *v_new, int num) {
  uint vv[12], v_quad[4];
  for (int i = 0; i < 12; i++) {
    vv[i] = v_new[i];
  }
  v_quad[0] = v[0];
  v_quad[1] = vv[0];
  v_quad[2] = vv[8];
  v_quad[3] = vv[6];
  AddQuad(v_quad, num);
  v_quad[0] = v[1];
  v_quad[1] = vv[2];
  v_quad[2] = vv[10];
  v_quad[3] = vv[1];
  AddQuad(v_quad, num);
  v_quad[0] = v[2];
  v_quad[1] = vv[4];
  v_quad[2] = vv[11];
  v_quad[3] = vv[3];
  AddQuad(v_quad, num);
  v_quad[0] = v[3];
  v_quad[1] = vv[7];
  v_quad[2] = vv[9];
  v_quad[3] = vv[5];
  AddQuad(v_quad, num);
  v_quad[0] = vv[0];
  v_quad[1] = vv[1];
  v_quad[2] = vv[10];
  v_quad[3] = vv[8];
  AddQuad(v_quad, num);
  v_quad[0] = vv[2];
  v_quad[1] = vv[3];
  v_quad[2] = vv[11];
  v_quad[3] = vv[10];
  AddQuad(v_quad, num);
  v_quad[0] = vv[4];
  v_quad[1] = vv[5];
  v_quad[2] = vv[9];
  v_quad[3] = vv[11];
  AddQuad(v_quad, num);
  v_quad[0] = vv[6];
  v_quad[1] = vv[8];
  v_quad[2] = vv[9];
  v_quad[3] = vv[7];
  AddQuad(v_quad, num);
  v_quad[0] = vv[8];
  v_quad[1] = vv[10];
  v_quad[2] = vv[11];
  v_quad[3] = vv[9];
  AddQuad(v_quad, num);
}

void Geoframe::AddVert_adaptive_3_1(uint *v, uint *v_new) {
  float pv0[3], pv1[3], nv0[3], nv1[3];
  v_new[0] = AddRefine_edgevtx(v[0], v[1]);
  v_new[2] = AddRefine_edgevtx(v[0], v[3]);
  for (int i = 0; i < 3; i++) {
    pv0[i] = (2 * verts[v[0]][i] + verts[v[1]][i]) / 3;
    pv1[i] = (2 * pv0[i] + (2 * verts[v[3]][i] + verts[v[2]][i]) / 3) / 3;
    nv0[i] = (2 * normals[v[0]][i] + normals[v[1]][i]) / 3;
    nv1[i] = (2 * nv0[i] + (2 * normals[v[3]][i] + normals[v[2]][i]) / 3) / 3;
  }
  v_new[1] = AddVert(pv1, nv1);
  AddBound(v_new[1], 1);
}

void Geoframe::AddQuad_adaptive_3_1(uint *v, uint *v_new, int num) {
  uint vv[3], v_quad[4];
  for (int i = 0; i < 3; i++) {
    vv[i] = v_new[i];
  }
  v_quad[0] = v[1];
  v_quad[1] = v[2];
  v_quad[2] = vv[1];
  v_quad[3] = vv[0];
  AddQuad(v_quad, num);
  v_quad[0] = v[2];
  v_quad[1] = v[3];
  v_quad[2] = vv[2];
  v_quad[3] = vv[1];
  AddQuad(v_quad, num);
  v_quad[0] = v[0];
  v_quad[1] = vv[0];
  v_quad[2] = vv[1];
  v_quad[3] = vv[2];
  AddQuad(v_quad, num);
}

void Geoframe::AddVert_adaptive_3_2a(uint *v, uint *v_new) {
  float pv0[3], pv1[3], pv2[3], pv3[3], pv4[3], pv5[3], pv6[3], pv7[3], pv8[3],
      pv9[3], pv10[3], pv11[3];
  float nv0[3], nv1[3], nv2[3], nv3[3], nv4[3], nv5[3], nv6[3], nv7[3], nv8[3],
      nv9[3], nv10[3], nv11[3];
  for (int i = 0; i < 3; i++) {
    pv0[i] = (2 * verts[v[0]][i] + verts[v[1]][i]) / 3;
    pv1[i] = (2 * verts[v[1]][i] + verts[v[0]][i]) / 3;
    pv2[i] = (2 * verts[v[1]][i] + verts[v[2]][i]) / 3;
    pv3[i] = (2 * verts[v[2]][i] + verts[v[1]][i]) / 3;
    pv4[i] = (2 * verts[v[2]][i] + verts[v[3]][i]) / 3;
    pv5[i] = (2 * verts[v[3]][i] + verts[v[2]][i]) / 3;
    pv6[i] = (2 * verts[v[0]][i] + verts[v[3]][i]) / 3;
    pv7[i] = (2 * verts[v[3]][i] + verts[v[0]][i]) / 3;
    pv8[i] = (2 * pv0[i] + pv5[i]) / 3;
    pv9[i] = (2 * pv5[i] + pv0[i]) / 3;
    pv10[i] = (2 * pv1[i] + pv4[i]) / 3;
    pv11[i] = (2 * pv4[i] + pv1[i]) / 3;
    nv0[i] = (2 * normals[v[0]][i] + normals[v[1]][i]) / 3;
    nv1[i] = (2 * normals[v[1]][i] + normals[v[0]][i]) / 3;
    nv2[i] = (2 * normals[v[1]][i] + normals[v[2]][i]) / 3;
    nv3[i] = (2 * normals[v[2]][i] + normals[v[1]][i]) / 3;
    nv4[i] = (2 * normals[v[2]][i] + normals[v[3]][i]) / 3;
    nv5[i] = (2 * normals[v[3]][i] + normals[v[2]][i]) / 3;
    nv6[i] = (2 * normals[v[0]][i] + normals[v[3]][i]) / 3;
    nv7[i] = (2 * normals[v[3]][i] + normals[v[0]][i]) / 3;
    nv8[i] = (2 * nv0[i] + nv5[i]) / 3;
    nv9[i] = (2 * nv5[i] + nv0[i]) / 3;
    nv10[i] = (2 * nv1[i] + nv4[i]) / 3;
    nv11[i] = (2 * nv4[i] + nv1[i]) / 3;
  }
  v_new[3] = AddVert(pv10, nv10);
  v_new[4] = AddVert(pv8, nv8);
  v_new[6] = AddVert(pv11, nv11);
  v_new[7] = AddVert(pv9, nv9);
  AddBound(v_new[3], 1);
  AddBound(v_new[4], 1);
  AddBound(v_new[6], 1);
  AddBound(v_new[7], 1);
  v_new[0] = AddRefine_edgevtx(v[0], v[1]);
  v_new[1] = AddRefine_edgevtx(v[1], v[0]);
  v_new[2] = AddRefine_edgevtx(v[1], v[2]);
  v_new[5] = AddRefine_edgevtx(v[0], v[3]);
}

void Geoframe::AddQuad_adaptive_3_2a(uint *v, uint *v_new, int num) {
  uint vv[8], v_quad[4];
  for (int i = 0; i < 8; i++) {
    vv[i] = v_new[i];
  }
  v_quad[0] = v[0];
  v_quad[1] = vv[0];
  v_quad[2] = vv[4];
  v_quad[3] = vv[5];
  AddQuad(v_quad, num);
  v_quad[0] = v[1];
  v_quad[1] = vv[2];
  v_quad[2] = vv[3];
  v_quad[3] = vv[1];
  AddQuad(v_quad, num);
  v_quad[0] = v[2];
  v_quad[1] = vv[6];
  v_quad[2] = vv[3];
  v_quad[3] = vv[2];
  AddQuad(v_quad, num);
  v_quad[0] = v[3];
  v_quad[1] = vv[5];
  v_quad[2] = vv[4];
  v_quad[3] = vv[7];
  AddQuad(v_quad, num);
  v_quad[0] = vv[0];
  v_quad[1] = vv[1];
  v_quad[2] = vv[3];
  v_quad[3] = vv[4];
  AddQuad(v_quad, num);
  v_quad[0] = vv[3];
  v_quad[1] = vv[6];
  v_quad[2] = vv[7];
  v_quad[3] = vv[4];
  AddQuad(v_quad, num);
  v_quad[0] = vv[7];
  v_quad[1] = vv[6];
  v_quad[2] = v[2];
  v_quad[3] = v[3];
  AddQuad(v_quad, num);
}

void Geoframe::AddVert_adaptive_3_2b(uint *v, uint *v_new) {
  float pv4[3];
  float nv4[3];
  for (int i = 0; i < 3; i++) {
    pv4[i] =
        (verts[v[0]][i] + verts[v[1]][i] + verts[v[2]][i] + verts[v[3]][i]) / 4;
    nv4[i] = (normals[v[0]][i] + normals[v[1]][i] + normals[v[2]][i] +
              normals[v[3]][i]) /
             43;
  }
  v_new[4] = AddVert(pv4, nv4);
  AddBound(v_new[4], 1);
  v_new[0] = AddRefine_edgevtx(v[0], v[1]);
  v_new[1] = AddRefine_edgevtx(v[2], v[1]);
  v_new[2] = AddRefine_edgevtx(v[2], v[3]);
  v_new[3] = AddRefine_edgevtx(v[0], v[3]);
}

void Geoframe::AddQuad_adaptive_3_2b(uint *v, uint *v_new, int num) {
  uint vv[5], v_quad[4];
  for (int i = 0; i < 5; i++) {
    vv[i] = v_new[i];
  }
  v_quad[0] = v[0];
  v_quad[1] = vv[0];
  v_quad[2] = vv[4];
  v_quad[3] = vv[3];
  AddQuad(v_quad, num);
  v_quad[0] = v[1];
  v_quad[1] = vv[1];
  v_quad[2] = vv[4];
  v_quad[3] = vv[0];
  AddQuad(v_quad, num);
  v_quad[0] = v[2];
  v_quad[1] = vv[2];
  v_quad[2] = vv[4];
  v_quad[3] = vv[1];
  AddQuad(v_quad, num);
  v_quad[0] = v[3];
  v_quad[1] = vv[3];
  v_quad[2] = vv[4];
  v_quad[3] = vv[2];
  AddQuad(v_quad, num);
}

void Geoframe::AddVert_adaptive_3_3(uint *v, uint *v_new) {
  float pv0[3], pv1[3], pv2[3], pv3[3], pv4[3], pv5[3], pv6[3], pv7[3], pv8[3],
      pv9[3], pv10[3], pv11[3];
  float nv0[3], nv1[3], nv2[3], nv3[3], nv4[3], nv5[3], nv6[3], nv7[3], nv8[3],
      nv9[3], nv10[3], nv11[3];
  for (int i = 0; i < 3; i++) {
    pv0[i] = (2 * verts[v[0]][i] + verts[v[1]][i]) / 3;
    pv1[i] = (2 * verts[v[1]][i] + verts[v[0]][i]) / 3;
    pv2[i] = (2 * verts[v[1]][i] + verts[v[2]][i]) / 3;
    pv3[i] = (2 * verts[v[2]][i] + verts[v[1]][i]) / 3;
    pv4[i] = (2 * verts[v[2]][i] + verts[v[3]][i]) / 3;
    pv5[i] = (2 * verts[v[3]][i] + verts[v[2]][i]) / 3;
    pv6[i] = (2 * verts[v[0]][i] + verts[v[3]][i]) / 3;
    pv7[i] = (2 * verts[v[3]][i] + verts[v[0]][i]) / 3;
    pv8[i] = (2 * pv0[i] + pv5[i]) / 3;
    pv9[i] = (2 * pv5[i] + pv0[i]) / 3;
    pv10[i] = (2 * pv1[i] + pv4[i]) / 3;
    pv11[i] = (2 * pv4[i] + pv1[i]) / 3;
    nv0[i] = (2 * normals[v[0]][i] + normals[v[1]][i]) / 3;
    nv1[i] = (2 * normals[v[1]][i] + normals[v[0]][i]) / 3;
    nv2[i] = (2 * normals[v[1]][i] + normals[v[2]][i]) / 3;
    nv3[i] = (2 * normals[v[2]][i] + normals[v[1]][i]) / 3;
    nv4[i] = (2 * normals[v[2]][i] + normals[v[3]][i]) / 3;
    nv5[i] = (2 * normals[v[3]][i] + normals[v[2]][i]) / 3;
    nv6[i] = (2 * normals[v[0]][i] + normals[v[3]][i]) / 3;
    nv7[i] = (2 * normals[v[3]][i] + normals[v[0]][i]) / 3;
    nv8[i] = (2 * nv0[i] + nv5[i]) / 3;
    nv9[i] = (2 * nv5[i] + nv0[i]) / 3;
    nv10[i] = (2 * nv1[i] + nv4[i]) / 3;
    nv11[i] = (2 * nv4[i] + nv1[i]) / 3;
  }
  v_new[5] = AddVert(pv11, nv11);
  v_new[7] = AddVert(pv10, nv10);
  v_new[8] = AddVert(pv8, nv8);
  v_new[9] = AddVert(pv9, nv9);
  AddBound(v_new[5], 1);
  AddBound(v_new[7], 1);
  AddBound(v_new[8], 1);
  AddBound(v_new[9], 1);
  v_new[0] = AddRefine_edgevtx(v[0], v[1]);
  v_new[1] = AddRefine_edgevtx(v[1], v[0]);
  v_new[2] = AddRefine_edgevtx(v[1], v[2]);
  v_new[3] = AddRefine_edgevtx(v[2], v[1]);
  v_new[4] = AddRefine_edgevtx(v[2], v[3]);
  v_new[6] = AddRefine_edgevtx(v[0], v[3]);
}

void Geoframe::AddQuad_adaptive_3_3(uint *v, uint *v_new, int num) {
  uint vv[10], v_quad[4];
  for (int i = 0; i < 10; i++) {
    vv[i] = v_new[i];
  }
  v_quad[0] = v[0];
  v_quad[1] = vv[0];
  v_quad[2] = vv[8];
  v_quad[3] = vv[6];
  AddQuad(v_quad, num);
  v_quad[0] = v[1];
  v_quad[1] = vv[2];
  v_quad[2] = vv[7];
  v_quad[3] = vv[1];
  AddQuad(v_quad, num);
  v_quad[0] = v[2];
  v_quad[1] = vv[4];
  v_quad[2] = vv[5];
  v_quad[3] = vv[3];
  AddQuad(v_quad, num);
  v_quad[0] = v[3];
  v_quad[1] = vv[6];
  v_quad[2] = vv[8];
  v_quad[3] = vv[9];
  AddQuad(v_quad, num);
  v_quad[0] = vv[0];
  v_quad[1] = vv[1];
  v_quad[2] = vv[7];
  v_quad[3] = vv[8];
  AddQuad(v_quad, num);
  v_quad[0] = vv[2];
  v_quad[1] = vv[3];
  v_quad[2] = vv[5];
  v_quad[3] = vv[7];
  AddQuad(v_quad, num);
  v_quad[0] = vv[4];
  v_quad[1] = v[3];
  v_quad[2] = vv[9];
  v_quad[3] = vv[5];
  AddQuad(v_quad, num);
  v_quad[0] = vv[9];
  v_quad[1] = vv[8];
  v_quad[2] = vv[7];
  v_quad[3] = vv[5];
  AddQuad(v_quad, num);
}

void Geoframe::AddVert_adaptive_4_2b(uint *v, uint *v_new) {
  float pv0[3], pv1[3], pv2[3], pv3[3], pv4[3], pv5[3], pv6[3], pv7[3], pv8[3],
      pv9[3], pv10[3], pv11[3];
  float nv0[3], nv1[3], nv2[3], nv3[3], nv4[3], nv5[3], nv6[3], nv7[3], nv8[3],
      nv9[3], nv10[3], nv11[3];
  for (int i = 0; i < 3; i++) {
    pv0[i] = (2 * verts[v[0]][i] + verts[v[1]][i]) / 3;
    pv1[i] = (2 * verts[v[1]][i] + verts[v[0]][i]) / 3;
    pv2[i] = (2 * verts[v[1]][i] + verts[v[2]][i]) / 3;
    pv3[i] = (2 * verts[v[2]][i] + verts[v[1]][i]) / 3;
    pv4[i] = (2 * verts[v[2]][i] + verts[v[3]][i]) / 3;
    pv5[i] = (2 * verts[v[3]][i] + verts[v[2]][i]) / 3;
    pv6[i] = (2 * verts[v[0]][i] + verts[v[3]][i]) / 3;
    pv7[i] = (2 * verts[v[3]][i] + verts[v[0]][i]) / 3;
    pv8[i] = (2 * pv0[i] + pv5[i]) / 3;
    pv9[i] = (2 * pv5[i] + pv0[i]) / 3;
    pv10[i] = (2 * pv1[i] + pv4[i]) / 3;
    pv11[i] = (2 * pv4[i] + pv1[i]) / 3;
    nv0[i] = (2 * normals[v[0]][i] + normals[v[1]][i]) / 3;
    nv1[i] = (2 * normals[v[1]][i] + normals[v[0]][i]) / 3;
    nv2[i] = (2 * normals[v[1]][i] + normals[v[2]][i]) / 3;
    nv3[i] = (2 * normals[v[2]][i] + normals[v[1]][i]) / 3;
    nv4[i] = (2 * normals[v[2]][i] + normals[v[3]][i]) / 3;
    nv5[i] = (2 * normals[v[3]][i] + normals[v[2]][i]) / 3;
    nv6[i] = (2 * normals[v[0]][i] + normals[v[3]][i]) / 3;
    nv7[i] = (2 * normals[v[3]][i] + normals[v[0]][i]) / 3;
    nv8[i] = (2 * nv0[i] + nv5[i]) / 3;
    nv9[i] = (2 * nv5[i] + nv0[i]) / 3;
    nv10[i] = (2 * nv1[i] + nv4[i]) / 3;
    nv11[i] = (2 * nv4[i] + nv1[i]) / 3;
  }
  v_new[4] = AddVert(pv8, nv8);
  v_new[5] = AddVert(pv10, nv10);
  v_new[6] = AddVert(pv11, nv11);
  v_new[7] = AddVert(pv9, nv9);
  AddBound(v_new[4], 1);
  AddBound(v_new[5], 1);
  AddBound(v_new[6], 1);
  AddBound(v_new[7], 1);
  v_new[0] = AddRefine_edgevtx(v[0], v[1]);
  v_new[1] = AddRefine_edgevtx(v[2], v[1]);
  v_new[2] = AddRefine_edgevtx(v[2], v[3]);
  v_new[3] = AddRefine_edgevtx(v[0], v[3]);
}

void Geoframe::AddQuad_adaptive_4_2b(uint *v, uint *v_new, int num) {
  uint vv[8], v_quad[4];
  for (int i = 0; i < 8; i++) {
    vv[i] = v_new[i];
  }
  v_quad[0] = v[0];
  v_quad[1] = vv[0];
  v_quad[2] = vv[4];
  v_quad[3] = vv[3];
  AddQuad(v_quad, num);
  v_quad[0] = v[1];
  v_quad[1] = vv[5];
  v_quad[2] = vv[4];
  v_quad[3] = vv[0];
  AddQuad(v_quad, num);
  v_quad[0] = v[1];
  v_quad[1] = vv[1];
  v_quad[2] = vv[6];
  v_quad[3] = vv[5];
  AddQuad(v_quad, num);
  v_quad[0] = v[2];
  v_quad[1] = vv[2];
  v_quad[2] = vv[6];
  v_quad[3] = vv[1];
  AddQuad(v_quad, num);
  v_quad[0] = v[3];
  v_quad[1] = vv[7];
  v_quad[2] = vv[6];
  v_quad[3] = vv[2];
  AddQuad(v_quad, num);
  v_quad[0] = v[3];
  v_quad[1] = vv[3];
  v_quad[2] = vv[4];
  v_quad[3] = vv[7];
  AddQuad(v_quad, num);
  v_quad[0] = vv[4];
  v_quad[1] = vv[5];
  v_quad[2] = vv[6];
  v_quad[3] = vv[7];
  AddQuad(v_quad, num);
}

void Geoframe::AddQuad_hexa(uint *v, int num) {
  float pv0[3], pv1[3], pv2[3], pv3[3], pv4[3], pv5[3], pv6[3], norm[3];
  uint vv[7], v_quad[4];
  for (int i = 0; i < 3; i++) {
    pv0[i] = (verts[v[0]][i] + verts[v[1]][i]) / 2;
    pv1[i] = (verts[v[1]][i] + verts[v[2]][i]) / 2;
    pv2[i] = (verts[v[2]][i] + verts[v[3]][i]) / 2;
    pv3[i] = (verts[v[3]][i] + verts[v[0]][i]) / 2;
    pv4[i] = (verts[v[0]][i] + verts[v[2]][i]) / 2;
    pv5[i] = (verts[v[0]][i] + verts[v[1]][i] + verts[v[2]][i]) / 3;
    pv6[i] = (verts[v[0]][i] + verts[v[2]][i] + verts[v[3]][i]) / 3;
  }
  vv[0] = AddVert(pv0, norm);
  vv[1] = AddVert(pv1, norm);
  vv[2] = AddVert(pv2, norm);
  vv[3] = AddVert(pv3, norm);
  vv[4] = AddVert(pv4, norm);
  vv[5] = AddVert(pv5, norm);
  vv[6] = AddVert(pv6, norm);
  AddBound(vv[0], 1);
  AddBound(vv[1], 1);
  AddBound(vv[2], 1);
  AddBound(vv[3], 1);
  AddBound(vv[4], 1);
  AddBound(vv[5], 1);
  AddBound(vv[6], 1);
  v_quad[0] = v[0];
  v_quad[1] = vv[0];
  v_quad[2] = vv[5];
  v_quad[3] = vv[4];
  AddQuad(v_quad, num);
  v_quad[0] = v[1];
  v_quad[1] = vv[1];
  v_quad[2] = vv[5];
  v_quad[3] = vv[0];
  AddQuad(v_quad, num);
  v_quad[0] = v[2];
  v_quad[1] = vv[4];
  v_quad[2] = vv[5];
  v_quad[3] = vv[1];
  AddQuad(v_quad, num);
  v_quad[0] = v[2];
  v_quad[1] = vv[2];
  v_quad[2] = vv[6];
  v_quad[3] = vv[4];
  AddQuad(v_quad, num);
  v_quad[0] = v[3];
  v_quad[1] = vv[3];
  v_quad[2] = vv[6];
  v_quad[3] = vv[2];
  AddQuad(v_quad, num);
  v_quad[0] = v[0];
  v_quad[1] = vv[4];
  v_quad[2] = vv[6];
  v_quad[3] = vv[3];
  AddQuad(v_quad, num);
}

void Geoframe::AddVert_hexa_adaptive_1_center(uint *v, uint *v_new) {
  float pv[8][3], nv[8][3], po[3], no[3];
  for (int i = 0; i < 3; i++) {
    po[i] = 0;
    no[i] = 0;
    for (int j = 0; j < 8; j++) {
      po[i] += verts[v[j]][i];
      no[i] += normals[v[j]][i];
    }
    po[i] /= 8;
    no[i] /= 8;
    for (int j = 0; j < 8; j++) {
      pv[j][i] = (2 * po[i] + verts[v[j]][i]) / 3;
      nv[j][i] = (2 * no[i] + normals[v[j]][i]) / 3;
    }
  }
  for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 3; j++) {
      po[j] = pv[i][j];
      no[j] = nv[i][j];
    }
    v_new[i] = AddVert(po, no);
  }
}

void Geoframe::AddVert_hexa_adaptive_1_top(uint *v, uint *v_new) {
  float pv[8][3], nv[8][3], po[3], no[3];
  for (int i = 0; i < 3; i++) {
    po[i] = 0;
    no[i] = 0;
    for (int j = 0; j < 4; j++) {
      po[i] += verts[v[j]][i];
      no[i] += normals[v[j]][i];
    }
    po[i] /= 4;
    no[i] /= 4;
    for (int j = 0; j < 4; j++) {
      pv[j][i] = (2 * po[i] + verts[v[j]][i]) / 3;
      nv[j][i] = (2 * no[i] + normals[v[j]][i]) / 3;
    }
    for (int j = 4; j < 8; j++) {
      pv[j][i] = (po[i] + verts[v[j]][i]) / 2;
      nv[j][i] = (no[i] + normals[v[j]][i]) / 2;
    }
  }
  for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 3; j++) {
      po[j] = pv[i][j];
      no[j] = nv[i][j];
    }
    v_new[i] = AddVert(po, no);
    if (i < 4) {
      AddBound(v_new[i], 1);
    }
  }
}

void Geoframe::AddVert_hexa_adaptive_2(uint *v, int *edge_id, uint *v_new) {
  float pv[64][3], nv[64][3], po[3], no[3];
  for (int i = 0; i < 3; i++) {
    pv[0][i] = verts[v[0]][i];
    pv[3][i] = verts[v[1]][i];
    pv[51][i] = verts[v[5]][i];
    pv[48][i] = verts[v[4]][i];
    pv[12][i] = verts[v[3]][i];
    pv[15][i] = verts[v[2]][i];
    pv[63][i] = verts[v[6]][i];
    pv[60][i] = verts[v[7]][i];
    nv[0][i] = normals[v[0]][i];
    nv[3][i] = normals[v[1]][i];
    nv[51][i] = normals[v[5]][i];
    nv[48][i] = normals[v[4]][i];
    nv[12][i] = normals[v[3]][i];
    nv[15][i] = normals[v[2]][i];
    nv[63][i] = normals[v[6]][i];
    nv[60][i] = normals[v[7]][i];
    pv[16][i] = (2 * pv[0][i] + pv[48][i]) / 3;
    pv[32][i] = (2 * pv[48][i] + pv[0][i]) / 3;
    pv[19][i] = (2 * pv[3][i] + pv[51][i]) / 3;
    pv[35][i] = (2 * pv[51][i] + pv[3][i]) / 3;
    pv[28][i] = (2 * pv[12][i] + pv[60][i]) / 3;
    pv[44][i] = (2 * pv[60][i] + pv[12][i]) / 3;
    pv[31][i] = (2 * pv[15][i] + pv[63][i]) / 3;
    pv[47][i] = (2 * pv[63][i] + pv[15][i]) / 3;
    nv[16][i] = (2 * nv[0][i] + nv[48][i]) / 3;
    nv[32][i] = (2 * nv[48][i] + nv[0][i]) / 3;
    nv[19][i] = (2 * nv[3][i] + nv[51][i]) / 3;
    nv[35][i] = (2 * nv[51][i] + nv[3][i]) / 3;
    nv[28][i] = (2 * nv[12][i] + nv[60][i]) / 3;
    nv[44][i] = (2 * nv[60][i] + nv[12][i]) / 3;
    nv[31][i] = (2 * nv[15][i] + nv[63][i]) / 3;
    nv[47][i] = (2 * nv[63][i] + nv[15][i]) / 3;
    for (int j = 0; j < 4; j++) {
      pv[4 + 16 * j][i] = (2 * pv[0 + 16 * j][i] + pv[12 + 16 * j][i]) / 3;
      pv[8 + 16 * j][i] = (2 * pv[12 + 16 * j][i] + pv[0 + 16 * j][i]) / 3;
      pv[7 + 16 * j][i] = (2 * pv[3 + 16 * j][i] + pv[15 + 16 * j][i]) / 3;
      pv[11 + 16 * j][i] = (2 * pv[15 + 16 * j][i] + pv[3 + 16 * j][i]) / 3;
      pv[1 + 16 * j][i] = (2 * pv[0 + 16 * j][i] + pv[3 + 16 * j][i]) / 3;
      pv[2 + 16 * j][i] = (2 * pv[3 + 16 * j][i] + pv[0 + 16 * j][i]) / 3;
      pv[13 + 16 * j][i] = (2 * pv[12 + 16 * j][i] + pv[15 + 16 * j][i]) / 3;
      pv[14 + 16 * j][i] = (2 * pv[15 + 16 * j][i] + pv[12 + 16 * j][i]) / 3;
      pv[5 + 16 * j][i] = (2 * pv[4 + 16 * j][i] + pv[7 + 16 * j][i]) / 3;
      pv[6 + 16 * j][i] = (2 * pv[7 + 16 * j][i] + pv[4 + 16 * j][i]) / 3;
      pv[9 + 16 * j][i] = (2 * pv[8 + 16 * j][i] + pv[11 + 16 * j][i]) / 3;
      pv[10 + 16 * j][i] = (2 * pv[11 + 16 * j][i] + pv[8 + 16 * j][i]) / 3;
      nv[4 + 16 * j][i] = (2 * nv[0 + 16 * j][i] + nv[12 + 16 * j][i]) / 3;
      nv[8 + 16 * j][i] = (2 * nv[12 + 16 * j][i] + nv[0 + 16 * j][i]) / 3;
      nv[7 + 16 * j][i] = (2 * nv[3 + 16 * j][i] + nv[15 + 16 * j][i]) / 3;
      nv[11 + 16 * j][i] = (2 * nv[15 + 16 * j][i] + nv[3 + 16 * j][i]) / 3;
      nv[1 + 16 * j][i] = (2 * nv[0 + 16 * j][i] + nv[3 + 16 * j][i]) / 3;
      nv[2 + 16 * j][i] = (2 * nv[3 + 16 * j][i] + nv[0 + 16 * j][i]) / 3;
      nv[13 + 16 * j][i] = (2 * nv[12 + 16 * j][i] + nv[15 + 16 * j][i]) / 3;
      nv[14 + 16 * j][i] = (2 * nv[15 + 16 * j][i] + nv[12 + 16 * j][i]) / 3;
      nv[5 + 16 * j][i] = (2 * nv[4 + 16 * j][i] + nv[7 + 16 * j][i]) / 3;
      nv[6 + 16 * j][i] = (2 * nv[7 + 16 * j][i] + nv[4 + 16 * j][i]) / 3;
      nv[9 + 16 * j][i] = (2 * nv[8 + 16 * j][i] + nv[11 + 16 * j][i]) / 3;
      nv[10 + 16 * j][i] = (2 * nv[11 + 16 * j][i] + nv[8 + 16 * j][i]) / 3;
    }
  }
  v_new[1] = AddRefine_edgevtx_hex(v[0], v[1]);
  v_new[2] = AddRefine_edgevtx_hex(v[1], v[0]);
  v_new[4] = AddRefine_edgevtx_hex(v[0], v[3]);
  v_new[8] = AddRefine_edgevtx_hex(v[3], v[0]);
  v_new[7] = AddRefine_edgevtx_hex(v[1], v[2]);
  v_new[11] = AddRefine_edgevtx_hex(v[2], v[1]);
  v_new[13] = AddRefine_edgevtx_hex(v[3], v[2]);
  v_new[14] = AddRefine_edgevtx_hex(v[2], v[3]);
  v_new[5] = AddRefine_facevtx(v[0], v[2], v[1], v[3]);
  v_new[10] = AddRefine_facevtx(v[2], v[0], v[1], v[3]);
  v_new[6] = AddRefine_facevtx(v[1], v[3], v[0], v[2]);
  v_new[9] = AddRefine_facevtx(v[3], v[1], v[0], v[2]);
  v_new[1 + 48] = AddRefine_edgevtx_hex(v[0 + 4], v[1 + 4]);
  v_new[2 + 48] = AddRefine_edgevtx_hex(v[1 + 4], v[0 + 4]);
  v_new[4 + 48] = AddRefine_edgevtx_hex(v[0 + 4], v[3 + 4]);
  v_new[8 + 48] = AddRefine_edgevtx_hex(v[3 + 4], v[0 + 4]);
  v_new[7 + 48] = AddRefine_edgevtx_hex(v[1 + 4], v[2 + 4]);
  v_new[11 + 48] = AddRefine_edgevtx_hex(v[2 + 4], v[1 + 4]);
  v_new[13 + 48] = AddRefine_edgevtx_hex(v[3 + 4], v[2 + 4]);
  v_new[14 + 48] = AddRefine_edgevtx_hex(v[2 + 4], v[3 + 4]);
  v_new[5 + 48] = AddRefine_facevtx(v[0 + 4], v[2 + 4], v[1 + 4], v[3 + 4]);
  v_new[10 + 48] = AddRefine_facevtx(v[2 + 4], v[0 + 4], v[1 + 4], v[3 + 4]);
  v_new[6 + 48] = AddRefine_facevtx(v[1 + 4], v[3 + 4], v[0 + 4], v[2 + 4]);
  v_new[9 + 48] = AddRefine_facevtx(v[3 + 4], v[1 + 4], v[0 + 4], v[2 + 4]);
  v_new[16] = AddRefine_edgevtx_hex(v[0], v[4]);
  v_new[32] = AddRefine_edgevtx_hex(v[4], v[0]);
  v_new[19] = AddRefine_edgevtx_hex(v[1], v[5]);
  v_new[35] = AddRefine_edgevtx_hex(v[5], v[1]);
  v_new[28] = AddRefine_edgevtx_hex(v[3], v[7]);
  v_new[44] = AddRefine_edgevtx_hex(v[7], v[3]);
  v_new[31] = AddRefine_edgevtx_hex(v[2], v[6]);
  v_new[47] = AddRefine_edgevtx_hex(v[6], v[2]);
  v_new[17] = AddRefine_facevtx(v[0], v[5], v[1], v[4]);
  v_new[34] = AddRefine_facevtx(v[5], v[0], v[1], v[4]);
  v_new[18] = AddRefine_facevtx(v[1], v[4], v[0], v[5]);
  v_new[33] = AddRefine_facevtx(v[4], v[1], v[0], v[5]);
  v_new[23] = AddRefine_facevtx(v[1], v[6], v[2], v[5]);
  v_new[43] = AddRefine_facevtx(v[6], v[1], v[2], v[5]);
  v_new[27] = AddRefine_facevtx(v[2], v[5], v[1], v[6]);
  v_new[39] = AddRefine_facevtx(v[5], v[2], v[1], v[6]);
  v_new[20] = AddRefine_facevtx(v[0], v[7], v[3], v[4]);
  v_new[40] = AddRefine_facevtx(v[7], v[0], v[3], v[4]);
  v_new[24] = AddRefine_facevtx(v[3], v[4], v[0], v[7]);
  v_new[36] = AddRefine_facevtx(v[4], v[3], v[0], v[7]);
  v_new[29] = AddRefine_facevtx(v[3], v[6], v[2], v[7]);
  v_new[46] = AddRefine_facevtx(v[6], v[3], v[2], v[7]);
  v_new[30] = AddRefine_facevtx(v[2], v[7], v[3], v[6]);
  v_new[45] = AddRefine_facevtx(v[7], v[2], v[3], v[6]);
  for (int i = 0; i < 64; i++) {
    if (i == 0) {
      v_new[i] = v[0];
    } else if (i == 3) {
      v_new[i] = v[1];
    } else if (i == 51) {
      v_new[i] = v[5];
    } else if (i == 48) {
      v_new[i] = v[4];
    } else if (i == 12) {
      v_new[i] = v[3];
    } else if (i == 15) {
      v_new[i] = v[2];
    } else if (i == 63) {
      v_new[i] = v[6];
    } else if (i == 60) {
      v_new[i] = v[7];
    } else {
      if (i == 21 || i == 22 || i == 25 || i == 26 || i == 37 || i == 38 ||
          i == 41 || i == 42) {
        for (int j = 0; j < 3; j++) {
          po[j] = pv[i][j];
          no[j] = nv[i][j];
        }
        v_new[i] = AddVert(po, no);
      }
    }
    if ((edge_id[0] > 0) && (i < 16)) {
      AddBound(v_new[i], 1);
    }
    if ((edge_id[1] > 0) && (i > 47)) {
      AddBound(v_new[i], 1);
    }
    if ((edge_id[2] > 0) && ((i % 4) == 0)) {
      AddBound(v_new[i], 1);
    }
    if ((edge_id[3] > 0) && ((i % 4) == 3)) {
      AddBound(v_new[i], 1);
    }
    if ((edge_id[4] > 0) && ((i % 16) < 4)) {
      AddBound(v_new[i], 1);
    }
    if ((edge_id[5] > 0) && ((i % 16) > 11)) {
      AddBound(v_new[i], 1);
    }
  }
  // define bound_edge ...
  if (CheckBound_edge(v_new[0], v_new[3]) == 1) {
    AddBound(v_new[1], 1);
    AddBound(v_new[2], 1);
  }
  if (CheckBound_edge(v_new[3], v_new[15]) == 1) {
    AddBound(v_new[7], 1);
    AddBound(v_new[11], 1);
  }
  if (CheckBound_edge(v_new[15], v_new[12]) == 1) {
    AddBound(v_new[13], 1);
    AddBound(v_new[14], 1);
  }
  if (CheckBound_edge(v_new[0], v_new[12]) == 1) {
    AddBound(v_new[4], 1);
    AddBound(v_new[8], 1);
  }
  if (CheckBound_edge(v_new[48], v_new[51]) == 1) {
    AddBound(v_new[49], 1);
    AddBound(v_new[50], 1);
  }
  if (CheckBound_edge(v_new[51], v_new[63]) == 1) {
    AddBound(v_new[55], 1);
    AddBound(v_new[59], 1);
  }
  if (CheckBound_edge(v_new[63], v_new[60]) == 1) {
    AddBound(v_new[61], 1);
    AddBound(v_new[62], 1);
  }
  if (CheckBound_edge(v_new[48], v_new[60]) == 1) {
    AddBound(v_new[52], 1);
    AddBound(v_new[56], 1);
  }
  if (CheckBound_edge(v_new[0], v_new[48]) == 1) {
    AddBound(v_new[16], 1);
    AddBound(v_new[32], 1);
  }
  if (CheckBound_edge(v_new[3], v_new[51]) == 1) {
    AddBound(v_new[19], 1);
    AddBound(v_new[35], 1);
  }
  if (CheckBound_edge(v_new[12], v_new[60]) == 1) {
    AddBound(v_new[28], 1);
    AddBound(v_new[44], 1);
  }
  if (CheckBound_edge(v_new[15], v_new[63]) == 1) {
    AddBound(v_new[31], 1);
    AddBound(v_new[47], 1);
  }
}

void Geoframe::AddVert_hexa_adaptive_2_1(uint *v, int *edge_id, uint *v_new) {
  float pv[7][3], nv[7][3], po[3], t0[3], t1[3];
  for (int i = 0; i < 3; i++) {
    pv[0][i] = (2 * verts[v[0]][i] + verts[v[1]][i]) / 3;
    po[i] = (2 * verts[v[3]][i] + verts[v[2]][i]) / 3;
    pv[1][i] = (2 * pv[0][i] + po[i]) / 3;
    pv[2][i] = (2 * verts[v[0]][i] + verts[v[3]][i]) / 3;
    pv[3][i] = (2 * verts[v[0]][i] + verts[v[4]][i]) / 3;
    po[i] = (2 * verts[v[1]][i] + verts[v[5]][i]) / 3;
    pv[4][i] = (2 * pv[3][i] + po[i]) / 3;
    t0[i] = (2 * verts[v[4]][i] + verts[v[5]][i]) / 3;
    t1[i] = (2 * verts[v[7]][i] + verts[v[6]][i]) / 3;
    po[i] = (2 * t0[i] + t1[i]) / 3;
    pv[5][i] = (2 * pv[1][i] + po[i]) / 3;
    po[i] = (2 * verts[v[3]][i] + verts[v[7]][i]) / 3;
    pv[6][i] = (2 * pv[3][i] + po[i]) / 3;
    // normal
    nv[0][i] = (2 * normals[v[0]][i] + normals[v[1]][i]) / 3;
    po[i] = (2 * normals[v[3]][i] + normals[v[2]][i]) / 3;
    nv[1][i] = (2 * nv[0][i] + po[i]) / 3;
    nv[2][i] = (2 * normals[v[0]][i] + normals[v[3]][i]) / 3;
    nv[3][i] = (2 * normals[v[0]][i] + normals[v[4]][i]) / 3;
    po[i] = (2 * normals[v[1]][i] + normals[v[5]][i]) / 3;
    nv[4][i] = (2 * nv[3][i] + po[i]) / 3;
    t0[i] = (2 * normals[v[4]][i] + normals[v[5]][i]) / 3;
    t1[i] = (2 * normals[v[7]][i] + normals[v[6]][i]) / 3;
    po[i] = (2 * t0[i] + t1[i]) / 3;
    nv[5][i] = (2 * nv[1][i] + po[i]) / 3;
    po[i] = (2 * normals[v[3]][i] + normals[v[7]][i]) / 3;
    nv[6][i] = (2 * nv[3][i] + po[i]) / 3;
  }
  v_new[0] = AddRefine_edgevtx_hex(v[0], v[1]);
  v_new[1] = AddRefine_facevtx(v[0], v[2], v[1], v[3]);
  v_new[2] = AddRefine_edgevtx_hex(v[0], v[3]);
  v_new[3] = AddRefine_edgevtx_hex(v[0], v[4]);
  v_new[4] = AddRefine_facevtx(v[0], v[5], v[1], v[4]);
  v_new[6] = AddRefine_facevtx(v[0], v[7], v[3], v[4]);
  for (int j = 0; j < 3; j++) {
    po[j] = pv[5][j];
    t0[j] = nv[5][j];
  }
  v_new[5] = AddVert(po, t0);
  if (edge_id[0] > 0) {
    AddBound(v_new[0], 1);
    AddBound(v_new[1], 1);
    AddBound(v_new[2], 1);
  }
  if (edge_id[2] > 0) {
    AddBound(v_new[2], 1);
    AddBound(v_new[3], 1);
    AddBound(v_new[6], 1);
  }
  if (edge_id[4] > 0) {
    AddBound(v_new[0], 1);
    AddBound(v_new[3], 1);
    AddBound(v_new[4], 1);
  }
  if (CheckBound_edge(v[0], v[1]) == 1) {
    AddBound(v_new[0], 1);
  }
  if (CheckBound_edge(v[0], v[3]) == 1) {
    AddBound(v_new[2], 1);
  }
  if (CheckBound_edge(v[0], v[4]) == 1) {
    AddBound(v_new[3], 1);
  }
}

void Geoframe::AddVert_hexa_adaptive_2_2(uint *v, int *edge_id, uint *v_new) {
  float pv[64][3], nv[64][3], po[3], no[3];
  for (int i = 0; i < 3; i++) {
    pv[0][i] = verts[v[0]][i];
    pv[3][i] = verts[v[1]][i];
    pv[51][i] = verts[v[5]][i];
    pv[48][i] = verts[v[4]][i];
    pv[12][i] = verts[v[3]][i];
    pv[15][i] = verts[v[2]][i];
    pv[63][i] = verts[v[6]][i];
    pv[60][i] = verts[v[7]][i];
    nv[0][i] = normals[v[0]][i];
    nv[3][i] = normals[v[1]][i];
    nv[51][i] = normals[v[5]][i];
    nv[48][i] = normals[v[4]][i];
    nv[12][i] = normals[v[3]][i];
    nv[15][i] = normals[v[2]][i];
    nv[63][i] = normals[v[6]][i];
    nv[60][i] = normals[v[7]][i];
    pv[16][i] = (2 * pv[0][i] + pv[48][i]) / 3;
    pv[32][i] = (2 * pv[48][i] + pv[0][i]) / 3;
    pv[19][i] = (2 * pv[3][i] + pv[51][i]) / 3;
    pv[35][i] = (2 * pv[51][i] + pv[3][i]) / 3;
    pv[28][i] = (2 * pv[12][i] + pv[60][i]) / 3;
    pv[44][i] = (2 * pv[60][i] + pv[12][i]) / 3;
    pv[31][i] = (2 * pv[15][i] + pv[63][i]) / 3;
    pv[47][i] = (2 * pv[63][i] + pv[15][i]) / 3;
    nv[16][i] = (2 * nv[0][i] + nv[48][i]) / 3;
    nv[32][i] = (2 * nv[48][i] + nv[0][i]) / 3;
    nv[19][i] = (2 * nv[3][i] + nv[51][i]) / 3;
    nv[35][i] = (2 * nv[51][i] + nv[3][i]) / 3;
    nv[28][i] = (2 * nv[12][i] + nv[60][i]) / 3;
    nv[44][i] = (2 * nv[60][i] + nv[12][i]) / 3;
    nv[31][i] = (2 * nv[15][i] + nv[63][i]) / 3;
    nv[47][i] = (2 * nv[63][i] + nv[15][i]) / 3;
    for (int j = 0; j < 4; j++) {
      pv[4 + 16 * j][i] = (2 * pv[0 + 16 * j][i] + pv[12 + 16 * j][i]) / 3;
      pv[8 + 16 * j][i] = (2 * pv[12 + 16 * j][i] + pv[0 + 16 * j][i]) / 3;
      pv[7 + 16 * j][i] = (2 * pv[3 + 16 * j][i] + pv[15 + 16 * j][i]) / 3;
      pv[11 + 16 * j][i] = (2 * pv[15 + 16 * j][i] + pv[3 + 16 * j][i]) / 3;
      pv[1 + 16 * j][i] = (2 * pv[0 + 16 * j][i] + pv[3 + 16 * j][i]) / 3;
      pv[2 + 16 * j][i] = (2 * pv[3 + 16 * j][i] + pv[0 + 16 * j][i]) / 3;
      pv[13 + 16 * j][i] = (2 * pv[12 + 16 * j][i] + pv[15 + 16 * j][i]) / 3;
      pv[14 + 16 * j][i] = (2 * pv[15 + 16 * j][i] + pv[12 + 16 * j][i]) / 3;
      pv[5 + 16 * j][i] = (2 * pv[4 + 16 * j][i] + pv[7 + 16 * j][i]) / 3;
      pv[6 + 16 * j][i] = (2 * pv[7 + 16 * j][i] + pv[4 + 16 * j][i]) / 3;
      pv[9 + 16 * j][i] = (2 * pv[8 + 16 * j][i] + pv[11 + 16 * j][i]) / 3;
      pv[10 + 16 * j][i] = (2 * pv[11 + 16 * j][i] + pv[8 + 16 * j][i]) / 3;
      nv[4 + 16 * j][i] = (2 * nv[0 + 16 * j][i] + nv[12 + 16 * j][i]) / 3;
      nv[8 + 16 * j][i] = (2 * nv[12 + 16 * j][i] + nv[0 + 16 * j][i]) / 3;
      nv[7 + 16 * j][i] = (2 * nv[3 + 16 * j][i] + nv[15 + 16 * j][i]) / 3;
      nv[11 + 16 * j][i] = (2 * nv[15 + 16 * j][i] + nv[3 + 16 * j][i]) / 3;
      nv[1 + 16 * j][i] = (2 * nv[0 + 16 * j][i] + nv[3 + 16 * j][i]) / 3;
      nv[2 + 16 * j][i] = (2 * nv[3 + 16 * j][i] + nv[0 + 16 * j][i]) / 3;
      nv[13 + 16 * j][i] = (2 * nv[12 + 16 * j][i] + nv[15 + 16 * j][i]) / 3;
      nv[14 + 16 * j][i] = (2 * nv[15 + 16 * j][i] + nv[12 + 16 * j][i]) / 3;
      nv[5 + 16 * j][i] = (2 * nv[4 + 16 * j][i] + nv[7 + 16 * j][i]) / 3;
      nv[6 + 16 * j][i] = (2 * nv[7 + 16 * j][i] + nv[4 + 16 * j][i]) / 3;
      nv[9 + 16 * j][i] = (2 * nv[8 + 16 * j][i] + nv[11 + 16 * j][i]) / 3;
      nv[10 + 16 * j][i] = (2 * nv[11 + 16 * j][i] + nv[8 + 16 * j][i]) / 3;
    }
  }
  v_new[1] = AddRefine_edgevtx_hex(v[0], v[1]);
  v_new[2] = AddRefine_edgevtx_hex(v[1], v[0]);
  v_new[4] = AddRefine_edgevtx_hex(v[0], v[3]);
  v_new[7] = AddRefine_edgevtx_hex(v[1], v[2]);
  v_new[5] = AddRefine_facevtx(v[0], v[2], v[1], v[3]);
  v_new[9] = AddRefine_facevtx(v[2], v[0], v[1], v[3]);
  v_new[6] = AddRefine_facevtx(v[1], v[3], v[0], v[2]);
  v_new[8] = AddRefine_facevtx(v[3], v[1], v[0], v[2]);
  v_new[12] = AddRefine_edgevtx_hex(v[0], v[4]);
  v_new[15] = AddRefine_edgevtx_hex(v[1], v[5]);
  v_new[13] = AddRefine_facevtx(v[0], v[5], v[1], v[4]);
  v_new[21] = AddRefine_facevtx(v[5], v[0], v[1], v[4]);
  v_new[14] = AddRefine_facevtx(v[1], v[4], v[0], v[5]);
  v_new[20] = AddRefine_facevtx(v[4], v[1], v[0], v[5]);
  v_new[16] = AddRefine_facevtx(v[0], v[7], v[3], v[4]);
  v_new[19] = AddRefine_facevtx(v[1], v[6], v[2], v[5]);
  for (int i = 0; i < 28; i++) {
    if (i == 0) {
      v_new[i] = v[0];
    } else if (i == 3) {
      v_new[i] = v[1];
    } else if (i == 11) {
      v_new[i] = v[2];
    } else if (i == 10) {
      v_new[i] = v[3];
    } else if (i == 24) {
      v_new[i] = v[4];
    } else if (i == 25) {
      v_new[i] = v[5];
    } else if (i == 27) {
      v_new[i] = v[6];
    } else if (i == 26) {
      v_new[i] = v[7];
    } else {
      if (i == 17 || i == 18 || i == 22 || i == 23) {
        int vv;
        if (i < 20) {
          vv = i + 4;
        } else {
          vv = i + 19;
        }
        for (int j = 0; j < 3; j++) {
          po[j] = pv[vv][j];
          no[j] = nv[vv][j];
        }
        v_new[i] = AddVert(po, no);
      }
    }
    if ((edge_id[0] > 0) && (i < 12)) {
      AddBound(v_new[i], 1);
    }
    if ((edge_id[2] > 0) && (i == 4 || i == 12 || i == 16)) {
      AddBound(v_new[i], 1);
    }
    if ((edge_id[3] > 0) && (i == 7 || i == 15 || i == 19)) {
      AddBound(v_new[i], 1);
    }
    if ((edge_id[4] > 0) && ((i % 12 < 4 && i < 16) || i == 20 || i == 21)) {
      AddBound(v_new[i], 1);
    }
  }
  if (CheckBound_edge(v_new[0], v_new[3]) == 1) {
    AddBound(v_new[1], 1);
    AddBound(v_new[2], 1);
  }
  if (CheckBound_edge(v_new[3], v_new[11]) == 1) {
    AddBound(v_new[7], 1);
  }
  if (CheckBound_edge(v_new[0], v_new[10]) == 1) {
    AddBound(v_new[4], 1);
  }
  if (CheckBound_edge(v_new[0], v_new[24]) == 1) {
    AddBound(v_new[12], 1);
  }
  if (CheckBound_edge(v_new[3], v_new[25]) == 1) {
    AddBound(v_new[15], 1);
  }
}

void Geoframe::AddVert_hexa_adaptive_2_4(uint *v, int *edge_id, uint *v_new) {
  float pv[68][3], nv[68][3], po[3], no[3];
  for (int i = 0; i < 3; i++) {
    pv[0][i] = verts[v[0]][i];
    pv[3][i] = verts[v[1]][i];
    pv[51][i] = verts[v[5]][i];
    pv[48][i] = verts[v[4]][i];
    pv[12][i] = verts[v[3]][i];
    pv[15][i] = verts[v[2]][i];
    pv[63][i] = verts[v[6]][i];
    pv[60][i] = verts[v[7]][i];
    nv[0][i] = normals[v[0]][i];
    nv[3][i] = normals[v[1]][i];
    nv[51][i] = normals[v[5]][i];
    nv[48][i] = normals[v[4]][i];
    nv[12][i] = normals[v[3]][i];
    nv[15][i] = normals[v[2]][i];
    nv[63][i] = normals[v[6]][i];
    nv[60][i] = normals[v[7]][i];
    pv[16][i] = (2 * pv[0][i] + pv[48][i]) / 3;
    pv[32][i] = (2 * pv[48][i] + pv[0][i]) / 3;
    pv[19][i] = (2 * pv[3][i] + pv[51][i]) / 3;
    pv[35][i] = (2 * pv[51][i] + pv[3][i]) / 3;
    pv[28][i] = (2 * pv[12][i] + pv[60][i]) / 3;
    pv[44][i] = (2 * pv[60][i] + pv[12][i]) / 3;
    pv[31][i] = (2 * pv[15][i] + pv[63][i]) / 3;
    pv[47][i] = (2 * pv[63][i] + pv[15][i]) / 3;
    nv[16][i] = (2 * nv[0][i] + nv[48][i]) / 3;
    nv[32][i] = (2 * nv[48][i] + nv[0][i]) / 3;
    nv[19][i] = (2 * nv[3][i] + nv[51][i]) / 3;
    nv[35][i] = (2 * nv[51][i] + nv[3][i]) / 3;
    nv[28][i] = (2 * nv[12][i] + nv[60][i]) / 3;
    nv[44][i] = (2 * nv[60][i] + nv[12][i]) / 3;
    nv[31][i] = (2 * nv[15][i] + nv[63][i]) / 3;
    nv[47][i] = (2 * nv[63][i] + nv[15][i]) / 3;
    for (int j = 0; j < 4; j++) {
      pv[4 + 16 * j][i] = (2 * pv[0 + 16 * j][i] + pv[12 + 16 * j][i]) / 3;
      pv[8 + 16 * j][i] = (2 * pv[12 + 16 * j][i] + pv[0 + 16 * j][i]) / 3;
      pv[7 + 16 * j][i] = (2 * pv[3 + 16 * j][i] + pv[15 + 16 * j][i]) / 3;
      pv[11 + 16 * j][i] = (2 * pv[15 + 16 * j][i] + pv[3 + 16 * j][i]) / 3;
      pv[1 + 16 * j][i] = (2 * pv[0 + 16 * j][i] + pv[3 + 16 * j][i]) / 3;
      pv[2 + 16 * j][i] = (2 * pv[3 + 16 * j][i] + pv[0 + 16 * j][i]) / 3;
      pv[13 + 16 * j][i] = (2 * pv[12 + 16 * j][i] + pv[15 + 16 * j][i]) / 3;
      pv[14 + 16 * j][i] = (2 * pv[15 + 16 * j][i] + pv[12 + 16 * j][i]) / 3;
      pv[5 + 16 * j][i] = (2 * pv[4 + 16 * j][i] + pv[7 + 16 * j][i]) / 3;
      pv[6 + 16 * j][i] = (2 * pv[7 + 16 * j][i] + pv[4 + 16 * j][i]) / 3;
      pv[9 + 16 * j][i] = (2 * pv[8 + 16 * j][i] + pv[11 + 16 * j][i]) / 3;
      pv[10 + 16 * j][i] = (2 * pv[11 + 16 * j][i] + pv[8 + 16 * j][i]) / 3;
      nv[4 + 16 * j][i] = (2 * nv[0 + 16 * j][i] + nv[12 + 16 * j][i]) / 3;
      nv[8 + 16 * j][i] = (2 * nv[12 + 16 * j][i] + nv[0 + 16 * j][i]) / 3;
      nv[7 + 16 * j][i] = (2 * nv[3 + 16 * j][i] + nv[15 + 16 * j][i]) / 3;
      nv[11 + 16 * j][i] = (2 * nv[15 + 16 * j][i] + nv[3 + 16 * j][i]) / 3;
      nv[1 + 16 * j][i] = (2 * nv[0 + 16 * j][i] + nv[3 + 16 * j][i]) / 3;
      nv[2 + 16 * j][i] = (2 * nv[3 + 16 * j][i] + nv[0 + 16 * j][i]) / 3;
      nv[13 + 16 * j][i] = (2 * nv[12 + 16 * j][i] + nv[15 + 16 * j][i]) / 3;
      nv[14 + 16 * j][i] = (2 * nv[15 + 16 * j][i] + nv[12 + 16 * j][i]) / 3;
      nv[5 + 16 * j][i] = (2 * nv[4 + 16 * j][i] + nv[7 + 16 * j][i]) / 3;
      nv[6 + 16 * j][i] = (2 * nv[7 + 16 * j][i] + nv[4 + 16 * j][i]) / 3;
      nv[9 + 16 * j][i] = (2 * nv[8 + 16 * j][i] + nv[11 + 16 * j][i]) / 3;
      nv[10 + 16 * j][i] = (2 * nv[11 + 16 * j][i] + nv[8 + 16 * j][i]) / 3;
    }
    pv[64][i] = (pv[21][i] + pv[37][i]) / 2;
    pv[65][i] = (pv[22][i] + pv[38][i]) / 2;
    pv[66][i] = (pv[25][i] + pv[41][i]) / 2;
    pv[67][i] = (pv[26][i] + pv[41][i]) / 2;
    nv[64][i] = (nv[21][i] + nv[37][i]) / 2;
    nv[65][i] = (nv[22][i] + nv[38][i]) / 2;
    nv[66][i] = (nv[25][i] + nv[41][i]) / 2;
    nv[67][i] = (nv[26][i] + nv[41][i]) / 2;
  }
  v_new[1] = AddRefine_edgevtx_hex(v[0], v[1]);
  v_new[2] = AddRefine_edgevtx_hex(v[1], v[0]);
  v_new[4] = AddRefine_edgevtx_hex(v[0], v[3]);
  v_new[8] = AddRefine_edgevtx_hex(v[3], v[0]);
  v_new[7] = AddRefine_edgevtx_hex(v[1], v[2]);
  v_new[11] = AddRefine_edgevtx_hex(v[2], v[1]);
  v_new[13] = AddRefine_edgevtx_hex(v[3], v[2]);
  v_new[14] = AddRefine_edgevtx_hex(v[2], v[3]);
  v_new[5] = AddRefine_facevtx(v[0], v[2], v[1], v[3]);
  v_new[10] = AddRefine_facevtx(v[2], v[0], v[1], v[3]);
  v_new[6] = AddRefine_facevtx(v[1], v[3], v[0], v[2]);
  v_new[9] = AddRefine_facevtx(v[3], v[1], v[0], v[2]);
  v_new[16] = AddRefine_edgevtx_hex(v[0], v[4]);
  v_new[19] = AddRefine_edgevtx_hex(v[1], v[5]);
  v_new[28] = AddRefine_edgevtx_hex(v[3], v[7]);
  v_new[31] = AddRefine_edgevtx_hex(v[2], v[6]);
  v_new[17] = AddRefine_facevtx(v[0], v[5], v[1], v[4]);
  v_new[37] = AddRefine_facevtx(v[5], v[0], v[1], v[4]);
  v_new[18] = AddRefine_facevtx(v[1], v[4], v[0], v[5]);
  v_new[36] = AddRefine_facevtx(v[4], v[1], v[0], v[5]);
  v_new[23] = AddRefine_facevtx(v[1], v[6], v[2], v[5]);
  v_new[41] = AddRefine_facevtx(v[6], v[1], v[2], v[5]);
  v_new[27] = AddRefine_facevtx(v[2], v[5], v[1], v[6]);
  v_new[39] = AddRefine_facevtx(v[5], v[2], v[1], v[6]);
  v_new[20] = AddRefine_facevtx(v[0], v[7], v[3], v[4]);
  v_new[40] = AddRefine_facevtx(v[7], v[0], v[3], v[4]);
  v_new[24] = AddRefine_facevtx(v[3], v[4], v[0], v[7]);
  v_new[38] = AddRefine_facevtx(v[4], v[3], v[0], v[7]);
  v_new[29] = AddRefine_facevtx(v[3], v[6], v[2], v[7]);
  v_new[43] = AddRefine_facevtx(v[6], v[3], v[2], v[7]);
  v_new[30] = AddRefine_facevtx(v[2], v[7], v[3], v[6]);
  v_new[42] = AddRefine_facevtx(v[7], v[2], v[3], v[6]);
  for (int i = 0; i < 48; i++) {
    if (i == 0) {
      v_new[i] = v[0];
    } else if (i == 3) {
      v_new[i] = v[1];
    } else if (i == 15) {
      v_new[i] = v[2];
    } else if (i == 12) {
      v_new[i] = v[3];
    } else if (i == 44) {
      v_new[i] = v[4];
    } else if (i == 45) {
      v_new[i] = v[5];
    } else if (i == 47) {
      v_new[i] = v[6];
    } else if (i == 46) {
      v_new[i] = v[7];
    } else {
      if (i == 21 || i == 22 || i == 25 || i == 26 || (i > 31 && i < 36)) {
        int vv;
        if (i > 31) {
          vv = i + 32;
        } else {
          vv = i;
        }
        for (int j = 0; j < 3; j++) {
          po[j] = pv[vv][j];
          no[j] = nv[vv][j];
        }
        v_new[i] = AddVert(po, no);
      }
    }
    if ((edge_id[0] > 0) && (i < 16)) {
      AddBound(v_new[i], 1);
    }
    if ((edge_id[2] > 0) && ((i % 4 == 0 && i < 29) || i == 38 || i == 40)) {
      AddBound(v_new[i], 1);
    }
    if ((edge_id[3] > 0) && ((i % 4 == 3 && i < 32) || i == 39 || i == 41)) {
      AddBound(v_new[i], 1);
    }
    if ((edge_id[4] > 0) && ((i % 16 < 4 && i < 20) || i == 36 || i == 37)) {
      AddBound(v_new[i], 1);
    }
    if ((edge_id[5] > 0) && ((i % 16 > 11 && i < 32) || i == 42 || i == 43)) {
      AddBound(v_new[i], 1);
    }
  }
  if (CheckBound_edge(v_new[0], v_new[3]) == 1) {
    AddBound(v_new[1], 1);
    AddBound(v_new[2], 1);
  }
  if (CheckBound_edge(v_new[3], v_new[15]) == 1) {
    AddBound(v_new[7], 1);
    AddBound(v_new[11], 1);
  }
  if (CheckBound_edge(v_new[15], v_new[12]) == 1) {
    AddBound(v_new[13], 1);
    AddBound(v_new[14], 1);
  }
  if (CheckBound_edge(v_new[0], v_new[12]) == 1) {
    AddBound(v_new[4], 1);
    AddBound(v_new[8], 1);
  }
  if (CheckBound_edge(v_new[0], v_new[44]) == 1) {
    AddBound(v_new[16], 1);
  }
  if (CheckBound_edge(v_new[3], v_new[45]) == 1) {
    AddBound(v_new[19], 1);
  }
  if (CheckBound_edge(v_new[12], v_new[46]) == 1) {
    AddBound(v_new[28], 1);
  }
  if (CheckBound_edge(v_new[15], v_new[47]) == 1) {
    AddBound(v_new[31], 1);
  }
}

int Geoframe::AddTri(uint v1, uint v2, uint v3) {
  if (numTris + 1 >= tsize) {
    tsize <<= 1;
    triangles = (uint(*)[3])reallocSafely(triangles, sizeof(uint[3]) * tsize);
    bound_tri = (uint(*))reallocSafely(bound_tri, sizeof(uint) * tsize);
  }
  bound_tri[numTris] = 0;
  triangles[numTris][0] = v1;
  triangles[numTris][1] = v2;
  triangles[numTris][2] = v3;
  return numTris++;
}

float Geoframe::get_aspect_ratio(uint v0, uint v1, uint v2) {
  float p, s, r_in, r_out;
  float a = 0.0;
  float b = 0.0;
  float c = 0.0;
  for (int i = 0; i < 3; i++) {
    a += (verts[v1][i] - verts[v0][i]) * (verts[v1][i] - verts[v0][i]);
    b += (verts[v2][i] - verts[v1][i]) * (verts[v2][i] - verts[v1][i]);
    c += (verts[v0][i] - verts[v2][i]) * (verts[v0][i] - verts[v2][i]);
  }
  a = (float)sqrt(a);
  b = (float)sqrt(b);
  c = (float)sqrt(c);
  p = (a + b + c) / 2;
  s = (float)sqrt(p * (p - a) * (p - b) * (p - c));
  r_in = s / p;
  r_out = a * b * c / (4 * s);
  return (r_in / r_out);
}

void Geoframe::Add_2_Tetra(uint *v, uint my_vertex) {
  float aspect_ratio_0, aspect_ratio_1, temp;
  // dectect duplicate vertices
  if (v[0] == v[1]) {
    AddTetra(v[1], v[3], v[2], my_vertex);
  } else if (v[1] == v[2]) {
    AddTetra(v[0], v[3], v[1], my_vertex);
  } else if (v[2] == v[3] || v[3] == v[0]) {
    AddTetra(v[0], v[2], v[1], my_vertex);
  } else {
    aspect_ratio_0 = get_aspect_ratio(v[0], v[2], v[1]);
    temp = get_aspect_ratio(v[0], v[3], v[2]);
    if (temp < aspect_ratio_0) {
      aspect_ratio_0 = temp;
    }
    aspect_ratio_1 = get_aspect_ratio(v[0], v[3], v[1]);
    temp = get_aspect_ratio(v[1], v[3], v[2]);
    if (temp < aspect_ratio_1) {
      aspect_ratio_1 = temp;
    }
    if (aspect_ratio_0 > aspect_ratio_1) {
      AddTetra(v[0], v[2], v[1], my_vertex);
      AddTetra(v[0], v[3], v[2], my_vertex);
    } else {
      AddTetra(v[0], v[3], v[1], my_vertex);
      AddTetra(v[1], v[3], v[2], my_vertex);
    }
  }
}

void Geoframe::Extend_Tri(uint v0, uint v1, uint v2) {
  uint vv0, vv1, vv2, t0, t1, t2, t, i, nstep;
  float v_pos[3], norm[3], r, step, radius_o, center;
  vv0 = v0;
  vv1 = v1;
  vv2 = v2;
  if (vv1 < vv0 && vv1 < vv2) {
    v0 = vv1;
    v1 = vv2;
    v2 = vv0;
  }
  if (vv2 < vv0 && vv2 < vv1) {
    v0 = vv2;
    v1 = vv0;
    v2 = vv1;
  }
  t0 = v0;
  t1 = v1;
  t2 = v2;
  center = 64;
  radius_o = 64;
  step = (30 * 20 - radius_o) / 10; // 40 times -- 10
  nstep = 5;
  if (vtx_idx_arr_extend[v0] == -1) {
    for (i = 1; i < nstep; i++) {
      r = (float)sqrt((verts[v0][0] - center) * (verts[v0][0] - center) +
                      (verts[v0][1] - center) * (verts[v0][1] - center) +
                      (verts[v0][2] - center) * (verts[v0][2] - center));
      for (t = 0; t < 3; t++) {
        v_pos[t] =
            (verts[v0][t] - center) * (radius_o + step * i * (i + 1) / 2) / r +
            center;
        norm[t] = normals[v0][t];
      }
      vv0 = AddVert(v_pos, norm);
      AddBound(vv0, 1);
      if (i == 1) {
        vtx_idx_arr_extend[v0] = vv0;
      }
    }
  }
  if (vtx_idx_arr_extend[v1] == -1) {
    for (i = 1; i < nstep; i++) {
      r = (float)sqrt((verts[v1][0] - center) * (verts[v1][0] - center) +
                      (verts[v1][1] - center) * (verts[v1][1] - center) +
                      (verts[v1][2] - center) * (verts[v1][2] - center));
      for (t = 0; t < 3; t++) {
        v_pos[t] =
            (verts[v1][t] - center) * (radius_o + step * i * (i + 1) / 2) / r +
            center;
        norm[t] = normals[v1][t];
      }
      vv1 = AddVert(v_pos, norm);
      AddBound(vv1, 1);
      if (i == 1) {
        vtx_idx_arr_extend[v1] = vv1;
      }
    }
  }
  if (vtx_idx_arr_extend[v2] == -1) {
    for (i = 1; i < nstep; i++) {
      r = (float)sqrt((verts[v2][0] - center) * (verts[v2][0] - center) +
                      (verts[v2][1] - center) * (verts[v2][1] - center) +
                      (verts[v2][2] - center) * (verts[v2][2] - center));
      for (t = 0; t < 3; t++) {
        v_pos[t] =
            (verts[v2][t] - center) * (radius_o + step * i * (i + 1) / 2) / r +
            center;
        norm[t] = normals[v2][t];
      }
      vv2 = AddVert(v_pos, norm);
      AddBound(vv2, 1);
      if (i == 1) {
        vtx_idx_arr_extend[v2] = vv2;
      }
    }
  }
  for (i = 1; i < nstep; i++) {
    vv0 = vtx_idx_arr_extend[v0] + (i - 1);
    vv1 = vtx_idx_arr_extend[v1] + (i - 1);
    vv2 = vtx_idx_arr_extend[v2] + (i - 1);
    if (v0 < v1 && v1 < v2) {
      AddTetra(vv0, vv1, vv2, t0);
      AddTetra(vv2, vv1, t1, t0);
      AddTetra(vv2, t1, t2, t0);
    } else // if(v0 < v2 && v2 < v1)
    {
      AddTetra(vv0, vv1, vv2, t0);
      AddTetra(vv2, vv1, t2, t0);
      AddTetra(vv1, t1, t2, t0);
    }
    t0 = vv0;
    t1 = vv1;
    t2 = vv2;
  }
}

void Geoframe::Extend_Tetra(uint *v) {
  float aspect_ratio_0, aspect_ratio_1, temp, radius, center;
  center = 64;
  radius = (float)sqrt((verts[v[0]][0] - center) * (verts[v[0]][0] - center) +
                       (verts[v[0]][1] - center) * (verts[v[0]][1] - center) +
                       (verts[v[0]][2] - center) * (verts[v[0]][2] - center));
  if (radius > 40.0) {
    // dectect duplicate vertices
    if (v[0] == v[1]) {
      Extend_Tri(v[1], v[3], v[2]);
    } else if (v[1] == v[2]) {
      Extend_Tri(v[0], v[3], v[1]);
    } else if (v[2] == v[3] || v[3] == v[0]) {
      Extend_Tri(v[0], v[2], v[1]);
    } else {
      aspect_ratio_0 = get_aspect_ratio(v[0], v[2], v[1]);
      temp = get_aspect_ratio(v[0], v[3], v[2]);
      if (temp < aspect_ratio_0) {
        aspect_ratio_0 = temp;
      }
      aspect_ratio_1 = get_aspect_ratio(v[0], v[3], v[1]);
      temp = get_aspect_ratio(v[1], v[3], v[2]);
      if (temp < aspect_ratio_1) {
        aspect_ratio_1 = temp;
      }
      if (aspect_ratio_0 > aspect_ratio_1) {
        Extend_Tri(v[0], v[2], v[1]);
        Extend_Tri(v[0], v[3], v[2]);
      } else {
        Extend_Tri(v[0], v[3], v[1]);
        Extend_Tri(v[1], v[3], v[2]);
      }
    }
  }
}

void Geoframe::Add_2_Tri(uint *v) {
  float aspect_ratio_0, aspect_ratio_1, temp;
  // dectect duplicate vertices
  if (v[0] == v[1]) {
    AddTri(v[1], v[2], v[3]);
  } else if (v[1] == v[2]) {
    AddTri(v[0], v[1], v[3]);
  } else if (v[2] == v[3] || v[3] == v[0]) {
    AddTri(v[0], v[1], v[2]);
  } else {
    aspect_ratio_0 = get_aspect_ratio(v[0], v[2], v[1]);
    temp = get_aspect_ratio(v[0], v[3], v[2]);
    if (temp < aspect_ratio_0) {
      aspect_ratio_0 = temp;
    }
    aspect_ratio_1 = get_aspect_ratio(v[0], v[3], v[1]);
    temp = get_aspect_ratio(v[1], v[3], v[2]);
    if (temp < aspect_ratio_1) {
      aspect_ratio_1 = temp;
    }
    if (aspect_ratio_0 > aspect_ratio_1) {
      AddTri(v[0], v[1], v[2]);
      AddTri(v[2], v[3], v[0]);
    } else {
      AddTri(v[0], v[1], v[3]);
      AddTri(v[1], v[2], v[3]);
    }
  }
}

float Geoframe::getRadius(float a[3][3], float b[3], float v0[3]) {
  float temp, x, y, z;
  if (fabs(a[0][0]) < fabs(a[1][0]) && fabs(a[2][0]) <= fabs(a[1][0])) {
    for (int i = 0; i < 3; i++) {
      temp = a[0][i];
      a[0][i] = a[1][i];
      a[1][i] = temp;
    }
    temp = b[0];
    b[0] = b[1];
    b[1] = temp;
  } else {
    if (fabs(a[0][0]) < fabs(a[2][0]) && fabs(a[1][0]) <= fabs(a[2][0])) {
      for (int i = 0; i < 3; i++) {
        temp = a[0][i];
        a[0][i] = a[2][i];
        a[2][i] = temp;
      }
      temp = b[0];
      b[0] = b[2];
      b[2] = temp;
    }
  }
  for (int i = 1; i < 3; i++) {
    if (fabs(a[i][0]) > powf(10.0, -12.0)) {
      a[i][1] = a[0][0] * a[i][1] / a[i][0] - a[0][1];
      a[i][2] = a[0][0] * a[i][2] / a[i][0] - a[0][2];
      b[i] = a[0][0] * b[i] / a[i][0] - b[0];
      a[i][0] = 0.0;
    }
  }
  if (fabs(a[1][1]) < fabs(a[2][1])) {
    for (int i = 0; i < 3; i++) {
      temp = a[1][i];
      a[1][i] = a[2][i];
      a[2][i] = temp;
    }
    temp = b[1];
    b[1] = b[2];
    b[2] = temp;
  }
  if (fabs(a[2][1]) > powf(10.0, -12.0)) {
    a[2][2] = a[1][1] * a[2][2] / a[2][1] - a[1][2];
    b[2] = a[1][1] * b[2] / a[2][1] - b[1];
    a[2][1] = 0.0;
  }
  if (fabs(a[0][0]) > powf(10.0, -12.0) && fabs(a[1][1]) > powf(10.0, -12.0) &&
      fabs(a[2][2]) > powf(10.0, -12.0)) {
    z = b[2] / a[2][2];
    y = (b[1] - a[1][2] * z) / a[1][1];
    x = (b[0] - a[0][1] * y - a[0][2] * z) / a[0][0];
  } else {
    printf("--- singular ---\n");
  }
  float delt_x = x - v0[0];
  float delt_y = y - v0[1];
  float delt_z = z - v0[2];
  float radius =
      (float)sqrt(delt_x * delt_x + delt_y * delt_y + delt_z * delt_z);
  return radius;
}

// return 1 -- bad quality;             return 0 -- good quality
int Geoframe::testTetrahedron(float v0[3], float v1[3], float v2[3],
                              float v3[3]) {
  // min/max angles
  float v01[3], v02[3], n[3], v[3], t, min, max, temp;
  float v03[3], v13[3], v23[3], v0p[3], v1p[3], v2p[3];
  for (int i = 0; i < 3; i++) {
    v01[i] = v1[i] - v0[i];
    v02[i] = v2[i] - v0[i];
  }
  n[0] = v01[1] * v02[2] - v01[2] * v02[1];
  n[1] = v01[2] * v02[0] - v01[0] * v02[2];
  n[2] = v01[0] * v02[1] - v01[1] * v02[0];
  temp = (float)sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  n[0] = -n[0] / temp;
  n[1] = -n[1] / temp;
  n[2] = -n[2] / temp;
  t = (v2[0] - v3[0]) * n[0] + (v2[1] - v3[1]) * n[1] + (v2[2] - v3[2]) * n[2];
  if (t < 0.0) {
    return 1; // RHS
  }
  for (int i = 0; i < 3; i++) {
    v[i] = v3[i] + n[i] * t;
    v03[i] = v3[i] - v0[i];
    v0p[i] = v[i] - v0[i];
    v13[i] = v3[i] - v1[i];
    v1p[i] = v[i] - v1[i];
    v23[i] = v3[i] - v2[i];
    v2p[i] = v[i] - v2[i];
  }
  float a_dot_b = v03[0] * v0p[0] + v03[1] * v0p[1] + v03[2] * v0p[2];
  float module_a =
      (float)sqrt(v03[0] * v03[0] + v03[1] * v03[1] + v03[2] * v03[2]);
  float module_b =
      (float)sqrt(v0p[0] * v0p[0] + v0p[1] * v0p[1] + v0p[2] * v0p[2]);
  float angle0;
  if (module_b < powf(10.0, -12.0)) {
    angle0 = 90;
  } else {
    angle0 = (float)acos(a_dot_b / (module_a * module_b)) * 180 / 3.1415926f;
  }
  a_dot_b = v13[0] * v1p[0] + v13[1] * v1p[1] + v13[2] * v1p[2];
  module_a = (float)sqrt(v13[0] * v13[0] + v13[1] * v13[1] + v13[2] * v13[2]);
  module_b = (float)sqrt(v1p[0] * v1p[0] + v1p[1] * v1p[1] + v1p[2] * v1p[2]);
  float angle1;
  if (module_b < powf(10.0, -12.0)) {
    angle1 = 90;
  } else {
    angle1 = (float)acos(a_dot_b / (module_a * module_b)) * 180 / 3.1415926f;
  }
  a_dot_b = v23[0] * v2p[0] + v23[1] * v2p[1] + v23[2] * v2p[2];
  module_a = (float)sqrt(v23[0] * v23[0] + v23[1] * v23[1] + v23[2] * v23[2]);
  module_b = (float)sqrt(v2p[0] * v2p[0] + v2p[1] * v2p[1] + v2p[2] * v2p[2]);
  float angle2;
  if (module_b < powf(10.0, -12.0)) {
    angle2 = 90;
  } else {
    angle2 = (float)acos(a_dot_b / (module_a * module_b)) * 180 / 3.1415926f;
  }
  min = angle0;
  max = angle0;
  if (angle1 < min) {
    min = angle1;
  } else {
    max = angle1;
  }
  if (angle2 < min) {
    min = angle2;
  }
  if (angle2 > max) {
    max = angle2;
  }
  if ((min < 10.0) || (max > 160.0)) {
    return 1;
  }
  // tetrahedral quality measure = volume of tetrahedron / volume of equilateral
  // tetrahedron with same circumsphere radius
  float h = (float)sqrt((v[0] - v3[0]) * (v[0] - v3[0]) +
                        (v[1] - v3[1]) * (v[1] - v3[1]) +
                        (v[2] - v3[2]) * (v[2] - v3[2]));
  float side_a = (float)sqrt((v1[0] - v0[0]) * (v1[0] - v0[0]) +
                             (v1[1] - v0[1]) * (v1[1] - v0[1]) +
                             (v1[2] - v0[2]) * (v1[2] - v0[2]));
  float side_b = (float)sqrt((v2[0] - v1[0]) * (v2[0] - v1[0]) +
                             (v2[1] - v1[1]) * (v2[1] - v1[1]) +
                             (v2[2] - v1[2]) * (v2[2] - v1[2]));
  float side_c = (float)sqrt((v2[0] - v0[0]) * (v2[0] - v0[0]) +
                             (v2[1] - v0[1]) * (v2[1] - v0[1]) +
                             (v2[2] - v0[2]) * (v2[2] - v0[2]));
  float periphery = side_a + side_b + side_c;
  float area =
      0.25f * (float)sqrt(periphery * (periphery - 2 * side_a) *
                          (periphery - 2 * side_b) * (periphery - 2 * side_c));
  float volume_1 = area * h / 3;
  float matrix_lhs[3][3], vector_rhs[3], radius;
  for (int i = 0; i < 3; i++) {
    matrix_lhs[0][i] = v1[i] - v0[i];
    matrix_lhs[1][i] = v2[i] - v0[i];
    matrix_lhs[2][i] = v3[i] - v0[i];
  }
  vector_rhs[0] = (v1[0] * v1[0] - v0[0] * v0[0] + v1[1] * v1[1] -
                   v0[1] * v0[1] + v1[2] * v1[2] - v0[2] * v0[2]) /
                  2;
  vector_rhs[1] = (v2[0] * v2[0] - v0[0] * v0[0] + v2[1] * v2[1] -
                   v0[1] * v0[1] + v2[2] * v2[2] - v0[2] * v0[2]) /
                  2;
  vector_rhs[2] = (v3[0] * v3[0] - v0[0] * v0[0] + v3[1] * v3[1] -
                   v0[1] * v0[1] + v3[2] * v3[2] - v0[2] * v0[2]) /
                  2;
  if (volume_1 > powf(10.0, -12.0)) {
    radius = getRadius(matrix_lhs, vector_rhs, v0);
  } else {
    return 1;
  }
  float volume_2 = 0.25f * (float)sqrt(6) * radius * radius * radius;
  if (volume_1 / volume_2 < 0.02) {
    return 1;
  } else {
    return 0;
  }
}

// return v1
double Geoframe::testTetrahedron1(float v0[3], float v1[3], float v2[3],
                                  float v3[3]) {
  // min/max angles
  double v01[3], v02[3], n[3], v[3], t, min, max, temp;
  double v03[3], v13[3], v23[3], v0p[3], v1p[3], v2p[3];
  for (int i = 0; i < 3; i++) {
    v01[i] = v1[i] - v0[i];
    v02[i] = v2[i] - v0[i];
  }
  n[0] = v01[1] * v02[2] - v01[2] * v02[1];
  n[1] = v01[2] * v02[0] - v01[0] * v02[2];
  n[2] = v01[0] * v02[1] - v01[1] * v02[0];
  temp = (float)sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  n[0] = -n[0] / temp;
  n[1] = -n[1] / temp;
  n[2] = -n[2] / temp;
  t = (v2[0] - v3[0]) * n[0] + (v2[1] - v3[1]) * n[1] + (v2[2] - v3[2]) * n[2];
  for (int i = 0; i < 3; i++) {
    v[i] = v3[i] + n[i] * t;
    v03[i] = v3[i] - v0[i];
    v0p[i] = v[i] - v0[i];
    v13[i] = v3[i] - v1[i];
    v1p[i] = v[i] - v1[i];
    v23[i] = v3[i] - v2[i];
    v2p[i] = v[i] - v2[i];
  }
  float a_dot_b = (float)(v03[0] * v0p[0] + v03[1] * v0p[1] + v03[2] * v0p[2]);
  float module_a =
      (float)sqrt(v03[0] * v03[0] + v03[1] * v03[1] + v03[2] * v03[2]);
  float module_b =
      (float)sqrt(v0p[0] * v0p[0] + v0p[1] * v0p[1] + v0p[2] * v0p[2]);
  float angle0;
  if (module_b < powf(10.0, -12.0)) {
    angle0 = 90;
  } else {
    angle0 = (float)acos(a_dot_b / (module_a * module_b)) * 180 / 3.1415926f;
  }
  a_dot_b = (float)(v13[0] * v1p[0] + v13[1] * v1p[1] + v13[2] * v1p[2]);
  module_a = (float)sqrt(v13[0] * v13[0] + v13[1] * v13[1] + v13[2] * v13[2]);
  module_b = (float)sqrt(v1p[0] * v1p[0] + v1p[1] * v1p[1] + v1p[2] * v1p[2]);
  float angle1;
  if (module_b < powf(10.0, -12.0)) {
    angle1 = 90;
  } else {
    angle1 = (float)acos(a_dot_b / (module_a * module_b)) * 180 / 3.1415926f;
  }
  a_dot_b = (float)(v23[0] * v2p[0] + v23[1] * v2p[1] + v23[2] * v2p[2]);
  module_a = (float)sqrt(v23[0] * v23[0] + v23[1] * v23[1] + v23[2] * v23[2]);
  module_b = (float)sqrt(v2p[0] * v2p[0] + v2p[1] * v2p[1] + v2p[2] * v2p[2]);
  float angle2;
  if (module_b < powf(10.0, -12.0)) {
    angle2 = 90;
  } else {
    angle2 = (float)acos(a_dot_b / (module_a * module_b)) * 180 / 3.1415926f;
  }
  min = angle0;
  max = angle0;
  if (angle1 < min) {
    min = angle1;
  } else {
    max = angle1;
  }
  if (angle2 < min) {
    min = angle2;
  }
  if (angle2 > max) {
    max = angle2;
  }
  // if((min < 10.0) || (max > 160.0)) return 1;
  //  tetrahedral quality measure = volume of tetrahedron / volume of
  //  equilateral tetrahedron with same circumsphere radius
  double h = 0.0;
  h = (double)sqrt((v[0] - v3[0]) * (v[0] - v3[0]) +
                   (v[1] - v3[1]) * (v[1] - v3[1]) +
                   (v[2] - v3[2]) * (v[2] - v3[2]));
  double side_a = (double)sqrt((v1[0] - v0[0]) * (v1[0] - v0[0]) +
                               (v1[1] - v0[1]) * (v1[1] - v0[1]) +
                               (v1[2] - v0[2]) * (v1[2] - v0[2]));
  double side_b = (double)sqrt((v2[0] - v1[0]) * (v2[0] - v1[0]) +
                               (v2[1] - v1[1]) * (v2[1] - v1[1]) +
                               (v2[2] - v1[2]) * (v2[2] - v1[2]));
  double side_c = (double)sqrt((v2[0] - v0[0]) * (v2[0] - v0[0]) +
                               (v2[1] - v0[1]) * (v2[1] - v0[1]) +
                               (v2[2] - v0[2]) * (v2[2] - v0[2]));
  double periphery = side_a + side_b + side_c;
  double area =
      0.25f * (double)sqrt(periphery * (periphery - 2 * side_a) *
                           (periphery - 2 * side_b) * (periphery - 2 * side_c));
  double volume_1 = area * t / 3;
  double matrix_lhs[3][3], vector_rhs[3];
  for (int i = 0; i < 3; i++) {
    matrix_lhs[0][i] = v1[i] - v0[i];
    matrix_lhs[1][i] = v2[i] - v0[i];
    matrix_lhs[2][i] = v3[i] - v0[i];
  }
  vector_rhs[0] = (v1[0] * v1[0] - v0[0] * v0[0] + v1[1] * v1[1] -
                   v0[1] * v0[1] + v1[2] * v1[2] - v0[2] * v0[2]) /
                  2;
  vector_rhs[1] = (v2[0] * v2[0] - v0[0] * v0[0] + v2[1] * v2[1] -
                   v0[1] * v0[1] + v2[2] * v2[2] - v0[2] * v0[2]) /
                  2;
  vector_rhs[2] = (v3[0] * v3[0] - v0[0] * v0[0] + v3[1] * v3[1] -
                   v0[1] * v0[1] + v3[2] * v3[2] - v0[2] * v0[2]) /
                  2;
  return volume_1;
}

void Geoframe::edge_contraction_tri(uint *v) {
  float aspect_ratio, a, b, c;
  aspect_ratio = get_aspect_ratio(v[0], v[1], v[2]);
  if (aspect_ratio < 0.1) {
    a = 0.0;
    b = 0.0;
    c = 0.0;
    for (int i = 0; i < 3; i++) {
      a +=
          (verts[v[1]][i] - verts[v[0]][i]) * (verts[v[1]][i] - verts[v[0]][i]);
      b +=
          (verts[v[2]][i] - verts[v[1]][i]) * (verts[v[2]][i] - verts[v[1]][i]);
      c +=
          (verts[v[0]][i] - verts[v[2]][i]) * (verts[v[0]][i] - verts[v[2]][i]);
    }
    if (a <= b && a <= c) {
      for (int i = 0; i < 3; i++) {
        verts[v[1]][i] = verts[v[0]][i];
      }
    } else if (b <= a && b <= c) {
      for (int i = 0; i < 3; i++) {
        verts[v[2]][i] = verts[v[1]][i];
      }
    } else {
      for (int i = 0; i < 3; i++) {
        verts[v[0]][i] = verts[v[2]][i];
      }
    }
  }
}

void Geoframe::edge_contraction_tetra(uint *v, int num) {
  float min, max, e[6];
  float v0[3], v1[3], v2[3], v3[3];
  int i, id_min, sign;
  int dummy = 0;
  dummy = num;
  for (int i = 0; i < 6; i++) {
    e[i] = 0.0;
  }
  for (int i = 0; i < 3; i++) {
    e[0] += (verts[v[1]][i] - verts[v[0]][i]) *
            (verts[v[1]][i] - verts[v[0]][i]); // e_01
    e[1] += (verts[v[2]][i] - verts[v[1]][i]) *
            (verts[v[2]][i] - verts[v[1]][i]); // e_12
    e[2] += (verts[v[0]][i] - verts[v[2]][i]) *
            (verts[v[0]][i] - verts[v[2]][i]); // e_20
    e[3] += (verts[v[3]][i] - verts[v[0]][i]) *
            (verts[v[3]][i] - verts[v[0]][i]); // e_03
    e[4] += (verts[v[3]][i] - verts[v[1]][i]) *
            (verts[v[3]][i] - verts[v[1]][i]); // e_13
    e[5] += (verts[v[3]][i] - verts[v[2]][i]) *
            (verts[v[3]][i] - verts[v[2]][i]); // e_23
    v0[i] = verts[v[0]][i];
    v1[i] = verts[v[1]][i];
    v2[i] = verts[v[2]][i];
    v3[i] = verts[v[3]][i];
  }
  min = e[0];
  max = e[0];
  id_min = 0;
  for (i = 1; i < 6; i++) {
    if (e[i] < min) {
      min = e[i];
      id_min = i;
    }
    if (e[i] > max) {
      max = e[i];
    }
  }
  sign = testTetrahedron(v0, v1, v2, v3);
  if (sign == 1) {
    if (id_min == 0) {
      for (int i = 0; i < 3; i++) {
        verts[v[1]][i] = verts[v[0]][i];
      }
    } else if (id_min == 1) {
      for (int i = 0; i < 3; i++) {
        verts[v[2]][i] = verts[v[1]][i];
      }
    } else if (id_min == 2) {
      for (int i = 0; i < 3; i++) {
        verts[v[2]][i] = verts[v[0]][i];
      }
    } else if (id_min == 3) {
      for (int i = 0; i < 3; i++) {
        verts[v[3]][i] = verts[v[0]][i];
      }
    } else if (id_min == 4) {
      for (int i = 0; i < 3; i++) {
        verts[v[3]][i] = verts[v[1]][i];
      }
    } else // id_min = 5
    {
      for (int i = 0; i < 3; i++) {
        verts[v[3]][i] = verts[v[2]][i];
      }
    }
  }
}

void Geoframe::edge_contraction(uint *v, int num) {
  if (num == 3) {
    edge_contraction_tri(v);
  } else // num == 4 or 5
  {
    edge_contraction_tetra(v, num);
  }
}

void Geoframe::AddPyramid(uint *v, uint my_vertex, int num) {
  float vv0[3], vv1[3], vv2[3], vv3[3];
  AddQuad(v, num);
  for (int t = 0; t < 3; t++) {
    vv0[t] = verts[v[0]][t];
    vv1[t] = verts[v[1]][t];
    vv2[t] = verts[v[2]][t];
    vv3[t] = verts[my_vertex][t];
  }
  int sign = testRHS(vv0, vv1, vv2, vv3);
  if (sign == 0) {
    AddTri(v[0], my_vertex, v[1]);
    AddTri(v[1], my_vertex, v[2]);
  } else {
    AddTri(v[1], my_vertex, v[0]);
    AddTri(v[2], my_vertex, v[1]);
  }
  for (int t = 0; t < 3; t++) {
    vv0[t] = verts[v[0]][t];
    vv1[t] = verts[v[2]][t];
    vv2[t] = verts[v[3]][t];
    vv3[t] = verts[my_vertex][t];
  }
  sign = testRHS(vv0, vv1, vv2, vv3);
  if (sign == 0) {
    AddTri(v[2], my_vertex, v[3]);
    AddTri(v[3], my_vertex, v[0]);
  } else {
    AddTri(v[3], my_vertex, v[2]);
    AddTri(v[0], my_vertex, v[3]);
  }
}

// RHS -- 0;    LHS -- 1;
int Geoframe::testRHS(float v0[3], float v1[3], float v2[3], float v3[3]) {
  float v01[3], v02[3], v03[3], n[3];
  for (int i = 0; i < 3; i++) {
    v01[i] = v1[i] - v0[i];
    v02[i] = v2[i] - v0[i];
    v03[i] = v3[i] - v0[i];
  }
  n[0] = v01[1] * v02[2] - v01[2] * v02[1];
  n[1] = v01[2] * v02[0] - v01[0] * v02[2];
  n[2] = v01[0] * v02[1] - v01[1] * v02[0];
  float sign = v03[0] * n[0] + v03[1] * n[1] + v03[2] * n[2];
  if (sign < 0.0) {
    return 1;
  } else if (sign == 0.0) {
    return -1;
  } else {
    return 0;
  }
}

void Geoframe::AddTetra(uint v1, uint v2, uint v3, uint v4) {
  float vv0[3], vv1[3], vv2[3], vv3[3];
  for (int t = 0; t < 3; t++) {
    vv0[t] = verts[v1][t];
    vv1[t] = verts[v2][t];
    vv2[t] = verts[v3][t];
    vv3[t] = verts[v4][t];
  }
  int sign = testRHS(vv0, vv1, vv2, vv3);
  if (sign == 1) {
    AddTri(v1, v3, v2);
    AddTri(v2, v3, v4);
    AddTri(v1, v4, v3);
    AddTri(v1, v2, v4);
  }
  if (sign == 0) {
    AddTri(v1, v2, v3);
    bound_tri[numTris - 1] = 1;
    AddTri(v3, v2, v4);
    AddTri(v1, v3, v4);
    AddTri(v1, v4, v2);
  }
}

int Geoframe::AddVert(float v_pos[3], float norm[3]) {
  // cout << "lbieErrorBegin"<< endl; // Craig: segfaults here
  // cout << "bound_sign: " << bound_sign << endl; // allocated somewhere
  // cout << "bound_sign[0]: " << bound_sign[0] << endl; // can't dereference
  if (numVerts + 1 > vsize) {
    cout << "Inside If statement" << endl;
    vsize <<= 1;
    verts = (float(*)[3])reallocSafely(verts, sizeof(float[3]) * vsize);
    funcs = (float(*)[1])reallocSafely(funcs, sizeof(float[1]) * vsize);
    normals = (float(*)[3])reallocSafely(normals, sizeof(float[3]) * vsize);
    curvatures = (float(*)[2])reallocSafely(normals, sizeof(float[2]) * vsize);
    bound_sign = (uint(*))reallocSafely(bound_sign, sizeof(uint) * vsize);
    vtxnew_sign = (uint(*))reallocSafely(vtxnew_sign, sizeof(uint) * vsize);
    bound_edge =
        (uint(*)[18])reallocSafely(bound_edge, sizeof(uint[18]) * vsize);
    refine_edge =
        (int(*)[18])reallocSafely(refine_edge, sizeof(int[18]) * vsize);
    refine_edgevtx =
        (uint(*)[18])reallocSafely(refine_edgevtx, sizeof(uint[18]) * vsize);
  }
  // cout << "vsize: " << vsize << endl; // it's 3000000
  // cout << "numVerts: " << numVerts << endl; // it's 0
  // cout << "bound_sign: " << bound_sign << endl; // reallocated somewhere
  // cout << "*bound_sign: " << *bound_sign << endl; // can't dereference
  // cout << "bound_sign[0]: " << bound_sign[0] << endl; // can't dereference
  bound_sign[numVerts] = 0; // here is the actual error
  // cout << "lbieErrorEnd " << endl;
  vtxnew_sign[numVerts] = 0;
  for (int i = 0; i < 18; i++) {

    bound_edge[numVerts][i] = 0;
  }
  for (int i = 0; i < 18; i++) {
    refine_edge[numVerts][i] = -1;
  }

  for (int i = 0; i < 18; i++) {
    refine_edgevtx[numVerts][i] = 0;
  }
  for (int i = 0; i < 3; i++) {
    verts[numVerts][i] = v_pos[i];
  }
  for (int i = 0; i < 3; i++) {
    normals[numVerts][i] = norm[i];
  }
  for (int i = 0; i < 2; i++) {
    curvatures[numVerts][i] = 0;
  }
  return numVerts++;
}

// sign = 1 -- outer    // sign = -1 -- inside
// sign = 0 -- not on boundary
void Geoframe::AddBound(int index, int sign) { bound_sign[index] = sign; }

void Geoframe::AddBound_edge(uint index_1, uint index_2) {
  if (index_2 < index_1) {
    int temp = index_1;
    index_1 = index_2;
    index_2 = temp;
  }
  for (int i = 0; i < 18; i++) {
    if (bound_edge[index_1][i] == index_2) {
      break;
    }
    if (bound_edge[index_1][i] == 0) {
      bound_edge[index_1][i] = index_2;
    }
  }
}

int Geoframe::CheckBound_edge(uint index_1, uint index_2) {
  if (index_2 < index_1) {
    int temp = index_1;
    index_1 = index_2;
    index_2 = temp;
  }
  for (int i = 0; i < 18; i++) {
    if (bound_edge[index_1][i] == index_2) {
      return 1;
    }
  }
  return 0;
}

void Geoframe::AddVtxNew(int index, int sign) { vtxnew_sign[index] = sign; }

int Geoframe::center_vtx(int v1, int v2, int v3) {
  float center_vtx[3], norm[3];
  for (int i = 0; i < 3; i++) {
    center_vtx[i] = (verts[v1][i] + verts[v2][i] + verts[v3][i]) / 3;
    // grad addition
    norm[i] = (normals[v1][i] + normals[v2][i] + normals[v3][i]) / 3;
  }
  return AddVert(center_vtx, norm);
}

void Geoframe::setSpan(float span0, float span1, float span2) {
  span[0] = span0;
  span[1] = span1;
  span[2] = span2;
}

void Geoframe::setMin(float min_x_, float min_y_, float min_z_) {
  min_x = min_x_;
  min_y = min_y_;
  min_z = min_z_;
}
