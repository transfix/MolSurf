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
#include <Utility/utility.h>
#include <HLevelSet/HLevelSet.h>
#include <SimpleVolumeData/SimpleVolumeData.h>
#include <PDBParser/Atom.h>
#include <PDBParser/FlattenGOA.h>
// #include <Blurmaps/BlurMapsDataManager.h>
#include <LBIE_lib/normalspline.h>
#include <boost/tuple/tuple.hpp>

void ConvertToInterpolationCoefficients_Qu(float *, int, float *, int, float);
void ConvertToInterpolationCoefficients_1D(float *, int, float);
void ConvertToInterpolationCoefficients_2D(float *, int, int, float);
void ConvertToInterpolationCoefficients_3D_Xu(float *, int, int, int, float);
void ConvertToInterpolationCoefficients_3D(float *, int, int, int, float);
void Tensor_333(float *, float *, float *, float *);
void Tensor_444(float *, float *, float *, float *);
void Take_27_Coefficients(float *, int, int, int, int, int, int, float *);
void Take_64_Coefficients(float *, int, int, int, int, int, int, float *);
void Evaluat_Four_Basis(float x, float *values);
void Evaluate_Four_Basis_First_Partial(float x, float *firstPartials);
void Evaluate_Four_Basis_Second_Partial(float x, float *firstPartials);
void EvaluateCubicSplineOrder1PartialsAtGridPoint(float *, float, float, float,
                                                  int, int, int, int, int, int,
                                                  float *);
void EvaluateCubicSplineOrder2PartialsAtGridPoint(float *, float, float, float,
                                                  int, int, int, int, int, int,
                                                  float *);
void Divided_DifferenceOrder2PartialsAtGridPoint(float *, float, float, float,
                                                 int, int, int, int, int, int,
                                                 float *);
void InitialData_Sphere(float *, float, float, float, float, float, float, int,
                        int, int);
void InitialData_Box(float *, float, float, float, float, float, float, int,
                     int, int);
void MeanCurvatureFlow(float *, float *, float *, float, float, float, float,
                       int, int, int, int);
void Constraint_MeanCurvatureFlow(float *, float *, float *, float *, float,
                                  float, float, float, float, int, int, int,
                                  int);
void ComputeTensorXYZ();
void ReInitilazation(float *, float *, float, float, float, int, int, int,
                     float *);
void ReInitilazation_Upwind_Eno_Engquist(float *, float *, float *, float,
                                         float, float, int, int, int);
void ReInitilazation_Upwind_Eno_Godunov(float *, float *, float *, float, float,
                                        float, int, int, int, int);
void ReInitilazation_Upwind_Eno_Godunov_Xu(float *, float *, float *, float,
                                           float, float, int, int, int, int);
void inline DiviededDifferencing_3j(float *, float *, float *, float *, float);
void inline DiviededDifferencing_2j(float *, float *, float *, float);
void inline DiviededDifferencing_1j(float *, float *, float);
void Cubic_Eno_Interpolation(float *f0, float *u_plus, float *u_minus,
                             float dx);
void Quadr_Eno_Interpolation(float *f0, float *u_plus, float *u_minus,
                             float dx);
void Linear_Eno_Interpolation(float *f0, float *u_plus, float *u_minus,
                              float dx);
void Get_Seven_Function_Values_X(float *coeff, int nx, int ny, int nz, int i,
                                 int j, int k, float *fx);
void Get_Seven_Function_Values_Y(float *coeff, int nx, int ny, int nz, int i,
                                 int j, int k, float *fy);
void Get_Seven_Function_Values_Z(float *coeff, int nx, int ny, int nz, int i,
                                 int j, int k, float *fz);
void ReSamplingCubicSpline(float *coeff, int nx, int ny, int nz,
                           float *funvalues, int Nx, int Ny, int Nz);
void ReSamplingCubicSpline1OrderPartials(float *coeff, int nx, int ny, int nz,
                                         float *dxyz, float *MinExt,
                                         float *firstParValue, float *p);
void ReSamplingCubicSpline1OrderPartials2(float *coeff, int nx, int ny, int nz,
                                          float *dxyz, float *minExtent,
                                          float *firstParValue, float *p);
void Fast_Evaluate_Gauss_Map_One_Atom(float *corner, float *dxyz, float *center,
                                      float radius, unsigned int *nxyz, float C,
                                      float *data, float epsilon);
void HLevel_set();
float TakeACoefficient_Fast(float *c, int nx, int ny, int nz, int u, int v,
                            int w);
float TakeACoefficient_Slow(float *c, int nx, int ny, int nz, int u, int v,
                            int w);
float EvaluateCubicSplineAtGridPoint(float *, int, int, int, int, int, int);
float EvaluateCubicSplineAtAnyGivenPoint(float *, int, int, int, float, float,
                                         float);
float Extreme_Positive(float a, float b);
float Extreme_Negative(float a, float b);
float Gradient_2(float *fun, float dx, float dy, float dz, int i, int j, int k,
                 int nx, int ny, int nz);
void getBoundingBox(const std::vector<float> Vertices,
                    const std::vector<float> Radius, float min[3], float max[3],
                    int probeRadius, double blobbiness);
static float TensorF[27], TensorFx[27], TensorFy[27], TensorFz[27],
    TensorFxx[27], TensorFxy[27], TensorFxz[27], TensorFyy[27], TensorFyz[27],
    TensorFzz[27];
static float Height0, Height, Height1;
float TotalTime;
float IsoValueofHLS;
float weight; // arand - added these for TexMol version..
float iterationNum;

#define CVC_DBL_EPSILON 1.0e-9f
#define XUGUO_Z1 sqrt(3.0) - 2.0
#define OneSix3 1.0 / 216.0
#define OneSix2 2.0 / 108.0
#define OneSix1 4.0 / 54.0
#define OneSix0 8.0 / 27.0

using namespace PDBParser;

HLevelSet::HLevelSet() {}

HLevelSet::~HLevelSet() {}

bool HLevelSet::computeFunction(const std::vector<float> Vertices,
                                const std::vector<float> Radius,
                                float *func_phi, unsigned int *dim,
                                float *minExtent, float *maxExtent,
                                int nIterations, double weight) {
  float center[3], corner[3], radius0;
  float C, epsilon;
  float dx, dy, dz;
  float minx, maxx, miny, maxy, minz, maxz;
  float dt, minr, maxr, grad, width, cx, cy, cz, rx, ry, rz;
  int i, k, nx, ny, nz, t, end, size;
  float *vertices, *radius;
  vertices = new float[Vertices.size()];
  radius = new float[Radius.size()];
  size = Radius.size();
  memcpy(vertices, &(Vertices[0]), Vertices.size() * sizeof(float));
  memcpy(radius, &(Radius[0]), Radius.size() * sizeof(float));
  float *func_h = new float[dim[0] * dim[1] * dim[2]];
  float *coefficent = new float[dim[0] * dim[1] * dim[2]];
  float *funcvalue_bak = new float[dim[0] * dim[1] * dim[2]];
  if (!func_phi || !dim || !minExtent || !maxExtent) {
    return false;
  }
  // currently, at each atom location, I place a Cube ...
  // you can replace this part with your own function call
  // zero out data
  for (i = 0; i < dim[0] * dim[1] * dim[2]; i++) {
    func_h[i] = 1.0;
  }
  // loop through the atoms
  maxr = -10.0;
  minr = 10.0;
  // define a box
  minx = 10000.0;
  miny = 10000.0;
  minz = 10000.0;
  maxx = -10000.0;
  maxy = -10000.0;
  maxz = -10000.0;
  // size = atomList.size();
  for (i = 0; i < size; i++) {
    center[0] = vertices[3 * i + 0];
    center[1] = vertices[3 * i + 1];
    center[2] = vertices[3 * i + 2];
    if (minx > center[0]) {
      minx = center[0];
    }
    if (miny > center[1]) {
      miny = center[1];
    }
    if (minz > center[2]) {
      minz = center[2];
    }
    if (maxx < center[0]) {
      maxx = center[0];
    }
    if (maxy < center[1]) {
      maxy = center[1];
    }
    if (maxz < center[2]) {
      maxz = center[2];
    }
    radius0 = radius[i];
    if (radius0 > maxr) {
      maxr = radius0;
    }
    if (radius0 < minr) {
      minr = radius0;
    }
  }
  minx = minx - 5.0;
  miny = miny - 5.0;
  minz = minz - 5.0;
  maxx = maxx + 5.0;
  maxy = maxy + 5.0;
  maxz = maxz + 5.0;
  width = maxx - minx;
  if (width < maxy - miny) {
    width = maxy - miny;
  }
  if (width < maxz - minz) {
    width = maxz - minz;
  }
  // printf("Externt = %f,%f,%f,%f,%f,%f\n", minx, miny, minz, maxx, maxy,
  // maxz);
  MinExt[0] = minx;
  MinExt[1] = miny;
  MinExt[2] = minz; // Added new by zq on Dec. 23
  cx = 0.5 * (minx + maxx);
  cy = 0.5 * (miny + maxy);
  cz = 0.5 * (minz + maxz);
  rx = width / (maxx - minx);
  ry = width / (maxy - miny);
  rz = width / (maxz - minz);
  minx = rx * (minx - cx) + cx;
  maxx = rx * (maxx - cx) + cx;
  miny = ry * (miny - cy) + cy;
  maxy = ry * (maxy - cy) + cy;
  minz = rz * (minz - cz) + cz;
  maxz = rz * (maxz - cz) + cz;
  minExtent[0] = minx;
  maxExtent[0] = maxx;
  minExtent[1] = miny;
  maxExtent[1] = maxy;
  minExtent[2] = minz;
  maxExtent[2] = maxz;
  // printf("Externt = %f,%f,%f,%f,%f,%f\n", minx, miny, minz, maxx, maxy,
  // maxz); printf("Externt = %f,%f,%f\n", maxx - minx, maxy - miny, maxz -
  // minz); printf("maxr, minr = %f, %f\n", maxr, minr);
  C = (1.0 / 0.05) * log(3.0) / (1.4 + minr);
  // Note that my x, y, z is z, y, x of the input.
  corner[2] = minExtent[0];
  corner[1] = minExtent[1];
  corner[0] = minExtent[2];
  dz = (maxExtent[0] - minExtent[0]) / ((double)(dim[0] - 1));
  dy = (maxExtent[1] - minExtent[1]) / ((double)(dim[1] - 1));
  dx = (maxExtent[2] - minExtent[2]) / ((double)(dim[2] - 1));
  // printf("coner = %f, %f, %f, C = %f\n", corner[0], corner[1], corner[2], C);
  dxyz[0] = dx;
  dxyz[1] = dy;
  dxyz[2] = dz;
  // loop through the atoms
  epsilon = 1.0e-5;
  for (i = 0; i < size; i++) {
    center[0] = vertices[3 * i + 2];
    center[1] = vertices[3 * i + 1];
    center[2] = vertices[3 * i + 0];
    radius0 = radius[i];
    radius0 = radius0 + 1.4; // add water molecule radius
    Fast_Evaluate_Gauss_Map_One_Atom(corner, dxyz, center, radius0, dim, C,
                                     func_h, epsilon);
  }
  // Cut
  // initial phi, take sqrt of func_h
  grad = 2 * C * maxr;
  grad = 1.0 / grad;
  for (i = 0; i < dim[0] * dim[1] * dim[2]; i++) {
    func_phi[i] = func_h[i] * grad;
    if (func_phi[i] < -2.0) {
      func_phi[i] = -2.0;
    }
    if (func_phi[i] > 2.0) {
      func_phi[i] = 2.0;
    }
    if (func_h[i] > 1.0) {
      func_h[i] = 1.0;
    }
    if (func_h[i] < -1.0) {
      func_h[i] = -1.0;
    }
  }
  nx = dim[0];
  ny = dim[1];
  nz = dim[2];
  ComputeTensorXYZ();
  k = log(dx * dx * 0.005 / 1.732) / log(0.26795);
  k = k + 2;
  // printf("k = %d\n", k);
  Height = 3.0 + 1.4;
  Height0 = Height / 2.0;
  Height1 = 3.0;
  // printf("Height0, Height, Height1 = %f, %f, %f, dx = %f, %f\n", Height0,
  // Height,Height1, dx, pow(88.0, 3));
  //  Compute initial data
  // InitialData_Sphere(funcvalue, minx, maxx, miny, maxy, minz, maxz, nx, ny,
  // nz); InitialData_Box(funcvalue, minx, maxx, miny, maxy, minz, maxz, nx, ny,
  // nz);
  //  Convert to spline function
  for (i = 0; i < nx * ny * nz; i++) {
    coefficent[i] = func_phi[i];
  }
  ConvertToInterpolationCoefficients_3D(coefficent, nx, ny, nz,
                                        CVC_DBL_EPSILON);
  ConvertToInterpolationCoefficients_3D(func_h, nx, ny, nz, CVC_DBL_EPSILON);
  // Re-initialization
  ReInitilazation_Upwind_Eno_Godunov_Xu(func_phi, funcvalue_bak, coefficent, dx,
                                        dy, dz, nx, ny, nz, t);
  TotalTime = 0.0;
  dt = 0.01;
  dt = dx * dx;
  end = nIterations;
  for (t = 0; t < end; t++) {
    Constraint_MeanCurvatureFlow(func_h, func_phi, funcvalue_bak, coefficent,
                                 weight, dt, dx, dy, dz, nx, ny, nz, t);
    ReInitilazation_Upwind_Eno_Godunov_Xu(func_phi, funcvalue_bak, coefficent,
                                          dx, dy, dz, nx, ny, nz, t);
  }
  // MeanCurvatureFlow(func_phi, funcvalue_bak, coefficent, dt, dx, dy, dz, nx,
  // ny, nz, t);
  //  ReInitilazation(func_phi, coefficent, dx, dy, dz, nx, ny, nz, minExtent);
  //  transform the function valueto [0, 1] range for display
  float minf, maxf, w;
  minf = 1000.0;
  maxf = -1000.0;
  for (i = 0; i < nx * ny * nz; i++) {
    if (func_phi[i] > maxf) {
      maxf = func_phi[i];
    }
    if (func_phi[i] < minf) {
      minf = func_phi[i];
    }
  }
  // printf("minf, maxf = %f, %f, Toltal Time = %f, -1.4 = %f, 0 = %f\n",  minf,
  // maxf, TotalTime, (-1.4-minf)/(maxf - minf), -minf/(maxf - minf));
  IsoValueofHLS = -1.4;
  maxf = 1.0 / (maxf - minf);
  coeff = new float[nx * ny * nz];
  for (int i = 0; i < nx * ny * nz; i++) {
    coeff[i] = coefficent[i];
  }
  // Test convert
  ReSamplingCubicSpline(coefficent, nx, ny, nz, funcvalue_bak, nx, ny, nz);
  maxf = -1000.0;
  minf = 1000.0;
  int ii, j;
  for (i = 0; i < nx * ny * nz; i++) {
    w = fabs(func_phi[i] - funcvalue_bak[i]);
    if (w > maxf) {
      maxf = w;
    }
    if (w < minf) {
      minf = w;
    }
    j = i / nx;
    k = i - j * nx;
    ii = j / nx;
    j = j - ii * nx;
  }
  // printf("Maximal error = %e, minimal error = %e\n", maxf, minf);
  delete[] func_h;
  delete[] coefficent;
  delete[] funcvalue_bak;
  return true;
}

// arand: added this function to reverse the normals of the mesh
// and fix triangle orientation to be consistent with this
bool HLevelSet::fixOrientation(Geometry *geo) {
  for (int i = 0; i < geo->m_NumTris; i++) {
    int tmp = geo->m_Tris[3 * i];
    geo->m_Tris[3 * i] = geo->m_Tris[3 * i + 1];
    geo->m_Tris[3 * i + 1] = tmp;
  }
  for (int i = 0; i < geo->m_NumTriVerts; i++) {
    geo->m_TriVertNormals[i * 3] = -1.0 * geo->m_TriVertNormals[i * 3];
    geo->m_TriVertNormals[i * 3 + 1] = -1.0 * geo->m_TriVertNormals[i * 3 + 1];
    geo->m_TriVertNormals[i * 3 + 2] = -1.0 * geo->m_TriVertNormals[i * 3 + 2];
  }
}

bool HLevelSet::computeGeometryNormal(Geometry *geo, unsigned int *dim) {
  computeGeometryDerivative(geo, coeff, dim);
}

bool HLevelSet::computeGeometryError(Geometry *geo, unsigned int *dim) {
  computeGeometryDerivative(geo, coeff, dim);
}

void HLevelSet::computeGeometryDerivative(Geometry *geo, float *coeff,
                                          unsigned int *dim) {
  unsigned int nx, ny, nz, i, j, k;
  nx = dim[0];
  ny = dim[1];
  nz = dim[2];
  float H, K;
  for (int i = 0; i < geo->m_NumTriVerts; i++) {
    int u, v, w;
    float *p = geo->m_TriVerts + i * 3;
    float q[8][3];
    float partials[10];
    float partialsOld[4];
    getLowerCorner(p, u, v, w);
    EvaluateCubicSplineOrder1PartialsAtAnyGivenPoint(coeff, nx, ny, nz, u, v, w,
                                                     p, partialsOld);
    ReSamplingCubicSpline1OrderPartials(coeff, nx, ny, nz, dxyz, minExtent,
                                        partials, p);
    float sum = partials[1] * partials[1] + partials[2] * partials[2] +
                partials[3] * partials[3];
    float sum1 = partialsOld[1] * partialsOld[1] +
                 partialsOld[2] * partialsOld[2] +
                 partialsOld[3] * partialsOld[3];
    if (sum < 0.00001) {
      sum = 0.00001;
    }
    sum1 = sqrt(sum1);
    geo->m_TriVertNormals[i * 3] = -partials[3];     // /sum;
    geo->m_TriVertNormals[i * 3 + 1] = -partials[2]; // /sum;
    geo->m_TriVertNormals[i * 3 + 2] = -partials[1]; // /sum;
    H = -partials[1] * (partials[4] * partials[1] + partials[5] * partials[2] +
                        partials[6] * partials[3]) -
        partials[2] * (partials[5] * partials[1] + partials[7] * partials[2] +
                       partials[8] * partials[3]) -
        partials[3] * (partials[6] * partials[1] + partials[8] * partials[2] +
                       partials[9] * partials[3]);
    H = H / sum + partials[4] + partials[7] + partials[9];
    K = (partials[3] * partials[3] *
             (partials[4] * partials[7] - partials[5] * partials[5]) +
         partials[2] * partials[2] *
             (partials[4] * partials[9] - partials[6] * partials[6]) +
         partials[1] * partials[1] *
             (partials[7] * partials[9] - partials[8] * partials[8]) +
         2 * partials[1] * partials[2] *
             (partials[6] * partials[8] - partials[5] * partials[9]) +
         2 * partials[1] * partials[3] *
             (partials[5] * partials[8] - partials[6] * partials[7]) +
         2 * partials[2] * partials[3] *
             (partials[5] * partials[6] - partials[8] * partials[4])) /
        (sum * sum);
    float red, green, blue;
    geo->m_TriVertColorsTransparent[i * 3] = H;     // red/255.0;
    geo->m_TriVertColorsTransparent[i * 3 + 1] = K; // green/255.0;
    geo->m_TriVertColorsTransparent[i * 3 + 2] = 1;
    blue / 255.0;
  }
}

void HLevelSet::getLowerCorner(float *p, int &u, int &v, int &w) {
  float dx = dxyz[2];
  float dy = dxyz[1];
  float dz = dxyz[0]; // zq modified the order on 11/19/2008
  u = (int)((p[0] - minExtent[0]) / dz);
  v = (int)((p[1] - minExtent[1]) / dy);
  w = (int)((p[2] - minExtent[2]) / dx);
}

bool HLevelSet::getAtomListAndExtent(GroupOfAtoms *molecule,
                                     std::vector<PDBParser::Atom *> &atomList,
                                     float *minExtent, float *maxExtent) {
  if (!molecule || !minExtent || !maxExtent) {
    return false;
  }
  CollectionData *collectionData = 0;
  if (molecule->type == COLLECTION_TYPE) {
    collectionData = molecule->m_CollectionData;
  }
  GroupOfAtoms::RADIUS_TYPE radiusType = GroupOfAtoms::VDW_RADIUS;
  FlattenGOA(molecule, atomList, collectionData, 0, 0, 0, radiusType, ATOM_TYPE,
             false);
  minExtent[0] = minExtent[1] = minExtent[2] = 0.;
  maxExtent[0] = maxExtent[1] = maxExtent[2] = 0.;
  double probeRadius = 1.4;
  std::vector<float> Vertices;
  std::vector<float> Radius;
  float center[3];
  int size = atomList.size();
  for (int i = 0; i < size; i++) {
    PDBParser::Atom *atom = atomList[i];
    if (!atom) {
      continue;
    }
    center[0] = atom->m_Position[0];
    center[1] = atom->m_Position[1];
    center[2] = atom->m_Position[2];
    Vertices.push_back(center[0]);
    Vertices.push_back(center[1]);
    Vertices.push_back(center[2]);
    float radius = atom->getRadius();
    Radius.push_back(radius);
  }
  getBoundingBox(Vertices, Radius, minExtent, maxExtent, probeRadius * 2 + 1,
                 0);
  return true;
}

SimpleVolumeData *HLevelSet::getHigherOrderLevelSetSurface(
    GroupOfAtoms *molecule, unsigned int *dim, int nIterations, double weight) {
  // get the atom list from the molecule
  std::vector<PDBParser::Atom *> atomList;
  double probeRadius = 1.4;
  double blobbiness = 0;
  if (!getAtomListAndExtent(molecule, atomList, minExtent, maxExtent)) {
    return 0;
  }
  std::vector<float> Vertices;
  std::vector<float> Radius;
  float center[3];
  int size = atomList.size();
  for (int i = 0; i < size; i++) {
    PDBParser::Atom *atom = atomList[i];
    if (!atom) {
      continue;
    }
    center[0] = atom->m_Position[0];
    center[1] = atom->m_Position[1];
    center[2] = atom->m_Position[2];
    Vertices.push_back(center[0]);
    Vertices.push_back(center[1]);
    Vertices.push_back(center[2]);
    float radius = atom->getRadius();
    Radius.push_back(radius);
  }
  float *data = new float[dim[0] * dim[1] * dim[2]];
  if (!computeFunction(Vertices, Radius, data, dim, minExtent, maxExtent,
                       nIterations, weight)) {
    delete[] data;
    data = 0;
    return NULL;
  }
  // create volume data and return it.
  SimpleVolumeData *sData = new SimpleVolumeData(dim);
  sData->setDimensions(dim);
  sData->setNumberOfVariables(1);
  sData->setData(0, data);
  sData->setType(0, SimpleVolumeData::FLOAT);
  sData->setName(0, "HOrderLevelSet");
  sData->setMinExtent(minExtent);
  sData->setMaxExtent(maxExtent);
  return sData;
}

void HLevelSet::getBoundingBox(const std::vector<float> Vertices,
                               const std::vector<float> Radius, float min[3],
                               float max[3], int probeRadius,
                               double blobbiness) {
  int i, j;
  int size;
  float maxRad = 0.0;
  float *vertices, *radius;
  vertices = new float[Vertices.size()];
  radius = new float[Radius.size()];
  size = Radius.size();
  memcpy(vertices, &(Vertices[0]), Vertices.size() * sizeof(float));
  memcpy(radius, &(Radius[0]), Radius.size() * sizeof(float));
  if (size == 0) {
    for (i = 0; i < 3; i++) {
      min[i] = max[i] = 0.0;
    }
    return;
  }
  for (i = 0; i < 3; i++) {
    min[i] = vertices[i];
    max[i] = vertices[i];
  }
  float radius0 = radius[0];
  radius0 += probeRadius;
  float tempRad = radius0;
  if (blobbiness > 0.0001) {
    tempRad = radius0 * sqrt(1.0 + log(1e-2f) / blobbiness);
  }
  if (maxRad < tempRad) {
    maxRad = tempRad;
  }
  for (j = 1; j < size; j++) {
    for (i = 0; i < 3; i++) {
      if (vertices[3 * j + i] < min[i]) {
        min[i] = vertices[3 * j + i];
      }
      if (vertices[3 * j + i] > max[i]) {
        max[i] = vertices[3 * j + i];
      }
    }
    float radius0 = radius[j];
    radius0 += probeRadius;
    float tempRad = radius0;
    if (blobbiness > 0.0001) {
      tempRad = radius0 * sqrt(1.0 + log(1e-2f) / blobbiness);
    }
    if (maxRad < tempRad) {
      maxRad = tempRad;
    }
  }
  for (i = 0; i < 3; i++) {
    min[i] -= 4 * maxRad;
    max[i] += 4 * maxRad;
  }
}

void HLevelSet::EvaluateCubicSplineOrder1PartialsAtAnyGivenPoint(
    float *c, int nx, int ny, int nz, int u, int v, int w, float *p,
    float *partials)
// float *c, /* the spline coefficients */
// int nx, /* number of samples or coefficients in x direction */
// int ny, /* number of samples or coefficients in y direction */
// int nz, /* number of samples or coefficients in z direction */
// int  u, v, w,        /* the indices of the lower left corner of the bounding
// box */ float p[3],  /* the given point, be careful if any transformation is
// done    */ float *partial /* partial derivatives */
//              The minExtent[3] and maxExtent[3] in HLevelSet.h are the
//              bounding box of the molecule; dx, dy, dz are the grid length.
{
  float dx = dxyz[2];
  float dy = dxyz[1];
  float dz = dxyz[0];
  float c64[64], indx, indy, indz;
  float TensorF_p[64], TensorFx_p[64], TensorFy_p[64], TensorFz_p[64];
  int i;
  Take_64_Coefficients(c, nx, ny, nz, u, v, w, c64);
  ComputeTensorForGivenPoint(nx, ny, nz, u, v, w, p, TensorF_p, TensorFx_p,
                             TensorFy_p, TensorFz_p);
  for (i = 0; i < 4; i++) {
    partials[i] = 0.0;
  }
  for (i = 0; i < 64; i++) {
    partials[0] = partials[0] + c64[i] * TensorF_p[i];
    partials[1] = partials[1] + c64[i] * TensorFx_p[i];
    partials[2] = partials[2] + c64[i] * TensorFy_p[i];
    partials[3] = partials[3] + c64[i] * TensorFz_p[i];
  }
  indx = 1.0 / dx;
  indy = 1.0 / dy;
  indz = 1.0 / dz;
  partials[1] = partials[1] * indx;
  partials[2] = partials[2] * indy;
  partials[3] = partials[3] * indz;
}

void HLevelSet::ComputeTensorForGivenPoint(int nx, int ny, int nz, int u, int v,
                                           int w, float *p, float *TensorF_p,
                                           float *TensorFx_p, float *TensorFy_p,
                                           float *TensorFz_p)

// int nx, /* number of samples or coefficients in x direction */
// int ny, /* number of samples or coefficients in y direction */
// int nz, /* number of samples or coefficients in z direction */
{
  float f1[4], f2[4], f3[4];
  float f1x[4], f2x[4], f3x[4];
  float f1xx[4], f2xx[4], f3xx[4];
  int i, j, k, l;
  float x, y, z;
  float dx = dxyz[2];
  float dy = dxyz[1];
  float dz = dxyz[0];
  l = 0;
  if ((u > 0 && u < nx - 2) && (v > 0 && v < ny - 2) && (w > 0 && w < nz - 2)) {
    x = p[0] - (minExtent[0] + (u - 1) * dx);
    x = x / dx;
    f1[0] = (2.0 - x) * (2.0 - x) * (2.0 - x) / 6.0;
    f1x[0] = -0.5 * (2.0 - x) * (2.0 - x);
    x = p[0] - (minExtent[0] + u * dx);
    x = x / dx;
    f1[1] = 2.0 / 3.0 - x * x + 0.5 * x * x * x;
    f1x[1] = -2.0 * x + 1.5 * x * x;
    x = minExtent[0] + (u + 1) * dx - p[0];
    x = x / dx;
    f1[2] = 2.0 / 3.0 - x * x + 0.5 * x * x * x;
    f1x[2] = 2.0 * x - 1.5 * x * x;
    x = minExtent[0] + (u + 2) * dx - p[0];
    x = x / dx;
    f1[3] = (2.0 - x) * (2.0 - x) * (2.0 - x) / 6.0;
    f1x[3] = 0.5 * (2.0 - x) * (2.0 - x);
    y = p[1] - (minExtent[1] + (v - 1) * dy);
    y = y / dy;
    f2[0] = (2.0 - y) * (2.0 - y) * (2.0 - y) / 6.0;
    f2x[0] = -0.5 * (2.0 - y) * (2.0 - y);
    y = p[1] - (minExtent[1] + v * dy);
    y = y / dy;
    f2[1] = 2.0 / 3.0 - y * y + 0.5 * y * y * y;
    f2x[1] = -2.0 * y + 1.5 * y * y;
    y = minExtent[1] + (v + 1) * dy - p[1];
    y = y / dy;
    f2[2] = 2.0 / 3.0 - y * y + 0.5 * y * y * y;
    f2x[2] = 2.0 * y - 1.5 * y * y;
    y = minExtent[1] + (v + 2) * dy - p[1];
    y = y / dy;
    f2[3] = (2.0 - y) * (2.0 - y) * (2.0 - y) / 6.0;
    f2x[3] = 0.5 * (2.0 - y) * (2.0 - y);
    z = p[2] - (minExtent[2] + (w - 1) * dz);
    z = z / dz;
    f3[0] = (2.0 - z) * (2.0 - z) * (2.0 - z) / 6.0;
    f3x[0] = -0.5 * (2.0 - z) * (2.0 - z);
    z = p[2] - (minExtent[2] + w * dz);
    z = z / dz;
    f3[1] = 2.0 / 3.0 - z * z + 0.5 * z * z * z;
    f3x[1] = -2.0 * z + 1.5 * z * z;
    z = minExtent[2] + (w + 1) * dz - p[2];
    z = z / dz;
    f3[2] = 2.0 / 3.0 - z * z + 0.5 * z * z * z;
    f3x[2] = 2.0 * z - 1.5 * z * z;
    z = minExtent[2] + (w + 2) * dz - p[2];
    z = z / dz;
    f3[3] = (2.0 - z) * (2.0 - z) * (2.0 - z) / 6.0;
    f3x[3] = 0.5 * (2.0 - z) * (2.0 - z);
    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
        for (k = 0; k < 4; k++) {
          TensorF_p[l] = f1[i] * f2[j] * f3[k];
          TensorFx_p[l] = f1x[i] * f2[j] * f3[k];
          TensorFy_p[l] = f1[i] * f2x[j] * f3[k];
          TensorFz_p[l] = f1[i] * f2[j] * f3x[k];
          l = l + 1;
        }
      }
    }
    return;
  }
  x = minExtent[0] + (u - 1) * dx;
  if (x >= minExtent[0] && x <= maxExtent[0]) {
    x = p[0] - x;
    x = x / dx;
    f1[0] = (2.0 - x) * (2.0 - x) * (2.0 - x) / 6.0;
    f1x[0] = -0.5 * (2.0 - x) * (2.0 - x);
  } else {
    f1[0] = 0.0;
    f1x[0] = 0.0;
  }
  x = minExtent[0] + u * dx;
  if (x >= minExtent[0] && x <= maxExtent[0]) {
    x = p[0] - x;
    x = x / dx;
    f1[1] = 2.0 / 3.0 - x * x + 0.5 * x * x * x;
    f1x[1] = -2.0 * x + 1.5 * x * x;
  } else {
    f1[1] = 0.0;
    f1x[1] = 0.0;
  }
  x = minExtent[0] + (u + 1) * dx;
  if (x >= minExtent[0] && x <= maxExtent[0]) {
    x = x - p[0];
    x = x / dx;
    f1[2] = 2.0 / 3.0 - x * x + 0.5 * x * x * x;
    f1x[2] = 2.0 * x - 1.5 * x * x;
  } else {
    f1[2] = 0.0;
    f1x[2] = 0.0;
  }
  x = minExtent[0] + (u + 2) * dx;
  if (x >= minExtent[0] && x <= maxExtent[0]) {
    x = x - p[0];
    x = x / dx;
    f1[3] = (2.0 - x) * (2.0 - x) * (2.0 - x) / 6.0;
    f1x[3] = 0.5 * (2.0 - x) * (2.0 - x);
  } else {
    f1[3] = 0.0;
    f1x[3] = 0.0;
  }
  // y
  y = minExtent[1] + (v - 1) * dy;
  if (y >= minExtent[1] && y <= maxExtent[1]) {
    y = p[1] - y;
    y = y / dy;
    f2[0] = (2.0 - y) * (2.0 - y) * (2.0 - y) / 6.0;
    f2x[0] = -0.5 * (2.0 - y) * (2.0 - y);
  } else {
    f2[0] = 0.0;
    f2x[0] = 0.0;
  }
  y = minExtent[1] + v * dy;
  if (y >= minExtent[1] && y <= maxExtent[1]) {
    y = p[1] - y;
    y = y / dy;
    f2[1] = 2.0 / 3.0 - y * y + 0.5 * y * y * y;
    f2x[1] = -2.0 * y + 1.5 * y * y;
  } else {
    f2[1] = 0.0;
    f2x[1] = 0.0;
  }
  y = minExtent[1] + (v + 1) * dy;
  if (y >= minExtent[1] && y <= maxExtent[1]) {
    y = y - p[1];
    y = y / dy;
    f2[2] = 2.0 / 3.0 - y * y + 0.5 * y * y * y;
    f2x[2] = 2.0 * y - 1.5 * y * y;
  } else {
    f2[2] = 0.0;
    f2x[2] = 0.0;
  }
  y = minExtent[1] + (v + 2) * dy;
  if (y >= minExtent[1] && y <= maxExtent[1]) {
    y = y - p[1];
    y = y / dy;
    f2[3] = (2.0 - y) * (2.0 - y) * (2.0 - y) / 6.0;
    f2x[3] = 0.5 * (2.0 - y) * (2.0 - y);
  } else {
    f2[3] = 0.0;
    f2x[3] = 0.0;
  }
  // z
  z = minExtent[2] + (w - 1) * dz;
  if (z >= minExtent[2] && z <= maxExtent[2]) {
    z = p[2] - z;
    z = z / dz;
    f3[0] = (2.0 - z) * (2.0 - z) * (2.0 - z) / 6.0;
    f3x[0] = -0.5 * (2.0 - z) * (2.0 - z);
  } else {
    f3[0] = 0.0;
    f3x[0] = 0.0;
  }
  z = minExtent[2] + w * dz;
  if (z >= minExtent[2] && z <= maxExtent[2]) {
    z = p[2] - z;
    z = z / dz;
    f3[1] = 2.0 / 3.0 - z * z + 0.5 * z * z * z;
    f3x[1] = -2.0 * z + 1.5 * z * z;
  } else {
    f3[1] = 0.0;
    f3x[1] = 0.0;
  }
  z = minExtent[2] + (w + 1) * dz;
  if (z >= minExtent[2] && z <= maxExtent[2]) {
    z = z - p[2];
    z = z / dz;
    f3[2] = 2.0 / 3.0 - z * z + 0.5 * z * z * z;
    f3x[2] = 2.0 * z - 1.5 * z * z;
  } else {
    f3[2] = 0.0;
    f3x[2] = 0.0;
  }
  z = minExtent[2] + (w + 2) * dz;
  if (z >= minExtent[2] && z <= maxExtent[2]) {
    z = z - p[2];
    z = z / dz;
    f3[3] = (2.0 - z) * (2.0 - z) * (2.0 - z) / 6.0;
    f3x[3] = 0.5 * (2.0 - z) * (2.0 - z);
  } else {
    f3[3] = 0.0;
    f3x[3] = 0.0;
  }
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      for (k = 0; k < 4; k++) {
        TensorF_p[l] = f1[i] * f2[j] * f3[k];
        TensorFx_p[l] = f1x[i] * f2[j] * f3[k];
        TensorFy_p[l] = f1[i] * f2x[j] * f3[k];
        TensorFz_p[l] = f1[i] * f2[j] * f3x[k];
        l = l + 1;
      }
    }
  }
  return;
}

void ConvertToInterpolationCoefficients_1D(float *s, int DataLength,
                                           float EPSILON)
//              float   *s,             /* input samples --> output coefficients
//              */ int     DataLength,     /* number of samples or coefficients
//              */ float   EPSILON         /* admissible relative error */
{
  int i, n, ni, ni1, K;
  float sum, z1, w1, w2;
  n = DataLength + 1;
  z1 = sqrt(3.0) - 2.0;
  K = log(EPSILON) / log(fabs(z1));
  sum = 0.0;
  w2 = pow(z1, 2 * n);
  if (n < K) {
    for (i = 1; i < n; i++) {
      w1 = pow(z1, i);
      sum = sum + s[i - 1] * (w1 - w2 / w1);
    }
  } else {
    for (i = 1; i < n; i++) {
      sum = (sum + s[n - i - 1]) * z1;
    }
  }
  sum = -sum / (1.0 - w2);
  n = DataLength;
  s[0] = s[0] + z1 * sum;
  for (i = 1; i < n; i++) {
    s[i] = s[i] + z1 * s[i - 1];
  }
  // compute c^-
  s[n - 1] = -z1 * s[n - 1];
  for (i = 1; i < n; i++) {
    ni = n - i;
    ni1 = ni - 1;
    s[ni1] = z1 * (s[ni] - s[ni1]);
  }
  for (i = 0; i < n; i++) {
    s[i] = 6.0 * s[i];
  }
}

void ConvertToInterpolationCoefficients_2D(float *s, int nx, int ny,
                                           float EPSILON)
//              float   *s, /* input samples --> output coefficients */
//              int     nx,     /* number of samples or coefficients in x
//              direction */ int     ny,     /* number of samples or
//              coefficients in y direction */ float   EPSILON /* admissible
//              relative error */
{
  float *d, z1;
  int i, l;
  d = (float *)malloc(nx * sizeof(float));
  z1 = sqrt(3.0) - 2.0;
  // x-direction interpolation
  for (l = 0; l < ny; l++) {
    for (i = 0; i < nx; i++) {
      d[i] = s[i * ny + l];
    }
    ConvertToInterpolationCoefficients_1D(d, nx, CVC_DBL_EPSILON);
    for (i = 0; i < nx; i++) {
      s[i * ny + l] = d[i];
    }
  }
  // y-direction interpolation
  for (i = 0; i < nx; i++) {
    ConvertToInterpolationCoefficients_1D(s + i * ny, ny, CVC_DBL_EPSILON);
  }
  free(d);
}

void ConvertToInterpolationCoefficients_3D(float *s, int nx, int ny, int nz,
                                           float EPSILON)
// float *s, /* input samples --> output coefficients */
// int nx, /* number of samples or coefficients in x direction */
// int ny, /* number of samples or coefficients in y direction */
// int nz, /* number of samples or coefficients in z direction */
// float EPSILON/* admissible relative error */
{
  float *d, z1;
  int u, v, w, k, kk;
  d = (float *)malloc(nx * sizeof(float));
  k = ny * nz;
  z1 = sqrt(3.0) - 2.0;
  // x-direction interpolation
  for (v = 0; v < ny; v++) {
    for (w = 0; w < nz; w++) {
      kk = v * nz + w;
      for (u = 0; u < nx; u++) {
        d[u] = s[u * k + kk];
      }
      ConvertToInterpolationCoefficients_1D(d, nx, CVC_DBL_EPSILON);
      for (u = 0; u < nx; u++) {
        s[u * k + kk] = d[u];
      }
    }
  }
  for (u = 0; u < nx; u++) {
    ConvertToInterpolationCoefficients_2D(s + u * k, ny, nz, CVC_DBL_EPSILON);
  }
  free(d);
}

void ConvertToInterpolationCoefficients_3D_zq(float *s, int nx, int ny, int nz,
                                              float EPSILON) {
  float *d = new float[nx * ny * nz];
  float sum;
  int i, j, k;
  float lbda[7] = {0.0};
  lbda[1] = sqrt(3.0) - 2.0;
  lbda[2] = lbda[1] * lbda[1];
  lbda[3] = lbda[2] * lbda[1];
  lbda[4] = lbda[2] * lbda[2];
  lbda[5] = lbda[3] * lbda[2];
  lbda[6] = lbda[3] * lbda[3];
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        sum = 0.0;
        for (int ii = 1; ii <= 4; ii++) {
          if (k + ii < nx) {
            sum += s[(i * ny + j) * nz + k + ii] * lbda[ii];
          }
          if (k - ii >= 0) {
            sum += s[(i * ny + j) * nz + k - ii] * lbda[ii];
          }
        }
        d[(i * ny + j) * nz + k] = sqrt(3.0) * (s[(i * ny + j) * nz + k] + sum);
      }
    }
  }
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        sum = 0.0;
        for (int ii = 1; ii <= 4; ii++) {
          if (j + ii < ny) {
            sum += d[(i * ny + j + ii) * nz + k] * lbda[ii];
          }
          if (j - ii >= 0) {
            sum += d[(i * ny + j - ii) * nz + k] * lbda[ii];
          }
        }
        s[(i * ny + j) * nz + k] = sqrt(3.0) * (d[(i * ny + j) * nz + k] + sum);
      }
    }
  }
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        sum = 0.0;
        for (int ii = 1; ii <= 4; ii++) {
          if (i + ii < nz) {
            sum += s[((i + ii) * ny + j) * nz + k] * lbda[ii];
          }
          if (i - ii >= 0) {
            sum += s[((i - ii) * ny + j) * nz + k] * lbda[ii];
          }
        }
        d[(i * ny + j) * nz + k] = sqrt(3.0) * (s[(i * ny + j) * nz + k] + sum);
      }
    }
  }
  for (i = 0; i < nx * ny * nz; i++) {
    s[i] = d[i];
  }
  delete[] d;
}

float EvaluateCubicSplineAtGridPoint(float *c, int nx, int ny, int nz, int u,
                                     int v, int w)
// float *c, /* the spline coefficients */
// int nx, /* number of samples or coefficients in x direction */
// int ny, /* number of samples or coefficients in y direction */
// int nz, /* number of samples or coefficients in z direction */
{
  float result;
  if ((u > 0 && u < nx - 1) && (v > 0 && v < ny - 1) && (w > 0 && w < nz - 1)) {
    result =
        OneSix3 * (TakeACoefficient_Fast(c, nx, ny, nz, u - 1, v - 1, w - 1) +
                   TakeACoefficient_Fast(c, nx, ny, nz, u - 1, v - 1, w + 1) +
                   TakeACoefficient_Fast(c, nx, ny, nz, u - 1, v + 1, w + 1) +
                   TakeACoefficient_Fast(c, nx, ny, nz, u - 1, v + 1, w - 1) +
                   TakeACoefficient_Fast(c, nx, ny, nz, u + 1, v - 1, w - 1) +
                   TakeACoefficient_Fast(c, nx, ny, nz, u + 1, v - 1, w + 1) +
                   TakeACoefficient_Fast(c, nx, ny, nz, u + 1, v + 1, w + 1) +
                   TakeACoefficient_Fast(c, nx, ny, nz, u + 1, v + 1, w - 1));
    result = OneSix2 * (TakeACoefficient_Fast(c, nx, ny, nz, u, v - 1, w - 1) +
                        TakeACoefficient_Fast(c, nx, ny, nz, u, v - 1, w + 1) +
                        TakeACoefficient_Fast(c, nx, ny, nz, u, v + 1, w - 1) +
                        TakeACoefficient_Fast(c, nx, ny, nz, u, v + 1, w + 1) +
                        TakeACoefficient_Fast(c, nx, ny, nz, u - 1, v, w - 1) +
                        TakeACoefficient_Fast(c, nx, ny, nz, u - 1, v, w + 1) +
                        TakeACoefficient_Fast(c, nx, ny, nz, u + 1, v, w - 1) +
                        TakeACoefficient_Fast(c, nx, ny, nz, u + 1, v, w + 1) +
                        TakeACoefficient_Fast(c, nx, ny, nz, u - 1, v - 1, w) +
                        TakeACoefficient_Fast(c, nx, ny, nz, u - 1, v + 1, w) +
                        TakeACoefficient_Fast(c, nx, ny, nz, u + 1, v - 1, w) +
                        TakeACoefficient_Fast(c, nx, ny, nz, u + 1, v + 1, w)) +
             result;
    result = OneSix1 * (TakeACoefficient_Fast(c, nx, ny, nz, u, v, w - 1) +
                        TakeACoefficient_Fast(c, nx, ny, nz, u, v, w + 1) +
                        TakeACoefficient_Fast(c, nx, ny, nz, u, v - 1, w) +
                        TakeACoefficient_Fast(c, nx, ny, nz, u, v + 1, w) +
                        TakeACoefficient_Fast(c, nx, ny, nz, u - 1, v, w) +
                        TakeACoefficient_Fast(c, nx, ny, nz, u + 1, v, w)) +
             result;
    result = OneSix0 * TakeACoefficient_Fast(c, nx, ny, nz, u, v, w) + result;
    return (result);
  }
  result =
      OneSix3 * (TakeACoefficient_Slow(c, nx, ny, nz, u - 1, v - 1, w - 1) +
                 TakeACoefficient_Slow(c, nx, ny, nz, u - 1, v - 1, w + 1) +
                 TakeACoefficient_Slow(c, nx, ny, nz, u - 1, v + 1, w + 1) +
                 TakeACoefficient_Slow(c, nx, ny, nz, u - 1, v + 1, w - 1) +
                 TakeACoefficient_Slow(c, nx, ny, nz, u + 1, v - 1, w - 1) +
                 TakeACoefficient_Slow(c, nx, ny, nz, u + 1, v - 1, w + 1) +
                 TakeACoefficient_Slow(c, nx, ny, nz, u + 1, v + 1, w + 1) +
                 TakeACoefficient_Slow(c, nx, ny, nz, u + 1, v + 1, w - 1));
  result = OneSix2 * (TakeACoefficient_Slow(c, nx, ny, nz, u, v - 1, w - 1) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u, v - 1, w + 1) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u, v + 1, w - 1) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u, v + 1, w + 1) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u - 1, v, w - 1) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u - 1, v, w + 1) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u + 1, v, w - 1) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u + 1, v, w + 1) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u - 1, v - 1, w) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u - 1, v + 1, w) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u + 1, v - 1, w) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u + 1, v + 1, w)) +
           result;
  result = OneSix1 * (TakeACoefficient_Slow(c, nx, ny, nz, u, v, w - 1) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u, v, w + 1) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u, v - 1, w) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u, v + 1, w) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u - 1, v, w) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u + 1, v, w)) +
           result;
  result = OneSix0 * TakeACoefficient_Slow(c, nx, ny, nz, u, v, w) + result;
  return (result);
}

float TakeACoefficient_Fast(float *c, int nx, int ny, int nz, int u, int v,
                            int w) {
  return (c[(u * ny + v) * nz + w]);
}

float TakeACoefficient_Slow(float *c, int nx, int ny, int nz, int u, int v,
                            int w) {
  float result;
  result = 0.0;
  if ((u >= 0 && u < nx) && (v >= 0 && v < ny) && (w >= 0 && w < nz)) {
    result = c[(u * ny + v) * nz + w];
  }
  return (result);
}

void Tensor_333(float *xx, float *yy, float *zz, float *result) {
  int i, j, k, l;
  l = 0;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        result[l] = xx[i] * yy[j] * zz[k];
        l = l + 1;
      }
    }
  }
}

void Tensor_444(float *xx, float *yy, float *zz, float *result) {
  int i, j, k, l;
  l = 0;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      for (k = 0; k < 4; k++) {
        result[l] = xx[i] * yy[j] * zz[k];
        l = l + 1;
      }
    }
  }
}

void Take_27_Coefficients(float *c, int nx, int ny, int nz, int u, int v, int w,
                          float *c27)
// float *c, /* the spline coefficients */
// int nx, /* number of samples or coefficients in x direction */
// int ny, /* number of samples or coefficients in y direction */
// int nz, /* number of samples or coefficients in z direction */
// float *c27 /* 27 coefficients */
{
  if ((u > 0 && u < nx - 1) && (v > 0 && v < ny - 1) && (w > 0 && w < nz - 1)) {
    c27[0] = TakeACoefficient_Fast(c, nx, ny, nz, u - 1, v - 1, w - 1);
    c27[1] = TakeACoefficient_Fast(c, nx, ny, nz, u - 1, v - 1, w);
    c27[2] = TakeACoefficient_Fast(c, nx, ny, nz, u - 1, v - 1, w + 1);
    c27[3] = TakeACoefficient_Fast(c, nx, ny, nz, u - 1, v, w - 1);
    c27[4] = TakeACoefficient_Fast(c, nx, ny, nz, u - 1, v, w);
    c27[5] = TakeACoefficient_Fast(c, nx, ny, nz, u - 1, v, w + 1);
    c27[6] = TakeACoefficient_Fast(c, nx, ny, nz, u - 1, v + 1, w - 1);
    c27[7] = TakeACoefficient_Fast(c, nx, ny, nz, u - 1, v + 1, w);
    c27[8] = TakeACoefficient_Fast(c, nx, ny, nz, u - 1, v + 1, w + 1);
    c27[9] = TakeACoefficient_Fast(c, nx, ny, nz, u, v - 1, w - 1);
    c27[10] = TakeACoefficient_Fast(c, nx, ny, nz, u, v - 1, w);
    c27[11] = TakeACoefficient_Fast(c, nx, ny, nz, u, v - 1, w + 1);
    c27[12] = TakeACoefficient_Fast(c, nx, ny, nz, u, v, w - 1);
    c27[13] = TakeACoefficient_Fast(c, nx, ny, nz, u, v, w);
    c27[14] = TakeACoefficient_Fast(c, nx, ny, nz, u, v, w + 1);
    c27[15] = TakeACoefficient_Fast(c, nx, ny, nz, u, v + 1, w - 1);
    c27[16] = TakeACoefficient_Fast(c, nx, ny, nz, u, v + 1, w);
    c27[17] = TakeACoefficient_Fast(c, nx, ny, nz, u, v + 1, w + 1);
    c27[18] = TakeACoefficient_Fast(c, nx, ny, nz, u + 1, v - 1, w - 1);
    c27[19] = TakeACoefficient_Fast(c, nx, ny, nz, u + 1, v - 1, w);
    c27[20] = TakeACoefficient_Fast(c, nx, ny, nz, u + 1, v - 1, w + 1);
    c27[21] = TakeACoefficient_Fast(c, nx, ny, nz, u + 1, v, w - 1);
    c27[22] = TakeACoefficient_Fast(c, nx, ny, nz, u + 1, v, w);
    c27[23] = TakeACoefficient_Fast(c, nx, ny, nz, u + 1, v, w + 1);
    c27[24] = TakeACoefficient_Fast(c, nx, ny, nz, u + 1, v + 1, w - 1);
    c27[25] = TakeACoefficient_Fast(c, nx, ny, nz, u + 1, v + 1, w);
    c27[26] = TakeACoefficient_Fast(c, nx, ny, nz, u + 1, v + 1, w + 1);
    return;
  }
  c27[0] = TakeACoefficient_Slow(c, nx, ny, nz, u - 1, v - 1, w - 1);
  c27[1] = TakeACoefficient_Slow(c, nx, ny, nz, u - 1, v - 1, w);
  c27[2] = TakeACoefficient_Slow(c, nx, ny, nz, u - 1, v - 1, w + 1);
  c27[3] = TakeACoefficient_Slow(c, nx, ny, nz, u - 1, v, w - 1);
  c27[4] = TakeACoefficient_Slow(c, nx, ny, nz, u - 1, v, w);
  c27[5] = TakeACoefficient_Slow(c, nx, ny, nz, u - 1, v, w + 1);
  c27[6] = TakeACoefficient_Slow(c, nx, ny, nz, u - 1, v + 1, w - 1);
  c27[7] = TakeACoefficient_Slow(c, nx, ny, nz, u - 1, v + 1, w);
  c27[8] = TakeACoefficient_Slow(c, nx, ny, nz, u - 1, v + 1, w + 1);
  c27[9] = TakeACoefficient_Slow(c, nx, ny, nz, u, v - 1, w - 1);
  c27[10] = TakeACoefficient_Slow(c, nx, ny, nz, u, v - 1, w);
  c27[11] = TakeACoefficient_Slow(c, nx, ny, nz, u, v - 1, w + 1);
  c27[12] = TakeACoefficient_Slow(c, nx, ny, nz, u, v, w - 1);
  c27[13] = TakeACoefficient_Slow(c, nx, ny, nz, u, v, w);
  c27[14] = TakeACoefficient_Slow(c, nx, ny, nz, u, v, w + 1);
  c27[15] = TakeACoefficient_Slow(c, nx, ny, nz, u, v + 1, w - 1);
  c27[16] = TakeACoefficient_Slow(c, nx, ny, nz, u, v + 1, w);
  c27[17] = TakeACoefficient_Slow(c, nx, ny, nz, u, v + 1, w + 1);
  c27[18] = TakeACoefficient_Slow(c, nx, ny, nz, u + 1, v - 1, w - 1);
  c27[19] = TakeACoefficient_Slow(c, nx, ny, nz, u + 1, v - 1, w);
  c27[20] = TakeACoefficient_Slow(c, nx, ny, nz, u + 1, v - 1, w + 1);
  c27[21] = TakeACoefficient_Slow(c, nx, ny, nz, u + 1, v, w - 1);
  c27[22] = TakeACoefficient_Slow(c, nx, ny, nz, u + 1, v, w);
  c27[23] = TakeACoefficient_Slow(c, nx, ny, nz, u + 1, v, w + 1);
  c27[24] = TakeACoefficient_Slow(c, nx, ny, nz, u + 1, v + 1, w - 1);
  c27[25] = TakeACoefficient_Slow(c, nx, ny, nz, u + 1, v + 1, w);
  c27[26] = TakeACoefficient_Slow(c, nx, ny, nz, u + 1, v + 1, w + 1);
}

void Take_64_Coefficients(float *c, int nx, int ny, int nz, int u, int v, int w,
                          float *c64)
// float *c, /* the spline coefficients */
// int nx, /* number of samples or coefficients in x direction */
// int ny, /* number of samples or coefficients in y direction */
// int nz, /* number of samples or coefficients in z direction */
// float *c64 /* 64 coefficients */
{
  int i, j, k, l;
  if ((u > 0 && u < nx - 2) && (v > 0 && v < ny - 2) && (w > 0 && w < nz - 2)) {
    l = 0;
    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
        for (k = 0; k < 4; k++) {
          c64[l] = TakeACoefficient_Fast(c, nx, ny, nz, u - 1 + i, v - 1 + j,
                                         w - 1 + k);
          l = l + 1;
        }
      }
    }
    return;
  }
  l = 0;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      for (k = 0; k < 4; k++) {
        c64[l] = TakeACoefficient_Slow(c, nx, ny, nz, u - 1 + i, v - 1 + j,
                                       w - 1 + k);
        l = l + 1;
      }
    }
  }
}

void EvaluateCubicSplineOrder2PartialsAtGridPoint(float *c, float dx, float dy,
                                                  float dz, int nx, int ny,
                                                  int nz, int u, int v, int w,
                                                  float *partials)
// float *c, /* the spline coefficients */
// float dx /* spacing in x direction */
// float dy /* spacing in y direction */
// float dz /* spacing in z direction */
// int nx, /* number of samples or coefficients in x direction */
// int ny, /* number of samples or coefficients in y direction */
// int nz, /* number of samples or coefficients in z direction */
// float *partial /* partial derivatives */
{
  float c27[27], indx, indy, indz;
  Take_27_Coefficients(c, nx, ny, nz, u, v, w, c27);
  memset(partials, 0, sizeof(float) * 10);
  /*for(int i = 0; i < 10; i++)
  {
          partials[i] = 0.0;
  }*/
  for (int i = 0; i < 27; i++) {
    partials[0] = partials[0] + c27[i] * TensorF[i];
    partials[1] = partials[1] + c27[i] * TensorFx[i];
    partials[2] = partials[2] + c27[i] * TensorFy[i];
    partials[3] = partials[3] + c27[i] * TensorFz[i];
    partials[4] = partials[4] + c27[i] * TensorFxx[i];
    partials[5] = partials[5] + c27[i] * TensorFxy[i];
    partials[6] = partials[6] + c27[i] * TensorFxz[i];
    partials[7] = partials[7] + c27[i] * TensorFyy[i];
    partials[8] = partials[8] + c27[i] * TensorFyz[i];
    partials[9] = partials[9] + c27[i] * TensorFzz[i];
  }
  indx = 1.0 / dx;
  indy = 1.0 / dy;
  indz = 1.0 / dz;
  partials[1] = partials[1] * indx;
  partials[2] = partials[2] * indy;
  partials[3] = partials[3] * indz;
  partials[4] = partials[4] * indx * indx;
  partials[5] = partials[5] * indx * indy;
  partials[6] = partials[6] * indx * indz;
  partials[7] = partials[7] * indy * indy;
  partials[8] = partials[8] * indy * indz;
  partials[9] = partials[9] * indz * indz;
}

void EvaluateCubicSplineOrder1PartialsAtGridPoint(float *c, float dx, float dy,
                                                  float dz, int nx, int ny,
                                                  int nz, int u, int v, int w,
                                                  float *partials)
// float *c, /* the spline coefficients */
// float dx /* spacing in x direction */
// float dy /* spacing in y direction */
// float dz /* spacing in z direction */
// int nx, /* number of samples or coefficients in x direction */
// int ny, /* number of samples or coefficients in y direction */
// int nz, /* number of samples or coefficients in z direction */
// float *partial /* partial derivatives */
{
  float c27[27], indx, indy, indz;
  int i;
  Take_27_Coefficients(c, nx, ny, nz, u, v, w, c27);
  for (i = 0; i < 4; i++) {
    partials[i] = 0.0;
  }
  for (i = 0; i < 27; i++) {
    partials[0] = partials[0] + c27[i] * TensorF[i];
    partials[1] = partials[1] + c27[i] * TensorFx[i];
    partials[2] = partials[2] + c27[i] * TensorFy[i];
    partials[3] = partials[3] + c27[i] * TensorFz[i];
  }
  indx = 1.0 / dx;
  indy = 1.0 / dy;
  indz = 1.0 / dz;
  partials[1] = partials[1] * indx;
  partials[2] = partials[2] * indy;
  partials[3] = partials[3] * indz;
}

void InitialData_Sphere(float *v, float minx, float maxx, float miny,
                        float maxy, float minz, float maxz, int nx, int ny,
                        int nz)
// float *v, /* the volume data of initial signed distance function */
// float minx /* left end-point of x diection */
// float maxx /* right end-point of x diection */
// float miny /* left end-point of y diection */
// float maxy /* right end-point of y diection */
// float minz /* left end-point of z diection */
// float maxz /* right end-point of z diection */
// int nx, /* number of points in x direction */
// int ny, /* number of points in y direction */
// int nz, /* number of points in z direction */
{
  int i, j, k, l;
  float dx, dy, dz, x, y, z;
  float cx, cy, cz, radius, dis;
  dx = (maxx - minx) / (nx - 1.0);
  dy = (maxy - miny) / (ny - 1.0);
  dz = (maxz - minz) / (nz - 1.0);
  // determine the center and radius
  cx = (maxx + minx) / 2.0;
  cy = (maxy + miny) / 2.0;
  cz = (maxz + minz) / 2.0;
  radius = (maxx - minx) / 2.0;
  if (radius > (maxy - miny) / 2.0) {
    radius = (maxy - miny) / 2.0;
  }
  if (radius > (maxz - minz) / 2.0) {
    radius = (maxz - minz) / 2.0;
  }
  radius = radius - 3 * (dx + dy + dz);
  // printf("radius = %f\n", radius);
  //  compute function value at grid points
  l = 0;
  for (i = 0; i < nx; i++) {
    x = minx + i * dx;
    for (j = 0; j < ny; j++) {
      y = miny + j * dy;
      for (k = 0; k < nz; k++) {
        z = minz + k * dz;
        dis = sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy) +
                   (z - cz) * (z - cz));
        dis = dis - radius;
        // cut
        if (dis < -Height) {
          dis = -Height1;
        }
        if (dis > Height) {
          dis = Height1;
        }
        v[l] = dis;
        l = l + 1;
      }
    }
  }
}

void InitialData_Box(float *v, float minx, float maxx, float miny, float maxy,
                     float minz, float maxz, int nx, int ny, int nz)
// float *v, /* the volume data of initial signed distance function */
// float minx /* left end-point of x diection */
// float maxx /* right end-point of x diection */
// float miny /* left end-point of y diection */
// float maxy /* right end-point of y diection */
// float minz /* left end-point of z diection */
// float maxz /* right end-point of z diection */
// int nx, /* number of points in x direction */
// int ny, /* number of points in y direction */
// int nz, /* number of points in z direction */
{
  int i, j, k, l;
  float dx, dy, dz, x, y, z;
  float cx, cy, cz, radius, dis;
  dx = (maxx - minx) / (nx - 1.0);
  dy = (maxy - miny) / (ny - 1.0);
  dz = (maxz - minz) / (nz - 1.0);
  // determine the center and radius
  cx = (maxx + minx) / 2.0;
  cy = (maxy + miny) / 2.0;
  cz = (maxz + minz) / 2.0;
  radius = (maxx - minx) / 2.0;
  if (radius > (maxy - miny) / 2.0) {
    radius = (maxy - miny) / 2.0;
  }
  if (radius > (maxz - minz) / 2.0) {
    radius = (maxz - minz) / 2.0;
  }
  radius = radius - 3 * (dx + dy + dz);
  // radius = radius - dx;
  // printf("radius = %f\n", radius);
  //  compute function value at grid points
  l = 0;
  for (i = 0; i < nx; i++) {
    x = fabs(minx + i * dx);
    for (j = 0; j < ny; j++) {
      y = fabs(miny + j * dy);
      for (k = 0; k < nz; k++) {
        z = fabs(minz + k * dz);
        if (x >= radius && y >= radius && z >= radius) {
          dis = sqrt((x - radius) * (x - radius) + (y - radius) * (y - radius) +
                     (z - radius) * (z - radius));
        } else {
          dis = 0.0;
          if (x >= radius) {
            dis = (x - radius) * (x - radius);
          }
          if (y >= radius) {
            dis = (y - radius) * (y - radius) + dis;
          }
          if (z >= radius) {
            dis = (z - radius) * (z - radius) + dis;
          }
          dis = sqrt(dis);
          if (x >= radius && y <= radius && z <= radius) {
            dis = x - radius;
          }
          if (x <= radius && y >= radius && z <= radius) {
            dis = y - radius;
          }
          if (x <= radius && y <= radius && z >= radius) {
            dis = z - radius;
          }
          if (x <= radius && y <= radius && z <= radius) {
            dis = radius - x;
            if (dis > radius - y) {
              dis = radius - y;
            }
            if (dis > radius - z) {
              dis = radius - z;
            }
            dis = -dis;
          }
        }
        // cut
        if (dis < -Height) {
          dis = -Height1;
        }
        if (dis > Height) {
          dis = Height1;
        }
        v[l] = dis;
        l = l + 1;
      }
    }
  }
}

void MeanCurvatureFlow(float *fun, float *fun_bak, float *coeff, float dt,
                       float dx, float dy, float dz, int nx, int ny, int nz,
                       int t)
// float *fun /* the volume data of initial signed distance function */
// float *coeff /* spline coefficients */
// float dt /* time step size */
// float dx /* x direction spacing */
// float dy /* y direction spacing */
// float dz /* z direction spacing */
// int nx, /* number of points in x direction */
// int ny, /* number of points in y direction */
// int nz, /* number of points in z direction */
{
  int i, j, k, index;
  float partials[10], grad1, grad2, H, w, ww, ww1, ww2;
  float x, y, z, r, fvalue, maxf, maxt;
  float radius;
  radius = 5.0 - 6 * dx;
  ww1 = Height - 3 * Height0;
  ww2 = (Height - Height0) * (Height - Height0) * (Height - Height0);
  ww2 = 1.0 / ww2;
  maxf = -1000.0;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        index = (i * ny + j) * nz + k;
        fvalue = fun[index];
        fun_bak[index] = 0.0;
        if (fabs(fvalue) < Height) {
          EvaluateCubicSplineOrder2PartialsAtGridPoint(
              coeff, dx, dy, dz, nx, ny, nz, i, j, k, partials);
          grad2 = partials[1] * partials[1] + partials[2] * partials[2] +
                  partials[3] * partials[3];
          if (grad2 < 0.0001) {
            grad2 = 0.0001;
          }
          grad1 = sqrt(grad2);
          H = -partials[1] *
                  (partials[4] * partials[1] + partials[5] * partials[2] +
                   partials[6] * partials[3]) -
              partials[2] *
                  (partials[5] * partials[1] + partials[7] * partials[2] +
                   partials[8] * partials[3]) -
              partials[3] *
                  (partials[6] * partials[1] + partials[8] * partials[2] +
                   partials[9] * partials[3]);
          H = H / grad2 + partials[4] + partials[7] + partials[9];
          x = -5.0 + i * dx;
          y = -5.0 + j * dy;
          z = -5.0 + k * dz;
          r = sqrt(x * x + y * y + z * z);
          w = H * grad1;
          ww = fabs(w);
          if (fabs(fvalue) < Height0) {
            if (ww > maxf) {
              maxf = ww;
            }
          }
          fun_bak[index] = w;
        }
      }
    }
  }
  maxt = dx / maxf;
  if (maxt > dt) {
    maxt = dt;
  }
  // printf("maxf = %f, dx = %f, maxt = %f\n", maxf, dx, maxt);
  TotalTime = TotalTime + maxt;
  for (i = 0; i < nx * ny * nz; i++) {
    fun[i] = fun[i] + maxt * fun_bak[i];
    // cut the function by Height
    if (fun[i] < -Height) {
      fun[i] = -Height1;
    }
    if (fun[i] > Height) {
      fun[i] = Height1;
    }
  }
}

void Constraint_MeanCurvatureFlow(float *func_h, float *fun, float *fun_bak,
                                  float *coeff, float weight, float dt,
                                  float dx, float dy, float dz, int nx, int ny,
                                  int nz, int t)
// float *fun /* the volume data of initial signed distance function */
// float *coeff /* spline coefficients */
// float dt /* time step size */
// float dx /* x direction spacing */
// float dy /* y direction spacing */
// float dz /* z direction spacing */
// int nx, /* number of points in x direction */
// int ny, /* number of points in y direction */
// int nz, /* number of points in z direction */
{
  int i, j, k, index;
  float partials[10], partials_h[10], grad1, grad2, H, w, ww, ww1, ww2;
  float x, y, z, r, fvalue, maxf, maxt;
  float radius;
  radius = 5.0 - 6 * dx;
  ww1 = Height - 3 * Height0;
  ww2 = (Height - Height0) * (Height - Height0) * (Height - Height0);
  ww2 = 1.0 / ww2;
  maxf = -1000.0;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        index = (i * ny + j) * nz + k;
        fvalue = fun[index];
        fun_bak[index] = 0.0;
        if (fabs(fvalue) < Height) {
          EvaluateCubicSplineOrder2PartialsAtGridPoint(
              coeff, dx, dy, dz, nx, ny, nz, i, j, k, partials);
          grad2 = partials[1] * partials[1] + partials[2] * partials[2] +
                  partials[3] * partials[3];
          if (grad2 < 0.0001) {
            grad2 = 0.0001;
          }
          grad1 = sqrt(grad2);
          H = -partials[1] *
                  (partials[4] * partials[1] + partials[5] * partials[2] +
                   partials[6] * partials[3]) -
              partials[2] *
                  (partials[5] * partials[1] + partials[7] * partials[2] +
                   partials[8] * partials[3]) -
              partials[3] *
                  (partials[6] * partials[1] + partials[8] * partials[2] +
                   partials[9] * partials[3]);
          H = H / grad2 + partials[4] + partials[7] + partials[9];
          x = -5.0 + i * dx;
          y = -5.0 + j * dy;
          z = -5.0 + k * dz;
          r = sqrt(x * x + y * y + z * z);
          EvaluateCubicSplineOrder2PartialsAtGridPoint(
              func_h, dx, dy, dz, nx, ny, nz, i, j, k, partials_h);
          w = (partials_h[0] * partials_h[0] + weight) * H * grad1 +
              4 * partials_h[0] *
                  (partials_h[1] * partials[1] + partials_h[2] * partials[2] +
                   partials_h[3] * partials[3]);
          ww = fabs(w);
          if (fabs(fvalue) < Height0) {
            if (ww > maxf) {
              maxf = ww;
            }
          }
          fun_bak[index] = w;
        }
      }
    }
  }
  maxt = dx / maxf;
  if (maxt > dt) {
    maxt = dt;
  }
  // printf("maxf = %f, dx = %f, maxt = %f\n", maxf, dx, maxt);
  TotalTime = TotalTime + maxt;
  // update fun
  for (i = 0; i < nx * ny * nz; i++) {
    fun[i] = fun[i] + maxt * fun_bak[i];
    // cut the function by Height
    if (fun[i] < -Height) {
      fun[i] = -Height1;
    }
    if (fun[i] > Height) {
      fun[i] = Height1;
    }
  }
}

void ComputeTensorXYZ() {
  float f[3], fx[3], fxx[3];
  f[0] = 1.0 / 6.0;
  f[1] = 2.0 / 3.0;
  f[2] = f[0];
  fx[0] = -0.5;
  fx[1] = 0.0;
  fx[2] = 0.5;
  fxx[0] = 1.0;
  fxx[1] = -2.0;
  fxx[2] = 1.0;
  Tensor_333(f, f, f, TensorF);
  Tensor_333(fx, f, f, TensorFx);
  Tensor_333(f, fx, f, TensorFy);
  Tensor_333(f, f, fx, TensorFz);
  Tensor_333(fxx, f, f, TensorFxx);
  Tensor_333(fx, fx, f, TensorFxy);
  Tensor_333(fx, f, fx, TensorFxz);
  Tensor_333(f, fxx, f, TensorFyy);
  Tensor_333(f, fx, fx, TensorFyz);
  Tensor_333(f, f, fxx, TensorFzz);
}

void ReInitilazation(float *fun, float *coeff, float dx, float dy, float dz,
                     int nx, int ny, int nz, float *minExtent)
// float *fun /* the volume data of initial signed distance function */
// float *coeff /* spline coefficients */
// float dx /* x direction spacing */
// float dy /* y direction spacing */
// float dz /* z direction spacing */
// int nx, /* number of points in x direction */
// int ny, /* number of points in y direction */
// int nz, /* number of points in z direction */
{
  int i, j, k, index, t;
  float partials[4], grad1, grad2;
  float x, y, z, r;
  float dt, sd;
  float ddxyz[3];
  ddxyz[0] = dx;
  ddxyz[1] = dy;
  ddxyz[2] = dz;
  float p[3];
  dt = 0.01;
  FILE *fp = fopen("output.txt", "w");
  FILE *fp1 = fopen("outputzq.txt", "w");
  for (t = 0; t < 1; t++) {
    for (i = 0; i < nx; i++) {
      for (j = 0; j < ny; j++) {
        for (k = 0; k < nz; k++) {
          EvaluateCubicSplineOrder1PartialsAtGridPoint(
              coeff, dx, dy, dz, nx, ny, nz, i, j, k, partials);
          grad2 = partials[1] * partials[1] + partials[2] * partials[2] +
                  partials[3] * partials[3];
          grad1 = sqrt(grad2);
          sd = partials[0] / sqrt(partials[0] * partials[0] + grad2 * dx * dx);
          fprintf(fp, "i , j, k = %d, %d, %d Dxyz = %f, %f, %f\n", i, j, k,
                  partials[1] / grad1, partials[2] / grad1,
                  partials[3] / grad1);
          p[0] = minExtent[0] + i * dx;
          p[1] = minExtent[1] + j * dy;
          p[2] = minExtent[2] + k * dz;
          ReSamplingCubicSpline1OrderPartials(coeff, nx, ny, nz, ddxyz,
                                              minExtent, partials, p);
          grad2 = partials[1] * partials[1] + partials[2] * partials[2] +
                  partials[0] * partials[0];
          grad1 = sqrt(grad2);
          fprintf(fp1, "i , j, k = %d, %d, %d Dxyz = %f, %f, %f\n", i, j, k,
                  partials[0] / grad1, partials[1] / grad1,
                  partials[2] / grad1);
          x = -5.0 + i * dx;
          y = -5.0 + j * dy;
          z = -5.0 + k * dz;
          r = sqrt(x * x + y * y + z * z);
          /*if(i == 44 && j== 31 && k == 26)
          {
                  printf("t = %d, i , j, k = %d, %d, %d, grad1 = %f S(d) =
          %f\n", t, i, j, k, grad1, sd);
          }*/
          index = (i * ny + j) * nz + k;
          fun[index] = fun[index] - dt * sd * (grad1 - 1.0);
        }
      }
    }
    for (i = 0; i < nx * ny * nz; i++) {
      coeff[i] = fun[i];
    }
    ConvertToInterpolationCoefficients_3D(coeff, nx, ny, nz, CVC_DBL_EPSILON);
  }
  fclose(fp);
  fclose(fp1);
}

void ReInitilazation_Upwind_Eno_Engquist(float *fun, float *fun_bak,
                                         float *coeff, float dx, float dy,
                                         float dz, int nx, int ny, int nz)
// float *fun /* the volume data of initial signed distance function */
// float *coeff /* spline coefficients */
// float dx /* x direction spacing */
// float dy /* y direction spacing */
// float dz /* z direction spacing */
// int nx, /* number of points in x direction */
// int ny, /* number of points in y direction */
// int nz, /* number of points in z direction */
{
  int i, j, k, index, t;
  float fvalue, grad1, grad2, ww, ww1;
  float x, y, z, r, fx[7], fy[7], fz[7];
  float dt, sd, sdp, sdm;
  float a, b, c, d, e, f;
  float ap, bp, cp, dp, ep, fp;
  float am, bm, cm, dm, em, fm;
  dt = 0.5 * dx;
  dt = 0.05;
  for (t = 0; t < 100; t++) {
    for (i = 0; i < nx; i++) {
      for (j = 0; j < ny; j++) {
        for (k = 0; k < nz; k++) {
          Get_Seven_Function_Values_X(fun, nx, ny, nz, i, j, k, fx);
          Quadr_Eno_Interpolation(fx, &b, &a, dx);
          ap = a;
          if (a < 0.0) {
            ap = 0.0;
          }
          bp = b;
          if (b < 0.0) {
            bp = 0.0;
          }
          am = a;
          if (a > 0.0) {
            am = 0.0;
          }
          bm = b;
          if (b > 0.0) {
            bm = 0.0;
          }
          Get_Seven_Function_Values_Y(fun, nx, ny, nz, i, j, k, fy);
          Quadr_Eno_Interpolation(fy, &d, &c, dy);
          cp = c;
          if (c < 0.0) {
            cp = 0.0;
          }
          dp = d;
          if (d < 0.0) {
            dp = 0.0;
          }
          cm = c;
          if (c > 0.0) {
            cm = 0.0;
          }
          dm = d;
          if (d > 0.0) {
            dm = 0.0;
          }
          Get_Seven_Function_Values_Z(fun, nx, ny, nz, i, j, k, fz);
          Quadr_Eno_Interpolation(fz, &f, &e, dz);
          ep = e;
          if (e < 0.0) {
            ep = 0.0;
          }
          fp = f;
          if (f < 0.0) {
            fp = 0.0;
          }
          em = e;
          if (e > 0.0) {
            em = 0.0;
          }
          fm = f;
          if (f > 0.0) {
            fm = 0.0;
          }
          grad2 = 0.5 * (a * a + c * c + e * e + b * b + d * d + f * f);
          grad1 = sqrt(grad2);
          index = (i * ny + j) * nz + k;
          fvalue = fun[index];
          sd = fvalue / sqrt(fvalue * fvalue + grad2 * dx * dx);
          sdp = sd;
          if (sd < 0.0) {
            sdp = 0.0;
          }
          sdm = sd;
          if (sd > 0.0) {
            sdm = 0.0;
          }
          x = -5.0 + i * dx;
          y = -5.0 + j * dy;
          z = -5.0 + k * dz;
          r = sqrt(x * x + y * y + z * z);
          if (sdp > 0.0) {
            ap = ap * ap;
            cp = cp * cp;
            ep = ep * ep;
            bm = bm * bm;
            dm = dm * dm;
            fm = fm * fm;
            // Engquist-Osher scheme
            ww = sqrt(ap + cp + ep + bm + dm + fm);
            fun_bak[index] = fvalue - dt * sdp * (ww - 1.0);
          }
          if (sdm < 0.0) {
            am = am * am;
            cm = cm * cm;
            em = em * em;
            bp = bp * bp;
            dp = dp * dp;
            fp = fp * fp;
            // Engquist-Osher scheme
            ww = sqrt(am + cm + em + bp + dp + fp);
            fun_bak[index] = fvalue - dt * sdm * (ww - 1.0);
          }
          /*if(r < 0.4)
          {
                  printf("t = %d, i , j, k = %d, %d, %d, grad1 = %f S(d) = %f,
          fun = %f, ww = %f\n", t, i, j, k, grad1, sd, fvalue, ww);
          }*/
        }
      }
    }
    for (i = 0; i < nx * ny * nz; i++) {
      fun[i] = fun_bak[i];
    }
  }
  for (i = 0; i < nx * ny * nz; i++) {
    coeff[i] = fun[i];
  }
  ConvertToInterpolationCoefficients_3D(coeff, nx, ny, nz, CVC_DBL_EPSILON);
}

void ReInitilazation_Upwind_Eno_Godunov(float *fun, float *fun_bak,
                                        float *coeff, float dx, float dy,
                                        float dz, int nx, int ny, int nz,
                                        int tt)
// float *fun /* the volume data of initial signed distance function */
// float *coeff /* spline coefficients */
// float dx /* x direction spacing */
// float dy /* y direction spacing */
// float dz /* z direction spacing */
// int nx, /* number of points in x direction */
// int ny, /* number of points in y direction */
// int nz, /* number of points in z direction */
{
  cout << " In ReInitilazation_Upwind_Eno_Godunov" << endl;
  float x, y, z, fx[7], fy[7], fz[7];
  float a, b, c, d, e, f;
  for (int t = 0; t < 20; t++) {
    float maxf = -100000.0;
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
          int index = (i * ny + j) * nz + k;
          float fvalue = fun[index];
          fun_bak[index] = 0.0;
          Get_Seven_Function_Values_X(fun, nx, ny, nz, i, j, k, fx);
          Get_Seven_Function_Values_Y(fun, nx, ny, nz, i, j, k, fy);
          Get_Seven_Function_Values_Z(fun, nx, ny, nz, i, j, k, fz);
          if (fabs(fx[2]) < Height || fabs(fx[4]) < Height ||
              fabs(fy[2]) < Height || fabs(fy[4]) < Height ||
              fabs(fz[2]) < Height || fabs(fz[4]) < Height ||
              fabs(fvalue) < Height) {
            Quadr_Eno_Interpolation(fx, &b, &a, dx);
            Quadr_Eno_Interpolation(fy, &d, &c, dy);
            Quadr_Eno_Interpolation(fz, &f, &e, dz);
            float grad = 0.5 * (a * a + c * c + e * e + b * b + d * d + f * f);
            float sd = fvalue / sqrt(fvalue * fvalue + grad * dx * dx);
            float sdp = sd;
            if (sd < 0.0) {
              sdp = 0.0;
            }
            float sdm = sd;
            if (sd > 0.0) {
              sdm = 0.0;
            }
            x = -5.0 + i * dx;
            y = -5.0 + j * dy;
            z = -5.0 + k * dz;
            float r = sqrt(x * x + y * y + z * z);
            float fvalue_new = 0;
            float dt = 0.5 * dx;
            if (sdp >= 0.0) {
              float ww = sqrt(Extreme_Positive(a, b) + Extreme_Positive(c, d) +
                              Extreme_Positive(e, f));
              fvalue_new = -dt * sdp * (ww - 1.0); // Zhao's method
            }
            if (sdm < 0.0) {
              float ww = sqrt(Extreme_Negative(a, b) + Extreme_Negative(c, d) +
                              Extreme_Negative(e, f));
              fvalue_new = -dt * sdm * (ww - 1.0); // Zhao's method
            }
            fun_bak[index] = fvalue_new;
            if (fabs(fun_bak[index]) > maxf) {
              maxf = fabs(fun_bak[index]);
            }
          }
        }
      }
    }
    float maxt = dx / maxf;
    if (maxt > t) {
      maxt = t;
    }
    for (int i = 0; i < nx * ny * nz; i++) {
      fun[i] = fun[i] + maxt * fun_bak[i];
      // cut the function by Height
      if (fun[i] < -Height) {
        fun[i] = -Height1;
      }
      if (fun[i] > Height) {
        fun[i] = Height1;
      }
    }
  }
  for (int i = 0; i < nx * ny * nz; i++) {
    coeff[i] = fun[i];
  }
  ConvertToInterpolationCoefficients_3D(coeff, nx, ny, nz, CVC_DBL_EPSILON);
}

void ReInitilazation_Upwind_Eno_Godunov_Xu(float *fun, float *fun_bak,
                                           float *coeff, float dx, float dy,
                                           float dz, int nx, int ny, int nz,
                                           int tt)
// float *fun /* the volume data of initial signed distance function */
// float *coeff /* spline coefficients */
// float dx /* x direction spacing */
// float dy /* y direction spacing */
// float dz /* z direction spacing */
// int nx, /* number of points in x direction */
// int ny, /* number of points in y direction */
// int nz, /* number of points in z direction */
{

  // anchor
  float x, y, z, r, fx[7], fy[7], fz[7];
  float maxt;
  float a, b, c, d, e, f;
  float dt = 0.5 * dx;
  int last = nx * ny * nz;
  for (int t = 0; t < 20; t++) {
    int index = 0;
    float maxf = -100000.0;
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
          float fvalue = fun[index];
          fun_bak[index] = 0.0;
          Get_Seven_Function_Values_X(fun, nx, ny, nz, i, j, k, fx);
          Get_Seven_Function_Values_Y(fun, nx, ny, nz, i, j, k, fy);
          Get_Seven_Function_Values_Z(fun, nx, ny, nz, i, j, k, fz);
          if (fabs(fx[2]) < Height || fabs(fx[4]) < Height ||
              fabs(fy[2]) < Height || fabs(fy[4]) < Height ||
              fabs(fz[2]) < Height || fabs(fz[4]) < Height ||
              fabs(fvalue) < Height) {
            Quadr_Eno_Interpolation(fx, &b, &a, dx);
            Quadr_Eno_Interpolation(fy, &d, &c, dy);
            Quadr_Eno_Interpolation(fz, &f, &e, dz);
            float grad = 0.5 * (a * a + c * c + e * e + b * b + d * d + f * f);
            float sd = fvalue / sqrt(fvalue * fvalue + grad * dx * dx);
            float sdp = sd;
            if (sd < 0.0) {
              sdp = 0.0;
            }
            float sdm = sd;
            if (sd > 0.0) {
              sdm = 0.0;
            }
            // x = -5.0 + i*dx;
            // y = -5.0 + j*dy;
            // z = -5.0 + k*dz;
            if (sdp >= 0.0) {
              float ww = Extreme_Positive(a, b) + Extreme_Positive(c, d) +
                         Extreme_Positive(e, f);
              fun_bak[index] = -0.7 * dt * sdp * (ww - 1.0); // Zhao's method
            }
            if (sdm < 0.0) {
              float ww = Extreme_Negative(a, b) + Extreme_Negative(c, d) +
                         Extreme_Negative(e, f);
              fun_bak[index] = -0.7 * dt * sdm * (ww - 1.0); // Zhao's method
            }
            if (fabs(fun_bak[index]) > maxf) {
              maxf = fabs(fun_bak[index]);
            }
          }
          index++;
        }
      }
    }
    maxt = dx / maxf;
    if (maxt > t) {
      maxt = t;
    }
    for (int i = 0; i < last; i++) {
      fun[i] = fun[i] + maxt * fun_bak[i];
      // cut the function by Height
      if (fun[i] < -Height) {
        fun[i] = -Height1;
      }
      if (fun[i] > Height) {
        fun[i] = Height1;
      }
    }
  }
  for (int i = 0; i < last; i++) {
    coeff[i] = fun[i];
  }
  ConvertToInterpolationCoefficients_3D(coeff, nx, ny, nz, CVC_DBL_EPSILON);
}

void DiviededDifferencing_3j(float *f0, float *f1, float *f2, float *f3,
                             float dx) {
  int i;
  float indx1, indx2, indx3;
  indx1 = 1.0 / dx;
  // The first order divided difference
  for (i = 0; i < 6; i++) {
    f1[i] = (f0[i + 1] - f0[i]) * indx1;
  }
  // The second order divided difference
  indx2 = indx1 / 2.0;
  for (i = 0; i < 5; i++) {
    f2[i] = (f1[i + 1] - f1[i]) * indx2;
  }
  // The third order divided difference
  indx3 = indx1 / 3.0;
  for (i = 0; i < 4; i++) {
    f3[i] = (f2[i + 1] - f2[i]) * indx3;
  }
}

void DiviededDifferencing_2j(float *f0, float *f1, float *f2, float dx) {
  float indx1 = 1.0 / dx;
  // The first order divided difference
  for (int i = 1; i < 5; i++) {
    f1[i] = (f0[i + 1] - f0[i]) * indx1;
  }
  // The second order divided difference
  float indx2 = indx1 / 2.0;
  for (int i = 1; i < 4; i++) {
    f2[i] = (f1[i + 1] - f1[i]) * indx2;
  }
}

void DiviededDifferencing_1j(float *f0, float *f1, float dx) {
  float indx1, indx2, indx3;
  indx1 = 1.0 / dx;
  // The first order divided difference
  for (int i = 2; i < 4; i++) {
    f1[i] = (f0[i + 1] - f0[i]) * indx1;
  }
}

void Cubic_Eno_Interpolation(float *f0, float *u_plus, float *u_minus,
                             float dx) {
  int k;
  float f1[6], f2[5], f3[4];
  float x1, x2, x3;
  DiviededDifferencing_3j(f0, f1, f2, f3, dx);
  // Compute u^-
  k = 2;
  if (fabs(f2[k - 1]) < fabs(f2[k])) {
    k = k - 1;
  }
  if (fabs(f3[k - 1]) < fabs(f3[k])) {
    k = k - 1;
  }
  x1 = (3 - k) * dx;
  x2 = (2 - k) * dx;
  x3 = (1 - k) * dx;
  *u_minus = f1[k] + f2[k] * (x1 + x2) + f3[k] * (x1 * x2 + x1 * x3 + x2 * x3);
  // Compute u^+
  k = 3;
  if (fabs(f2[k - 1]) < fabs(f2[k])) {
    k = k - 1;
  }
  if (fabs(f3[k - 1]) < fabs(f3[k])) {
    k = k - 1;
  }
  x1 = (3 - k) * dx;
  x2 = (2 - k) * dx;
  x3 = (1 - k) * dx;
  *u_plus = f1[k] + f2[k] * (x1 + x2) + f3[k] * (x1 * x2 + x1 * x3 + x2 * x3);
}

void Quadr_Eno_Interpolation(float *f0, float *u_plus, float *u_minus,
                             float dx) {
  float f1[6], f2[5];
  DiviededDifferencing_2j(f0, f1, f2, dx);
  float temp = fabs(f2[2]);
  if (fabs(f2[1]) < temp) {
    *u_minus = f1[1] + f2[1] * dx * 3;
  } else {
    *u_minus = f1[2] + f2[2] * dx;
  }
  if (temp < fabs(f2[3])) {
    *u_plus = f1[2] + f2[2] * dx;
  } else {
    *u_plus = f1[3] + f2[3] * (-dx);
  }
}

void Linear_Eno_Interpolation(float *f0, float *u_plus, float *u_minus,
                              float dx) {
  int k;
  float f1[6];
  float x1, x2, x3;
  DiviededDifferencing_1j(f0, f1, dx);
  // Compute u^-
  k = 2;
  *u_minus = f1[k];
  // Compute u^+
  k = 3;
  *u_plus = f1[k];
}

void inline Get_Seven_Function_Values_X(float *coeff, int nx, int ny, int nz,
                                        int u, int v, int w, float *fx) {
  if (u >= 3 && u <= nx - 4) {
    for (int j = 0; j < 7; j++) {
      int l = u - 3 + j;
      fx[j] = coeff[(l * ny + v) * nz + w];
    }
    return;
  }
  switch (u) {
  case 2: {
    for (int j = 1; j < 7; j++) {
      int l = j - 1;
      fx[j] = coeff[(l * ny + v) * nz + w];
    }
    fx[0] = 2 * fx[1] - fx[2];
    return;
  }
  case 1: {
    for (int j = 2; j < 7; j++) {
      int l = j - 2;
      fx[j] = coeff[(l * ny + v) * nz + w];
    }
    fx[1] = 2 * fx[2] - fx[3];
    fx[0] = 2 * fx[1] - fx[2];
    return;
  }
  case 0: {
    for (int j = 3; j < 7; j++) {
      int l = j - 3;
      fx[j] = coeff[(l * ny + v) * nz + w];
    }
    fx[2] = 2 * fx[3] - fx[4];
    fx[1] = 2 * fx[2] - fx[3];
    fx[0] = 2 * fx[1] - fx[2];
    return;
  }
  }
  switch (nx - u) {
  case 3: {
    for (int j = 0; j < 6; j++) {
      int l = u - 3 + j;
      fx[j] = coeff[(l * ny + v) * nz + w];
    }
    fx[6] = 2 * fx[5] - fx[4];
    return;
  }
  case 2: {
    for (int j = 0; j < 5; j++) {
      int l = u - 3 + j;
      fx[j] = coeff[(l * ny + v) * nz + w];
    }
    fx[5] = 2 * fx[4] - fx[3];
    fx[6] = 2 * fx[5] - fx[4];
    return;
  }
  case 1: {
    for (int j = 0; j < 4; j++) {
      int l = u - 3 + j;
      fx[j] = coeff[(l * ny + v) * nz + w];
    }
    fx[4] = 2 * fx[3] - fx[2];
    fx[5] = 2 * fx[4] - fx[3];
    fx[6] = 2 * fx[5] - fx[4];
    return;
  }
  }
}

void inline Get_Seven_Function_Values_Y(float *coeff, int nx, int ny, int nz,
                                        int u, int v, int w, float *fx) {
  int j, l, uj;
  if (v >= 3 && v <= ny - 4) {
    for (j = 0; j < 7; j++) {
      l = v - 3 + j;
      fx[j] = coeff[(u * ny + l) * nz + w];
    }
    return;
  }
  if (v == 2) {
    for (j = 1; j < 7; j++) {
      l = j - 1;
      fx[j] = coeff[(u * ny + l) * nz + w];
    }
    fx[0] = 2 * fx[1] - fx[2];
    return;
  }
  if (v == 1) {
    for (j = 2; j < 7; j++) {
      l = j - 2;
      fx[j] = coeff[(u * ny + l) * nz + w];
    }
    fx[1] = 2 * fx[2] - fx[3];
    fx[0] = 2 * fx[1] - fx[2];
    return;
  }
  if (v == 0) {
    for (j = 3; j < 7; j++) {
      l = j - 3;
      fx[j] = coeff[(u * ny + l) * nz + w];
    }
    fx[2] = 2 * fx[3] - fx[4];
    fx[1] = 2 * fx[2] - fx[3];
    fx[0] = 2 * fx[1] - fx[2];
    return;
  }
  if (v == ny - 3) {
    for (j = 0; j < 6; j++) {
      l = v - 3 + j;
      fx[j] = coeff[(u * ny + l) * nz + w];
    }
    fx[6] = 2 * fx[5] - fx[4];
    return;
  }
  if (v == ny - 2) {
    for (j = 0; j < 5; j++) {
      l = v - 3 + j;
      fx[j] = coeff[(u * ny + l) * nz + w];
    }
    fx[5] = 2 * fx[4] - fx[3];
    fx[6] = 2 * fx[5] - fx[4];
    return;
  }
  if (v == ny - 1) {
    for (j = 0; j < 4; j++) {
      l = v - 3 + j;
      fx[j] = coeff[(u * ny + l) * nz + w];
    }
    fx[4] = 2 * fx[3] - fx[2];
    fx[5] = 2 * fx[4] - fx[3];
    fx[6] = 2 * fx[5] - fx[4];
    return;
  }
}

void inline Get_Seven_Function_Values_Z(float *coeff, int nx, int ny, int nz,
                                        int u, int v, int w, float *fx) {
  int j, l, uj;
  if (w >= 3 && w <= nz - 4) {
    for (j = 0; j < 7; j++) {
      l = w - 3 + j;
      fx[j] = coeff[(u * ny + v) * nz + l];
    }
    return;
  }
  if (w == 2) {
    for (j = 1; j < 7; j++) {
      l = j - 1;
      fx[j] = coeff[(u * ny + v) * nz + l];
    }
    fx[0] = 2 * fx[1] - fx[2];
    return;
  }
  if (w == 1) {
    for (j = 2; j < 7; j++) {
      l = j - 2;
      fx[j] = coeff[(u * ny + v) * nz + l];
    }
    fx[1] = 2 * fx[2] - fx[3];
    fx[0] = 2 * fx[1] - fx[2];
    return;
  }
  if (w == 0) {
    for (j = 3; j < 7; j++) {
      l = j - 3;
      fx[j] = coeff[(u * ny + v) * nz + l];
    }
    fx[2] = 2 * fx[3] - fx[4];
    fx[1] = 2 * fx[2] - fx[3];
    fx[0] = 2 * fx[1] - fx[2];
    return;
  }
  if (w == nz - 3) {
    for (j = 0; j < 6; j++) {
      l = w - 3 + j;
      fx[j] = coeff[(u * ny + v) * nz + l];
    }
    fx[6] = 2 * fx[5] - fx[4];
    return;
  }
  if (w == nz - 2) {
    for (j = 0; j < 5; j++) {
      l = w - 3 + j;
      fx[j] = coeff[(u * ny + v) * nz + l];
    }
    fx[5] = 2 * fx[4] - fx[3];
    fx[6] = 2 * fx[5] - fx[4];
    return;
  }
  if (w == nz - 1) {
    for (j = 0; j < 4; j++) {
      l = w - 3 + j;
      fx[j] = coeff[(u * ny + v) * nz + l];
    }
    fx[4] = 2 * fx[3] - fx[2];
    fx[5] = 2 * fx[4] - fx[3];
    fx[6] = 2 * fx[5] - fx[4];
    return;
  }
}

float Extreme_Positive(float a, float b) {
  float aa, bb;
  aa = a * a;
  bb = b * b;
  if (a <= b) {
    if (b < 0.0) {
      return (bb);
    }
    if (a > 0.0) {
      return (aa);
    }
    return (0.0);
  }
  if (aa > bb) {
    return (aa);
  }
  return (bb);
}

float Extreme_Negative(float a, float b) {
  float aa, bb;
  aa = a * a;
  bb = b * b;
  if (a <= b) {
    if (aa > bb) {
      return (aa);
    }
    return (bb);
  }
  if (a < 0.0) {
    return (aa);
  }
  if (b > 0.0) {
    return (bb);
  }
  return (0.0);
}

float Gradient_2(float *fun, float dx, float dy, float dz, int i, int j, int k,
                 int nx, int ny, int nz) {
  float fx, fy, fz;
  if (i > 0 && i < nx - 1) {
    fx = (fun[((i + 1) * ny + j) * nz + k] - fun[((i - 1) * ny + j) * nz + k]) /
         (dx + dx);
  }
  if (i == 0) {
    fx = (fun[((1) * ny + j) * nz + k] - fun[j * nz + k]) / dx;
  }
  if (i == nx - 1) {
    fx = (fun[((nx - 1) * ny + j) * nz + k] -
          fun[((nx - 2) * ny + j) * nz + k]) /
         dx;
  }
  if (j > 0 && j < ny - 1) {
    fy = (fun[(i * ny + j + 1) * nz + k] - fun[(i * ny + j - 1) * nz + k]) /
         (dy + dy);
  }
  if (j == 0) {
    fy = (fun[(i * ny + 1) * nz + k] - fun[(i * ny) * nz + k]) / dy;
  }
  if (j == ny - 1) {
    fy = (fun[(i * ny + ny - 1) * nz + k] - fun[(i * ny + ny - 2) * nz + k]) /
         dy;
  }
  if (k > 0 && k < nz - 1) {
    fz = (fun[(i * ny + j) * nz + k + 1] - fun[(i * ny + j) * nz + k - 1]) /
         (dz + dz);
  }
  if (k == 0) {
    fz = (fun[(i * ny + j) * nz + 1] - fun[(i * ny + j) * nz]) / dz;
  }
  if (k == nz - 1) {
    fz = (fun[(i * ny + j) * nz + nz - 1] - fun[(i * ny + j) * nz + nz - 2]) /
         dz;
  }
  return (fx * fx + fy * fy + fz * fz);
}

void Divided_DifferenceOrder2PartialsAtGridPoint(float *c, float dx, float dy,
                                                 float dz, int nx, int ny,
                                                 int nz, int u, int v, int w,
                                                 float *partials)
// float *c, /* the spline coefficients */
// float dx /* spacing in x direction */
// float dy /* spacing in y direction */
// float dz /* spacing in z direction */
// int nx, /* number of samples or coefficients in x direction */
// int ny, /* number of samples or coefficients in y direction */
// int nz, /* number of samples or coefficients in z direction */
// float *partial /* partial derivatives */
{
  float c27[27], indx, indy, indz;
  int i;
  Take_27_Coefficients(c, nx, ny, nz, u, v, w, c27);
  indx = 1.0 / dx;
  indy = 1.0 / dy;
  indz = 1.0 / dz;
  partials[0] = c27[13];
  // fx, fy, fz
  partials[1] = (c27[22] - c27[4]) * indx * 0.5;
  partials[2] = (c27[16] - c27[10]) * indy * 0.5;
  partials[3] = (c27[14] - c27[12]) * indz * 0.5;
  // fxx, fxy, fxz
  partials[4] = (c27[22] + c27[4] - c27[13] - c27[13]) * indx * indx;
  partials[5] = (c27[25] + c27[1] - c27[19] - c27[7]) * indx * indy * 0.25;
  partials[6] = (c27[23] + c27[3] - c27[21] - c27[5]) * indx * indz * 0.25;
  // fyy, fyz
  partials[7] = (c27[16] + c27[10] - c27[13] - c27[13]) * indy * indy;
  partials[8] = (c27[17] + c27[9] - c27[15] - c27[11]) * indy * indz * 0.25;
  // fzz
  partials[9] = (c27[12] + c27[14] - c27[13] - c27[13]) * indz * indz;
}

void ConvertToInterpolationCoefficients_Qu(
    float *c,       /* input samples --> output coefficients */
    int DataLength, /* number of samples or coefficients */
    float *z,       /* poles */
    int NbPoles,    /* number of poles */
    float Tolerance /* admissible relative error */
) {
  /* begin ConvertToInterpolationCoefficients */
  float Lambda = 1.0;
  int n, k;
  /* special case required by mirror boundaries */
  if (DataLength == 1) {
    return;
  }
  /* compute the overall gain */
  for (k = 0; k < NbPoles; k++) {
    Lambda = Lambda * (1.0f - z[k]) * (1.0f - 1.0f / z[k]);
  }
  /* apply the gain */
  for (n = 0; n < DataLength; n++) {
    c[n] *= Lambda;
  }
  /* loop over all poles */
  for (k = 0; k < NbPoles; k++) {
    /* causal initialization */
    c[0] = InitialCausalCoefficient(c, DataLength, z[k], Tolerance);
    /* causal recursion */
    for (n = 1; n < DataLength; n++) {
      c[n] += z[k] * c[n - 1];
    }
    /* anticausal initialization */
    c[DataLength - 1] = InitialAntiCausalCoefficient(c, DataLength, z[k]);
    /* anticausal recursion */
    for (n = DataLength - 2; 0 <= n; n--) {
      c[n] = z[k] * (c[n + 1] - c[n]);
    }
  }
}

void Evaluat_Four_Basis(float x, float *values) {
  float y;
  values[0] = (1.0 - x) * (1.0 - x) * (1.0 - x) / 6.0;
  values[1] = 0.666666666 - x * x + 0.5 * x * x * x;
  y = 1.0 - x;
  values[2] = 0.666666666 - y * y + 0.5 * y * y * y;
  y = 2.0 - x;
  values[3] = x * x * x / 6.0;
}

void Evaluate_Four_Basis_First_Partial(float x, float *firstPartials) {
  float y;
  firstPartials[0] = -0.5 * (1.0 - x) * (1.0 - x);
  firstPartials[1] = 1.5 * x * x - 2 * x;
  y = 1.0 - x;
  firstPartials[2] = -1.5 * y * y + 2 * y;
  firstPartials[3] = 0.5 * x * x;
}

void Evaluate_Four_Basis_Second_Partial(float x, float *secondPartials) {
  secondPartials[0] = 1.0 - x;
  secondPartials[1] = 3 * x - 2;
  secondPartials[2] = 1 - 3 * x;
  secondPartials[3] = x;
}

void ReSamplingCubicSpline1OrderPartials(float *coeff, int nx, int ny, int nz,
                                         float *dxyz, float *minExtent,
                                         float *zeroToTwoPartialValue, float *p)
/* p the point to be evaluated on */
/* coeff the spline coefficients */
/* firstParValue, 3 directional derivatives */
{
  int i, j, k, l;
  int u, v, w;
  float x, y, z;
  float funcValue_x[4], funcValue_y[4], funcValue_z[4];
  float f1Partial_x[4], f1Partial_y[4], f1Partial_z[4], f2Partial_xx[4],
      f2Partial_xy[4], f2Partial_xz[4], f2Partial_yy[4], f2Partial_yz[4],
      f2Partial_zz[4];
  float c64[64], tensor[64], tensorx[64], tensory[64], tensorz[64],
      tensorxx[64], tensorxy[64], tensorxz[64], tensoryy[64], tensoryz[64],
      tensorzz[64];
  float result[10];
  float dx = dxyz[2];
  float dy = dxyz[1];
  float dz = dxyz[0];
  float indx = 1.0 / dx;
  float indy = 1.0 / dy;
  float indz = 1.0 / dz;
  u = (int)((p[0] - minExtent[0]) * indx);
  v = (int)((p[1] - minExtent[1]) * indy);
  w = (int)((p[2] - minExtent[2]) * indz);
  float fu = (p[0] - minExtent[0]) / dx;
  float fv = (p[1] - minExtent[1]) / dy;
  float fw = (p[2] - minExtent[2]) / dz;
  if (fu - u > 0.999) {
    u = u + 1;
  }
  if (fv - v > 0.999) {
    v = v + 1;
  }
  if (fw - w > 0.999) {
    w = w + 1;
  }
  x = (p[0] - minExtent[0] - u * dx) / dx;
  y = (p[1] - minExtent[1] - v * dy) / dy;
  z = (p[2] - minExtent[2] - w * dz) / dz;
  if (fabs(x) < 0.0001) {
    x = 0.0;
  }
  if (fabs(y) < 0.0001) {
    y = 0.0;
  }
  if (fabs(z) < 0.0001) {
    z = 0.0;
  }
  if (x > 0.9999) {
    x = 1.0;
  }
  if (y > 0.9999) {
    y = 1.0;
  }
  if (z > 0.9999) {
    z = 1.0;
  }
  Evaluat_Four_Basis(x, funcValue_x);
  Evaluat_Four_Basis(y, funcValue_y);
  Evaluat_Four_Basis(z, funcValue_z);
  Evaluate_Four_Basis_First_Partial(x, f1Partial_x);
  Evaluate_Four_Basis_First_Partial(y, f1Partial_y);
  Evaluate_Four_Basis_First_Partial(z, f1Partial_z);
  Evaluate_Four_Basis_Second_Partial(x, f2Partial_xx);
  Evaluate_Four_Basis_Second_Partial(y, f2Partial_yy);
  Evaluate_Four_Basis_Second_Partial(z, f2Partial_zz);
  Tensor_444(funcValue_x, funcValue_y, funcValue_z, tensor);
  Tensor_444(f1Partial_x, funcValue_y, funcValue_z, tensorx);
  Tensor_444(funcValue_x, f1Partial_y, funcValue_z, tensory);
  Tensor_444(funcValue_x, funcValue_y, f1Partial_z, tensorz);
  Tensor_444(f2Partial_xx, funcValue_y, funcValue_z, tensorxx);
  Tensor_444(f1Partial_x, f1Partial_y, funcValue_z, tensorxy);
  Tensor_444(f1Partial_x, funcValue_y, f1Partial_z, tensorxz);
  Tensor_444(funcValue_x, f2Partial_yy, funcValue_z, tensoryy);
  Tensor_444(funcValue_x, f1Partial_y, f1Partial_z, tensoryz);
  Tensor_444(funcValue_x, funcValue_y, f2Partial_zz, tensorzz);
  if ((u > 0 && u < nx - 2) && (v > 0 && v < ny - 2) && (w > 0 && w < nz - 2)) {
    Take_64_Coefficients(coeff, nx, ny, nz, w, v, u, c64);
    for (l = 0; l < 10; l++) {
      result[l] = 0.0;
    }
    for (l = 0; l < 64; l++) {
      result[0] += tensor[l] * c64[l];
      result[1] += tensorx[l] * c64[l];  //*indx;
      result[2] += tensory[l] * c64[l];  //*indy;
      result[3] += tensorz[l] * c64[l];  //*indz;
      result[4] += tensorxx[l] * c64[l]; //*indx*indx;
      result[5] += tensorxy[l] * c64[l]; //*indx*indy;
      result[6] += tensorxz[l] * c64[l]; //*indx*indz;
      result[7] += tensoryy[l] * c64[l]; //*indy*indy;
      result[8] += tensoryz[l] * c64[l]; //*indy*indz;
      result[9] += tensorzz[l] * c64[l]; //*indz*indz;
    }
    for (l = 0; l < 10; l++) {
      zeroToTwoPartialValue[l] = result[l];
    }
  }
}

void ReSamplingCubicSpline1OrderPartials2(float *coeff, int nx, int ny, int nz,
                                          float *dxyz, float *minExtent,
                                          float *zeroToTwoPartialValue,
                                          float *p)
/* p the point to be evaluated on */
/* coeff the spline coefficients */
/* firstParValue, 3 directional derivatives */
// finite difference method to compute derivatices//
{
  int i, j, k, l;
  int u, v, w;
  float x, y, z;
  float funcValue_x[4], funcValue_y[4], funcValue_z[4], fval0[4], fval1[4];
  float f1Partial_x[4], f1Partial_y[4], f1Partial_z[4], f2Partial_xx[4],
      f2Partial_xy[4], f2Partial_xz[4], f2Partial_yy[4], f2Partial_yz[4],
      f2Partial_zz[4];
  float c64[64], tensor[64], tensorx[64], tensory[64], tensorz[64],
      tensorxx[64], tensorxy[64], tensorxz[64], tensoryy[64], tensoryz[64],
      tensorzz[64];
  float resultfun, resultx, resulty, resultz, resultxx, resultxy, resultxz,
      resultyy, resultyz, resultzz;
  float fvalA[64], fvalB[64], fvalC[64], fvalD[64], fvalE[64], fvalF[64],
      fvalG[64], fvalH[64];
  float resultA, resultB, resultC, resultD, resultE, resultF, resultG, resultH;
  float dx = dxyz[2];
  float dy = dxyz[1];
  float dz = dxyz[0];
  w = (int)((p[0] - minExtent[0]) / dx);
  v = (int)((p[1] - minExtent[1]) / dy);
  u = (int)((p[2] - minExtent[2]) / dz);
  z = (p[0] - minExtent[0] - u * dx) / dx;
  y = (p[1] - minExtent[1] - v * dy) / dy;
  x = (p[2] - minExtent[2] - w * dz) / dz;
  Evaluat_Four_Basis(0, fval0);
  Evaluat_Four_Basis(1, fval1);
  if ((u > 0 && u < nx - 2) && (v > 0 && v < ny - 2) && (w > 0 && w < nz - 2)) {
    Evaluat_Four_Basis(x, funcValue_x);
    Evaluat_Four_Basis(y, funcValue_y);
    Evaluat_Four_Basis(z, funcValue_z);
    Evaluate_Four_Basis_First_Partial(x, f1Partial_x);
    Evaluate_Four_Basis_First_Partial(y, f1Partial_y);
    Evaluate_Four_Basis_First_Partial(z, f1Partial_z);
    Evaluate_Four_Basis_Second_Partial(x, f2Partial_xx);
    Evaluate_Four_Basis_Second_Partial(y, f2Partial_yy);
    Evaluate_Four_Basis_Second_Partial(z, f2Partial_zz);
    Tensor_444(funcValue_x, funcValue_y, funcValue_z, tensor);
    Tensor_444(f1Partial_x, funcValue_y, funcValue_z, tensorx);
    Tensor_444(funcValue_x, f1Partial_y, funcValue_z, tensory);
    Tensor_444(funcValue_x, funcValue_y, f1Partial_z, tensorz);
    Tensor_444(fval0, funcValue_y, funcValue_z, fvalA);
    Tensor_444(fval1, funcValue_y, funcValue_z, fvalB);
    Tensor_444(funcValue_x, fval0, funcValue_z, fvalC);
    Tensor_444(funcValue_x, fval1, funcValue_z, fvalD);
    Tensor_444(funcValue_x, funcValue_y, fval0, fvalE);
    Tensor_444(funcValue_x, funcValue_y, fval1, fvalF);
    Take_64_Coefficients(coeff, nx, ny, nz, u, v, w, c64);
    resultfun = 0.0;
    resultx = 0.0;
    resulty = 0.0;
    resultz = 0.0;
    resultA = 0.0;
    resultB = 0.0;
    resultC = 0.0;
    resultD = 0.0;
    resultE = 0.0;
    resultF = 0.0;
    for (l = 0; l < 64; l++) {
      resultfun = resultfun + tensor[l] * c64[l];
      resultx = resultx + tensorx[l] * c64[l];
      resulty = resulty + tensory[l] * c64[l];
      resultz = resultz + tensorz[l] * c64[l];
      resultA += fvalA[l] * c64[l];
      resultB += fvalB[l] * c64[l];
      resultC += fvalC[l] * c64[l];
      resultD += fvalD[l] * c64[l];
      resultE += fvalE[l] * c64[l];
      resultF += fvalF[l] * c64[l];
    }
    zeroToTwoPartialValue[0] = resultfun;
    if (((1 - x) >= 0.0001) && (x >= 0.0001)) {
      zeroToTwoPartialValue[1] =
          0.5 * ((resultB - resultfun) / (1 - x) + (resultfun - resultA) / x);
    } else if (x < 0.0001) {
      zeroToTwoPartialValue[1] = (resultB - resultfun) / (1 - x);
    } else {
      zeroToTwoPartialValue[1] = (resultfun - resultA) / x;
    }
    if (((1 - y) >= 0.0001) && (y >= 0.0001)) {
      zeroToTwoPartialValue[2] =
          0.5 * ((resultD - resultfun) / (1 - y) + (resultfun - resultC) / y);
    } else if (y < 0.0001) {
      zeroToTwoPartialValue[2] = (resultD - resultfun) / (1 - y);
    } else {
      zeroToTwoPartialValue[2] = (resultfun - resultC) / y;
    }
    if (((1 - z) >= 0.0001) && (z >= 0.0001)) {
      zeroToTwoPartialValue[3] =
          0.5 * ((resultF - resultfun) / (1 - z) + (resultfun - resultE) / z);
    } else if (z < 0.0001) {
      zeroToTwoPartialValue[3] = (resultF - resultfun) / (1 - z);
    } else {
      zeroToTwoPartialValue[3] = (resultfun - resultE) / z;
    }
    return;
  }
}

void ReSamplingCubicSpline(float *coeff, int nx, int ny, int nz,
                           float *funvalues, int Nx, int Ny, int Nz) {
  int i, j, k, l, ix, jy, kz;
  float Dx, Dy, Dz, x, y, z;
  float values_x[4], values_y[4], values_z[4];
  float c64[64], tensor[64], result;
  Dx = (nx - 1.0) / (Nx - 1.0);
  Dy = (ny - 1.0) / (Ny - 1.0);
  Dz = (nz - 1.0) / (Nz - 1.0);
  for (i = 0; i < Nx; i++) {
    x = i * Dx;
    ix = x;
    x = x - ix;
    for (j = 0; j < Ny; j++) {
      y = j * Dy;
      jy = y;
      y = y - jy;
      for (k = 0; k < Nz; k++) {
        z = k * Dz;
        kz = z;
        z = z - kz;
        Evaluat_Four_Basis(x, values_x);
        Evaluat_Four_Basis(y, values_y);
        Evaluat_Four_Basis(z, values_z);
        Tensor_444(values_x, values_y, values_z, tensor);
        Take_64_Coefficients(coeff, nx, ny, nz, ix, jy, kz, c64);
        result = 0.0;
        for (l = 0; l < 64; l++) {
          result = result + tensor[l] * c64[l];
        }
        funvalues[(i * Ny + j) * Nz + k] = result;
      }
    }
  }
}

void Fast_Evaluate_Gauss_Map_One_Atom(float *corner, float *dxyz, float *center,
                                      float radius, unsigned int *nxyz, float C,
                                      float *data, float epsilon) {
  float a, b, c, dx, dy, dz, cx, cy, cz, x, y, z;
  float rr2, rr, xmin, ymin, zmin, w;
  int i, j, k, nx, ny, nz, wx, wy, wz, ix, iy, iz;
  a = corner[0];
  b = corner[1];
  c = corner[2];
  dx = 1.0 / dxyz[0];
  dy = 1.0 / dxyz[1];
  dz = 1.0 / dxyz[2];
  cx = center[0];
  cy = center[1];
  cz = center[2];
  nx = nxyz[0];
  ny = nxyz[1];
  nz = nxyz[2];
  rr2 = radius * radius;
  rr = rr2 - log(epsilon) / C;
  rr = sqrt(rr);
  xmin = cx - rr;
  ymin = cy - rr;
  zmin = cz - rr;
  rr = rr + rr;
  wx = rr * dx;
  wy = rr * dy;
  wz = rr * dz;
  ix = (xmin - a) * dx + 1;
  iy = (ymin - b) * dy + 1;
  iz = (zmin - c) * dz + 1;
  for (i = ix; i <= ix + wx; i++) {
    if (i >= 0 && i < nx) {
      x = a + i * dxyz[0];
      for (j = iy; j <= iy + wy; j++) {
        if (j >= 0 && j < ny) {
          y = b + j * dxyz[1];
          for (k = iz; k <= iz + wz; k++) {
            if (k >= 0 && k < nz) {
              z = c + k * dxyz[2];
              w = exp(-C * ((x - cx) * (x - cx) + (y - cy) * (y - cy) +
                            (z - cz) * (z - cz) - rr2));
              data[(i * ny + j) * nz + k] = data[(i * ny + j) * nz + k] - w;
            }
          }
        }
      }
    }
  }
}
