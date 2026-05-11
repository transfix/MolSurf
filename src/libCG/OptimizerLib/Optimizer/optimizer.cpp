#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <iostream>
#include <libCG/OptimizerLib/Optimizer/optimizer.h>
#include <libCG/OptimizerLib/Optimizer/optimizationFunction.h>
#include <PDBParser/GroupOfAtoms.h>
#include <PDBParser/Atom.h>
#include <libCG/CoarseGrain/geometry.h>
#include <libCG/CoarseGrain/surface.h>
#include "../src/levmar-2.5/levmar.h"

#include <time.h>
#ifdef _WIN32
#include <sys/types.h>
#include <sys/timeb.h>
#else
#include <sys/time.h>
#endif

using namespace PDBParser;
using namespace std;
using namespace GEOMETRY;

Optimizer::Optimizer(std::vector<Atom *> &al, Atom *bd, Surface *surface) {
  vector<Atom *>::iterator it;
  for (it = al.begin(); it != al.end(); it++) {
    atomList.push_back(*it);
  }
  bead = bd;
  getSamplePoints(surface);
  numOfSamples = surface->numbpts;
  printf("numsamples:%d\n", numOfSamples);
  blurValue = new double[numOfSamples];
  adata[0].center[0] = bead->m_Position[0];
  adata[0].center[1] = bead->m_Position[1];
  adata[0].center[2] = bead->m_Position[2];
  adata[0].radius = bead->getRadius();
  m_blobby = -2.3442;
}

Optimizer::~Optimizer() {
  // delete blurValue;
  /*for (int i = 0; i < numOfSamples; i++)
          {
  printf("%e\n",samplePoints[i][0]);
          delete samplePoints[i];
          }
  delete samplePoints;*/
}

double Optimizer::getMyTime() {
#ifdef _WIN32
  time_t ltime;
  time(&ltime);
  return (double)ltime;
#else
  struct timeval t;
  gettimeofday(&t, NULL);
  return (double)(t.tv_sec + 1e-6 * t.tv_usec);
#endif
}

/*
std::vector<Atom*> Optimizer::getAtomList()
{
        vector<Atom*> atomList;
        for (int i = 0; i < bead->numOfChildren; i++)
        {
                int atomIndex = bead->child[i];
                Atom *atom = new Atom(*(molecule->m_Atoms.at(atomIndex)));
//              Atom *atom = Atom(molecule->m_Atoms[atomIndex]);
                atomList.push_back(atom);
        }
}
*/

/*
void generateGridPoints()
{
  double span = (maxPos - minPos)/(double)dim;
  for(int i=0;i<dim;i++)
    for(int j=0;j<dim;j++)
      for(int k=0;k<dim;k++)
      {
        int index = i*dim*dim + j*dim + k;
        samplePoints[index][0] = minPos + span * (double)i;
        samplePoints[index][1] = minPos + span * (double)j;
        samplePoints[index][2] = minPos + span * (double)k;
        //printf("index:%d\n",index);
      }
}
*/

void Optimizer::getSamplePoints(Surface *surface) {
  int i, j, N;
  N = surface->numbpts;
  // printf("N:%d\n",N);
  samplePoints = new double *[N];
  for (i = 0; i < N; i++) {
    samplePoints[i] = new double[3];
  }
  for (i = 0; i < N; i++)
    for (j = 0; j < 3; j++) {
      samplePoints[i][j] = surface->m_Vertices[i]->Point[j];
    }
}

double Optimizer::evalDensity(double blobby, double ctr[3], double pnt[3],
                              double radius) {
  double r = (ctr[0] - pnt[0]) * (ctr[0] - pnt[0]) +
             (ctr[1] - pnt[1]) * (ctr[1] - pnt[1]) +
             (ctr[2] - pnt[2]) * (ctr[2] - pnt[2]);
  double r0 = radius * radius;
  double expval = blobby * r / r0 - blobby; // 10/16/2008
  //  double expval = blobby*(r - r0);
  return exp(expval);
}

void Optimizer::blurAtomList(std::vector<Atom *> atomList, double *blurValue) {
  double x[3];
  int i, j;
  for (i = 0; i < numOfSamples; i++) {
    blurValue[i] = 0.0;
    for (j = 0; j < atomList.size(); j++) {
      Atom *atom = atomList[j];
      if (!atom) {
        continue;
      }
      // center and radii in angstroms
      x[0] = atomList[j]->m_Position[0];
      x[1] = atomList[j]->m_Position[1];
      x[2] = atomList[j]->m_Position[2];
      double r = atomList[j]->getRadius();
      blurValue[i] += evalDensity(m_blobby, x, samplePoints[i], r);
    }
  }
}

void Optimizer::initializeFunctionData(FunctionData *fdata) {
  fdata->adata[0].center[0] = adata[0].center[0];
  fdata->adata[0].center[1] = adata[0].center[1];
  fdata->adata[0].center[2] = adata[0].center[2];
  fdata->adata[0].radius = adata[0].radius;
  fdata->samplePoints = new double *[numOfSamples];
  fdata->blurValue = new double[numOfSamples];
  for (int i = 0; i < numOfSamples; i++) {
    fdata->samplePoints[i] = new double[3];
    fdata->samplePoints[i][0] = samplePoints[i][0];
    fdata->samplePoints[i][1] = samplePoints[i][1];
    fdata->samplePoints[i][2] = samplePoints[i][2];
    fdata->blurValue[i] = blurValue[i];
  }
  fdata->m_blobby = m_blobby;
  fdata->numOfSamples = numOfSamples;
}

void Optimizer::updateFunctionData(FunctionData *fdata, double *c, double p) {
  fdata->adata[0].center[0] = c[0];
  fdata->adata[0].center[1] = c[1];
  fdata->adata[0].center[2] = c[2];
  fdata->adata[0].radius = p;
}

void Optimizer::radiusOptimize(double *p, double *x, int m, int n, void *data) {
  register int i, j, k;
  atomData *dptr = (atomData *)data;
  /* 1D least square version (fast) */
  for (k = 0; k < n; k++) {
    x[k] = 0.0;
    for (j = 0; j < numOfSamples; j++) {
      double tmp = 0.0;
      for (i = 0; i < m; i++) {
        tmp += evalDensity(m_blobby, dptr[i].center, samplePoints[j],
                           p[i]); // - blurValue[j];
      }
      x[k] += (tmp - blurValue[j]) * (tmp - blurValue[j]);
    }
  }
  /* vector version, slow,  minima is close to the least square version
    for(j=0; j<n; j++)
    {
      x[j] = 0.0;
      for(i=0; i<m; i++)
      {
        x[j] += evalDensity(m_blobby, dptr[i].center, samplePoints[j], p[i]);//
    - blurValue[j];
      //x[i]=100.0*(p[1]-p[0]*p[0])*(p[1]-p[0]*p[0]) + (p[0]-1.0)*(p[0]-1.0);
      }
      //printf("x[%d] = %f\n",j,x[j]);
      x[j] -= blurValue[j];
    }
    //printf("x[71]=%f\n",x[71]);
  */
}

void Optimizer::centerOptimize(double *c, double *x, int mm, int n,
                               void *data) {
  register int i, j, k;
  atomData *dptr = (atomData *)data;
  double center[3];
  int m = mm / 3;
  /* 1D least square version (fast)*/
  for (k = 0; k < n; k++) {
    x[k] = 0.0;
    for (j = 0; j < numOfSamples; j++) {
      double tmp = 0.0;
      for (i = 0; i < m; i++) {
        center[0] = c[3 * i];
        center[1] = c[3 * i + 1];
        center[2] = c[3 * i + 2];
        tmp += evalDensity(m_blobby, center, samplePoints[j],
                           dptr[i].radius); // - blurValue[j];
      }
      x[k] += (tmp - blurValue[j]) * (tmp - blurValue[j]);
    }
  }
  /* vector version, slow, minima is close to the least square version
    for(j=0; j<n; j++)
    {
      x[j] = 0.0;
      for(i=0; i<m; i++)
      {
        center[0] = c[3*i];       center[1] = c[3*i+1];   center[2] = c[3*i+2];
        x[j] += evalDensity(m_blobby, center, samplePoints[j],
    dptr[i].radius);// - blurValue[j];
      }
      x[j] -= blurValue[j];
    }
  */
}

void Optimizer::radiuscenterOptimize(double *c, double *x, int m, int n,
                                     void *data) {
  register int i, j, k;
  double tmp = 0.0;
  //  atomData* dptr = (atomData *)data;
  double center[3];
  double radius;
  /* 1D least square version (fast)*/
  for (k = 0; k < n; k++) {
    x[k] = 0.0;
    for (j = 0; j < numOfSamples; j++) {
      double tmp = 0.0;
      for (i = 0; i < m; i++) {
        center[0] = c[3 * i];
        center[1] = c[3 * i + 1];
        center[2] = c[3 * i + 2];
        radius = c[3 * i + 3];
        tmp += evalDensity(m_blobby, center, samplePoints[j],
                           radius); // - blurValue[j];
      }
      x[k] += (tmp - blurValue[j]) * (tmp - blurValue[j]);
    }
  }
  /* vector version, slow, minima is close to the least square version
    for(j=0; j<n; j++)
    {
      x[j] = 0.0;
      for(i=0; i<m; i++)
      {
        center[0] = c[3*i];       center[1] = c[3*i+1];   center[2] = c[3*i+2];
    radius = c[3*i+3]; x[j] += evalDensity(m_blobby, center, samplePoints[j],
    radius);// - blurValue[j];
      }
      x[j] -= blurValue[j];
    }
  */
}

void Optimizer::optimizeCenterAndRadius(int iteration) {
  // m - num of CG atoms
  // n - num of surface pts
  // p - CG radii
  // x - f(CG) - f(AA), the chi function
  double p[1], c[3];
  int m, n;
  register int i, iter;
  double start, end;
  double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
  int ret;
  opts[0] = LM_INIT_MU;
  opts[1] = 1E-15;
  opts[2] = 1E-15;
  opts[3] = 1E-20;
  opts[4] = LM_DIFF_DELTA; // relevant only if the finite difference jacobian
                           // version is used
  //     std::vector<Atom*> atomList = getAtomList();
  blurAtomList(atomList, blurValue);
  /* vector Chi function
          n = numOfSamples;
          double x[n]; // chi function
          for (i = 0; i < n; i++)
                  x[i] = 0.0;
  */
  srand((unsigned int)time(NULL));
  p[0] = adata[0].radius;
  c[0] = adata[0].center[0];
  c[1] = adata[0].center[1];
  c[2] = adata[0].center[2];
  FunctionData *fdata = new FunctionData;
  initializeFunctionData(fdata);
  start = getMyTime();
  for (iter = 0; iter < iteration; iter++) {
    // optimize radius first since the initial centers are close to the optimal
    // centers
    m = 1;
    n = m;
    double x1[n]; // chi function
    for (i = 0; i < n; i++) {
      x1[i] = 0.0;
    }
    ret = dlevmar_dif(radiusOptimizeFunction, p, x1, m, n, 1000, opts, info,
                      NULL, NULL, (void *)fdata);
    updateFunctionData(fdata, c, p[0]);
    int mm = 3 * m;
    n = mm;
    double x2[n]; // chi function
    for (i = 0; i < n; i++) {
      x2[i] = 0.0;
    }
    ret = dlevmar_dif(centerOptimizeFunction, c, x2, mm, n, 1000, opts, info,
                      NULL, NULL, (void *)fdata);
    updateFunctionData(fdata, c, p[0]);
  }
  end = getMyTime();
  adata[0].center[0] = fdata->adata[0].center[0];
  adata[0].center[1] = fdata->adata[0].center[1];
  adata[0].center[2] = fdata->adata[0].center[2];
  adata[0].radius = fabs(fdata->adata[0].radius);
  for (int i = 0; i < numOfSamples; i++) {
    delete fdata->samplePoints[i];
  }
  delete fdata->samplePoints;
  delete fdata->blurValue;
  delete fdata;
  delete blurValue;
  for (int i = 0; i < numOfSamples; i++) {
    // printf("%e\n",samplePoints[i][0]);
    delete samplePoints[i];
  }
  delete samplePoints;
  //  ret = dlevmar_dif(radiuscenterOptimizeFunction, c, x, m, n, 50, opts,
  //  info, NULL, NULL, (void*)adata);
  // no jacobian
  /*
    printf("Levenberg-Marquardt returned %d in %g iter, takes %f seconds, reason
    %g\nSolution: \n", ret, info[5], end-start, info[6]); for(i=0; i<m; ++i)
    {
        if (p[i] < 0.0)   p[i] = -p[i];
        printf("%.7g %.7g %.7g %.7g\n", c[3*i], c[3*i+1], c[3*i+2], p[i]);
    }
    printf("\nMinimization info:\n");
    for(i=0; i<LM_INFO_SZ; ++i)
      printf("%g ", info[i]);
    printf("\n");
  */
  /*
   double info[LM_INFO_SZ],
   O: information regarding the minimization. Set to NULL if don't care
   info[0]= ||e||_2 at initial p.
   info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||Dp||_2, \mu/max[J^T J]_ii ], all
   computed at estimated p. info[5]= # iterations, info[6]=reason for
   terminating: 1 - stopped by small gradient J^T e 2 - stopped by small Dp 3 -
   stopped by itmax 4 - singular matrix. Restart from current p with increased
   \mu 5 - no further error reduction is possible. Restart with increased mu 6 -
   stopped by small ||e||_2 7 - stopped by invalid (i.e. NaN or Inf) "func"
   values. This is a user error info[7]= # function evaluations info[8]= #
   Jacobian evaluations
   */
}

// Following two use  ||x-x_0|| - r
void expfunc(double *p, double *x, int m, int n, void *data) {
  register int i;
  FunctionData *fdata = (FunctionData *)data;
  for (i = 0; i < n; i++) {
    x[i] = sqrt((fdata->samplePoints[i][0] - p[0]) *
                    (fdata->samplePoints[i][0] - p[0]) +
                (fdata->samplePoints[i][1] - p[1]) *
                    (fdata->samplePoints[i][1] - p[1]) +
                (fdata->samplePoints[i][2] - p[2]) *
                    (fdata->samplePoints[i][2] - p[2])) -
           p[3];
  }
}

void jacexpfunc(double *p, double *jac, int m, int n, void *data) {
  register int i, j;
  FunctionData *fdata = (FunctionData *)data;
  for (i = j = 0; i < n; i++) {
    jac[j++] = (p[0] - fdata->samplePoints[i][0]) /
               sqrt((fdata->samplePoints[i][0] - p[0]) *
                        (fdata->samplePoints[i][0] - p[0]) +
                    (fdata->samplePoints[i][1] - p[1]) *
                        (fdata->samplePoints[i][1] - p[1]) +
                    (fdata->samplePoints[i][2] - p[2]) *
                        (fdata->samplePoints[i][2] - p[2]));
    jac[j++] = (p[1] - fdata->samplePoints[i][1]) /
               sqrt((fdata->samplePoints[i][0] - p[0]) *
                        (fdata->samplePoints[i][0] - p[0]) +
                    (fdata->samplePoints[i][1] - p[1]) *
                        (fdata->samplePoints[i][1] - p[1]) +
                    (fdata->samplePoints[i][2] - p[2]) *
                        (fdata->samplePoints[i][2] - p[2]));
    jac[j++] = (p[2] - fdata->samplePoints[i][2]) /
               sqrt((fdata->samplePoints[i][0] - p[0]) *
                        (fdata->samplePoints[i][0] - p[0]) +
                    (fdata->samplePoints[i][1] - p[1]) *
                        (fdata->samplePoints[i][1] - p[1]) +
                    (fdata->samplePoints[i][2] - p[2]) *
                        (fdata->samplePoints[i][2] - p[2]));
    jac[j++] = -1.0;
  }
}

// Following two use ||x-x_0||^2 - r^2. After optimization, radius seems a
// little bigger than the above one.
void expfuncN(double *p, double *x, int m, int n, void *data) {
  register int i;
  FunctionData *fdata = (FunctionData *)data;
  for (i = 0; i < n; i++) {
    x[i] = (fdata->samplePoints[i][0] - p[0]) *
               (fdata->samplePoints[i][0] - p[0]) +
           (fdata->samplePoints[i][1] - p[1]) *
               (fdata->samplePoints[i][1] - p[1]) +
           (fdata->samplePoints[i][2] - p[2]) *
               (fdata->samplePoints[i][2] - p[2]) -
           p[3] * p[3];
  }
}

void jacexpfuncN(double *p, double *jac, int m, int n, void *data) {
  register int i, j;
  FunctionData *fdata = (FunctionData *)data;
  for (i = j = 0; i < n; i++) {
    jac[j++] = 2 * (p[0] - fdata->samplePoints[i][0]);
    jac[j++] = 2 * (p[1] - fdata->samplePoints[i][1]);
    jac[j++] = 2 * (p[2] - fdata->samplePoints[i][2]);
    jac[j++] = -2.0 * p[3];
  }
}

void Optimizer::updateFunctionDataN(FunctionData *fdata, double *p) {
  fdata->adata[0].center[0] = p[0];
  fdata->adata[0].center[1] = p[1];
  fdata->adata[0].center[2] = p[2];
  fdata->adata[0].radius = p[3];
}

void Optimizer::optimizeCenterAndRadiusN(int iteration) {
  int n = numOfSamples, m = 4;
  double p[m], x[n], opts[LM_OPTS_SZ], info[LM_INFO_SZ];
  register int i;
  int ret;
  p[0] = adata[0].center[0];
  p[1] = adata[0].center[1];
  p[2] = adata[0].center[2];
  p[3] = adata[0].radius;
  FunctionData *fdata = new FunctionData;
  initializeFunctionData(fdata);
  for (i = 0; i < n; i++) {
    x[i] = 0;
  }
  /* optimization control parameters; passing to levmar NULL instead of opts
   * reverts to defaults */
  opts[0] = LM_INIT_MU;
  opts[1] = 1E-15;
  opts[2] = 1E-15;
  opts[3] = 1E-20;
  opts[4] = LM_DIFF_DELTA; // relevant only if the finite difference Jacobian
                           // version is used
  //     for(int it =0; it<iteration; it++){
  /* invoke the optimization function */
  ret = dlevmar_der(expfuncN, jacexpfuncN, p, x, m, n, 1000, opts, info, NULL,
                    NULL, (void *)fdata); // with analytic Jacobian
  //  ret=dlevmar_dif(expfuncN, p, x, m, n, 1000, opts, info, NULL, NULL,
  //  (void*)fdata); // without analytic Jacobian
  // seems no difference for with or without analytic Jacobian.
  cout << "Iteration: " << ret << endl;
  updateFunctionDataN(fdata, p);
  // }
  adata[0].center[0] = fdata->adata[0].center[0];
  adata[0].center[1] = fdata->adata[0].center[1];
  adata[0].center[2] = fdata->adata[0].center[2];
  adata[0].radius = fdata->adata[0].radius;
  for (int i = 0; i < numOfSamples; i++) {
    delete fdata->samplePoints[i];
  }
  delete fdata->samplePoints;
  delete fdata->blurValue;
  delete fdata;
  delete blurValue;
  for (int i = 0; i < numOfSamples; i++) {
    // printf("%e\n",samplePoints[i][0]);
    delete samplePoints[i];
  }
  delete samplePoints;
}
