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
#include <MultiresolutionSurface/FastInverse.h>

extern double getTime();
extern bool fftshift3D(fftw_complex *data, int length);
extern bool compute3DFFT(fftw_complex *dest, fftw_complex *src, int dim1,
                         int dim2, int dim3, int type, unsigned int flags);
extern bool getNormalizedPositions(double *rbfCenterPositions, double *xk,
                                   double *yk, double *zk, int rbfWidth,
                                   int gridSize, int numCenters,
                                   int interpFuncExtent, int numFreq,
                                   double alpha, double *min, double *max);
extern bool multiplyWithLarger1DTensor(fftw_complex *frequencies, int numFreq,
                                       fftw_complex *largerbf, int gridSize);
extern bool computeFFT1d(fftw_complex *dest, fftw_complex *src, int dim,
                         int type, unsigned int flags);
double computeValueAt(int indexX, int indexY, int indexZ, int M, int N, int m,
                      int alpha, fftw_complex *g,
                      SmoothingFunction *smoothingFunction, double *iPhiVals,
                      double *jPhiVals, double *kPhiVals, int *iIndex,
                      int *jIndex, int *kIndex) {
  if (!g || m < 0 || M > N || M < 1 || N < 1 || indexX < 0 || indexX >= N ||
      alpha < 1) {
    return 0;
  }
  int n = M * alpha;
  double f = n / ((double)N);
  int gaussPosX = (int)(indexX * f);
  int gaussPosY = (int)(indexY * f);
  int gaussPosZ = (int)(indexZ * f);
  if (gaussPosX * N != indexX * n) {
    gaussPosX++;
  }
  if (gaussPosY * N != indexY * n) {
    gaussPosY++;
  }
  if (gaussPosZ * N != indexZ * n) {
    gaussPosZ++;
  }
  double distFromCenterX = gaussPosX - indexX * f;
  double distFromCenterY = gaussPosY - indexY * f;
  double distFromCenterZ = gaussPosZ - indexZ * f;
  int gIndexX = (gaussPosX - m) % n;
  if (gIndexX < 0) {
    gIndexX += n;
  }
  int gIndexY = (gaussPosY - m) % n;
  if (gIndexY < 0) {
    gIndexY += n;
  }
  int gIndexZ = (gaussPosZ - m) % n;
  if (gIndexZ < 0) {
    gIndexZ += n;
  }
  int i, j, k;
  double sumR = 0;
  int m2 = 2 * m + 1;
  {
    int c;
    for (c = 0; c <= m2; c++) {
      // init phi values
      double distZ = fabs(c - m + distFromCenterZ);
      iPhiVals[c] = smoothingFunction->getPhi(distZ);
      double distY = fabs(c - m + distFromCenterY);
      jPhiVals[c] = smoothingFunction->getPhi(distY);
      double distX = fabs(c - m + distFromCenterX);
      kPhiVals[c] = smoothingFunction->getPhi(distX);
      // init index values
      iIndex[c] = gIndexZ + c;
      iIndex[c] = iIndex[c] % n;
      iIndex[c] *= n * n;
      jIndex[c] = gIndexY + c;
      jIndex[c] = jIndex[c] % n;
      jIndex[c] *= n;
      kIndex[c] = gIndexX + c;
      kIndex[c] = kIndex[c] % n;
    }
  }
  for (i = 0; i <= m2; i++) {
    double iPhiVal = iPhiVals[i];
    int cur_iIndex = iIndex[i];
    for (j = 0; j <= m2; j++) {
      double ijPhiVal = iPhiVal * jPhiVals[j];
      int cur_ijIndex = cur_iIndex + jIndex[j];
      for (k = 0; k <= m2; k++) {
        double ijkPhiVal = ijPhiVal * kPhiVals[k];
        int cur_ijkIndex = cur_ijIndex + kIndex[k];
        sumR += ijkPhiVal * g[cur_ijkIndex][0];
      }
    }
  }
  return sumR;
}

bool slowInverse(fftw_complex *frequencies, int M, fftw_complex *opt, int N,
                 double alpha, int m, SmoothingFunction *smoothingFunction) {
  double bsplineGridTime = getTime();
  if (!frequencies || !opt || (M < 1) || (alpha < 1) || (m < 1) || (N < 1)) {
    return false;
  }
  int n = (int)(M * alpha);
  fftw_complex *g_hat =
      (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n * n * n);
  // three steps.
  // step 1. divide by phi
  {
    int i;
    for (i = 0; i < n * n * n; i++) {
      g_hat[i][0] = g_hat[i][1] = 0;
    }
  }
  {
    int c = 0;
    int i, j, k;
    double alphaM = (double)n;
    for (i = 0; i < M; i++) {
      double iPsiBar = smoothingFunction->getPhiHat(abs(i - M / 2));
      int newI = i - M / 2;
      if (newI < 0) {
        newI += n;
      }
      newI *= n * n;
      for (j = 0; j < M; j++) {
        double ijPsiBar =
            iPsiBar * smoothingFunction->getPhiHat(abs(j - M / 2));
        int newJ = j - M / 2;
        if (newJ < 0) {
          newJ += n;
        }
        newJ *= n;
        newJ += newI;
        for (k = 0; k < M; k++) {
          double psiBar =
              ijPsiBar * smoothingFunction->getPhiHat(abs(k - M / 2));
          int newK = k - M / 2;
          if (newK < 0) {
            newK += n;
          }
          newK += newJ;
          g_hat[newK][0] = frequencies[c][0] / psiBar;
          g_hat[newK][1] = frequencies[c][1] / psiBar;
          c++;
        }
      }
    }
  }
  // step 2. perform inverse FFT
  fftw_complex *g =
      (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n * n * n);
  {
    fftw_plan p;
    p = fftw_plan_dft_3d(n, n, n, g_hat, g, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }
  fftw_free(g_hat);
  printf("Inverse grid construction took %lf\n", getTime() - bsplineGridTime);
  // step 3. sum surrounding splines at each point
  {
    double *iPhiVals = new double[2 * m + 1];
    double *jPhiVals = new double[2 * m + 1];
    double *kPhiVals = new double[2 * m + 1];
    int *iIndex = new int[2 * m + 1];
    int *jIndex = new int[2 * m + 1];
    int *kIndex = new int[2 * m + 1];
    int i, j, k;
    int c = 0;
    for (i = 0; i < N; i++) {
      double lastStep = getTime();
      for (j = 0; j < N; j++) {
        for (k = 0; k < N; k++) {
          opt[c][0] = computeValueAt(k, j, i, M, N, m, (int)alpha, g,
                                     smoothingFunction, iPhiVals, jPhiVals,
                                     kPhiVals, iIndex, jIndex, kIndex);
          opt[c][1] = 0;
          c++;
        }
      }
      printf("Finished plane number %d of %d in %lf seconds\r", i, N,
             getTime() - lastStep);
    }
  }
  printf("\n");
  fftw_free(g);
  return true;
}

bool getBlurMapWithFastInverse(fftw_complex *blurmap,
                               fftw_complex *rbfCenterValues,
                               double *rbfCenterPositions, fftw_complex *rbf1D,
                               int numCenters, int numFreq, int gridSize,
                               double alpha, int interpFuncExtent, int rbfWidth,
                               SmoothingFunction *smoothingFunction,
                               double *min, double *max) {
  if (!blurmap || !rbfCenterValues || !rbfCenterPositions || !rbf1D ||
      (numCenters < 1) || (numFreq < 1) || (gridSize < 1) || (alpha < 1) ||
      (interpFuncExtent < 1) || !smoothingFunction) {
    return false;
  }
  int intAlpha = (int)alpha;
  fftw_complex *centerFrequencies = (fftw_complex *)fftw_malloc(
      sizeof(fftw_complex) * numFreq * numFreq * numFreq);
  double forwardTime = 0;
  double ndftTime;
  // get the frequencies of the centers.
  {
    double *xk = new double[numCenters];
    double *yk = new double[numCenters];
    double *zk = new double[numCenters];
    if (!getNormalizedPositions(rbfCenterPositions, xk, yk, zk, rbfWidth,
                                gridSize, numCenters, interpFuncExtent, numFreq,
                                alpha, min, max)) {
      return false;
    }
    double freqTime = getTime();
    forwardTime = getTime();
    ndftTime = getTime();
    if (!getCenterFrequencies(numCenters, alpha, interpFuncExtent, numFreq,
                              gridSize, xk, yk, zk, rbfCenterValues,
                              centerFrequencies, 0, smoothingFunction)) {
      return false;
    }
    if (!fftshift3D(centerFrequencies, numFreq)) {
      return false;
    }
    delete[] xk;
    delete[] yk;
    delete[] zk;
    printf("The time taken to do nfft was %f\n", getTime() - freqTime);
  }
  fftw_complex *RBFFrequencies =
      (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * gridSize);
  if (!computeFFT1d(RBFFrequencies, rbf1D, gridSize, FFTW_FORWARD,
                    FFTW_ESTIMATE)) {
    return false;
  }
  // perform multiplication with 1D tensor
  if (!multiplyWithLarger1DTensor(centerFrequencies, numFreq, RBFFrequencies,
                                  gridSize)) {
    return false;
  }
  fftw_free(RBFFrequencies);
  printf("Forward step took %lf \n", getTime() - forwardTime);
  if (!fftshift3D(centerFrequencies, numFreq)) {
    return false;
  }
  if (!slowInverse(centerFrequencies, numFreq, blurmap, gridSize, alpha,
                   interpFuncExtent, smoothingFunction)) {
    printf("Could not perform slow inverse\n");
  }
  fftw_free(centerFrequencies);
  printf("The Fast Inverse method took %lf seconds\n", getTime() - ndftTime);
  return true;
}
