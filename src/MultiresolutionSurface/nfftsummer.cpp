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
#include <MultiresolutionSurface/nfftsummer.h>

NfftSummer::NfftSummer(bool holder_trick) {
  RBFFreq = NULL;
  centerVals = NULL;
  centerFreq = NULL;
  blurmap = NULL;
  smoothingFunction = NULL;
  RBF1D = NULL;
  rbf_dud = NULL;
  centerFreqTmp = NULL;
  center_pos_d = NULL;
  this->holder_trick = holder_trick; // to size the atoms appropriately
}

void NfftSummer::initialize(int gridSize, int alpha, int numFreq,
                            int interpFuncExtent, float input_width,
                            float radius, float error, float decay, float *_min,
                            float *_max) {
  cout << "NfftwSummer initializing" << endl;
  this->gridSize = gridSize;
  this->alpha = alpha;
  this->numFreq = numFreq;
  this->interpFuncExtent = interpFuncExtent;
  this->radius = radius;
  this->error = error;
  this->decay = decay;
  this->input_width = input_width;
  copyVect(_min, min, 3);
  copyVect(_max, max, 3);
  if (holder_trick) {
    rbfWidth = 10;
  } else {
    rbfWidth = 10;
  }
  {
    int n = alpha * numFreq;
    int y = n * gridSize / ((float)(n + 2.0 * interpFuncExtent));
    int x = (int)(((float)interpFuncExtent) * ((float)y) / ((float)n));
    rgSize = 2 * rbfWidth + 2 * x + y;
  }
  int dud;
  if (RBFFreq != NULL) {
    fftw_free(RBFFreq);
  }
  RBFFreq = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * rgSize);
  if (rbf_dud != NULL) {
    fftw_free(rbf_dud);
  }
  rbf_dud = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * rgSize);
  if (RBF1D != NULL) {
    fftw_free(RBF1D);
  }
  RBF1D = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * rgSize);
  if (centerFreqTmp != NULL) {
    fftw_free(centerFreqTmp);
  }
  centerFreqTmp = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * numFreq *
                                              numFreq * numFreq);
  recompRBFFreq(decay);
  if (blurmap != NULL) {
    fftw_free(blurmap);
  }
  if (centerFreq != NULL) {
    fftw_free(centerFreq);
  }
  if (smoothingFunction != NULL) {
    delete smoothingFunction;
  }
  smoothingFunction = new CardinalBSpline(alpha, interpFuncExtent,
                                          (int)(alpha * numFreq), rgSize);
  blurmap = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * rgSize * rgSize *
                                        rgSize);
  for (int j = 0; j < rgSize * rgSize * rgSize; j++) {
    blurmap[j][0] = blurmap[j][1] = 0;
  }
  centerFreq = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * numFreq *
                                           numFreq * numFreq);
}

NfftSummer::~NfftSummer() {
  fftw_free(blurmap);
  fftw_free(RBFFreq);
  fftw_free(centerVals);
  fftw_free(centerFreq);
  delete smoothingFunction;
}

void NfftSummer::performSumWithNewCenters(vector<float *> in, float *sum_out,
                                          float *new_min, float *new_max) {
  int num_atoms = in.size();
  float *centers = new float[num_atoms * 3];
  fftw_complex *centerVals = new fftw_complex[num_atoms];
  copyVect(new_min, min, 3);
  copyVect(new_max, max, 3);
  int j;
  for (j = 0; j < num_atoms; j++) {
    centers[j * 3] = in[j][0];
    centers[j * 3 + 1] = in[j][1];
    centers[j * 3 + 2] = in[j][2];
    if (holder_trick && j > num_atoms - 3 && false) {
      centerVals[j][0] = 0.0;
      centerVals[j][1] = 0.0;
    } else {
      centerVals[j][0] = 1.0;
      centerVals[j][1] = 0.0;
    }
  }
  recompCenterFreq(centerVals, centers, num_atoms);
  combineStep();
  float scale = rgSize * rgSize * rgSize;
  int holder = rbfWidth;
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      for (int z = 0; z < gridSize; z++) {
        int ind = (x + holder) + (y + holder) * rgSize +
                  (z + holder) * rgSize * rgSize;
        int ind_j = x + y * gridSize + z * gridSize * gridSize;
        sum_out[ind_j] = (float)(blurmap[ind][0] / scale);
      }
    }
  }
  delete[] centers;
  delete[] centerVals;
}

void NfftSummer::performSumWithNewDecay(float decay, float *sum_out) {
  recompRBFFreq(decay);
  combineStep();
  float scale = rgSize * rgSize * rgSize;
  int holder = rbfWidth;
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      for (int z = 0; z < gridSize; z++) {
        int ind = (x + holder) + (y + holder) * rgSize +
                  (z + holder) * rgSize * rgSize;
        int ind_j = x + y * gridSize + z * gridSize * gridSize;
        sum_out[ind_j] = (float)(blurmap[ind][0] / scale);
      }
    }
  }
}

void NfftSummer::recompRBFFreq(float _decay) {
  this->decay = _decay;
  int dud;
  float new_width = input_width * rgSize / (float)gridSize;
  getGaussian(-_decay, rgSize, error, new_width, radius, rbf_dud, dud);
  int j;
  for (j = 0; j < rgSize; j++) {
    RBF1D[j][1] = 0;
    if (j < rbfWidth) {
      RBF1D[j][0] = rbf_dud[j][0];
      if (j < rbfWidth - 1) {
        RBF1D[rgSize - 1 - j][0] = rbf_dud[j + 1][0];
      }
    } else if (j >= rbfWidth && j < rgSize - rbfWidth - 1) {
      RBF1D[j][0] = 0;
    }
  }
  for (j = 0; j < rgSize; j++) {
    RBFFreq[j][0] = RBFFreq[j][1] = 0;
  }
  getRBFFreq(RBF1D, rgSize, RBFFreq);
}

void NfftSummer::recompCenterFreq(fftw_complex *centerVal, float *center_pos,
                                  int numCenters) {
  float dud[3];
  int j;
  for (j = 0; j < numFreq * numFreq * numFreq; j++) {
    centerFreq[j][0] = centerFreq[j][1] = 0;
  }
  center_pos_d =
      new double[numCenters *
                 3]; // this is very dangerous, is there a way around it?
  for (j = 0; j < numCenters * 3; j++) {
    center_pos_d[j] = (double)center_pos[j];
  }
  float fmin[3] = {(float)min[0], (float)min[1], (float)min[2]};
  float fmax[3] = {(float)max[0], (float)max[1], (float)max[2]};
  bool returned = getCenterFreq(
      centerFreq, centerVal, center_pos_d, numCenters, numFreq, gridSize,
      (double)alpha, interpFuncExtent, rbfWidth, fmin, fmax, smoothingFunction);
  if (!returned) {
    cout << "getCenterFreq not returned!!" << endl;
  }
  delete[] center_pos_d;
  center_pos_d = NULL;
}

void NfftSummer::combineStep() {
  int j;
  for (j = 0; j < numFreq * numFreq * numFreq; j++) {
    centerFreqTmp[j][0] = centerFreq[j][0];
    centerFreqTmp[j][1] = centerFreq[j][1];
  }
  for (j = 0; j < rgSize * rgSize * rgSize; j++) {
    blurmap[j][0] = blurmap[j][1] = 0;
  }
  multiplyWithTensor(centerFreqTmp, RBFFreq, rgSize, numFreq);
  inverseFFT(blurmap, centerFreqTmp, numFreq, rgSize);
}
