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
#include <fftw3.h>
#include <SmoothingFunction/CardinalBSpline.h>
#include <SmoothingFunction/SmoothingFunction.h>
#include <MultiresolutionSurface/blurring.h>
#include <AdaptiveIsocontouring/comp.h>
#include <fftw3.h>

void getGaussian(double b, int grid_size, double error, double width, double r,
                 fftw_complex *rbf, int &rbfWidth);

// a wrapper class for nfft.
class NfftSummer {
public:
  NfftSummer(bool holder_trick = false);
  ~NfftSummer();
  void performSumWithNewCenters(vector<float *> in, float *sum_out,
                                float *new_min, float *new_max);
  void performSumWithNewDecay(float decay, float *sum_out);
  void recompRBFFreq(float decay);
  void recompCenterFreq(fftw_complex *centerVal, float *center_pos,
                        int numCenters);
  void initialize(int gridSize, int alpha, int numFreq, int interpFuncExtent,
                  float input_width, float radius, float error, float decay,
                  float *_min, float *_max);
  void isOuter(bool set) { is_outer = set; }

private:
  void combineStep();
  int gridSize, alpha, numFreq, rbfWidth, interpFuncExtent;
  int rgSize;
  float input_width;
  float radius, error;
  float decay;
  float min[3];
  float max[3];
  bool holder_trick; // holder_trick is used to distinguish between the inner
                     // sum and the outer sum
  bool is_outer;
  SmoothingFunction *smoothingFunction;
  fftw_complex *RBFFreq;
  fftw_complex *centerVals;
  fftw_complex *centerFreq;
  fftw_complex *blurmap;
  fftw_complex *rbf_dud;
  fftw_complex *RBF1D;
  fftw_complex *centerFreqTmp;
  double *center_pos_d;
};
