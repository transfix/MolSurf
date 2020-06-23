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
#ifndef BLURRING_H
#define BLURRING_H

#include <Utility/utility.h>
#include <SmoothingFunction/Gaussian.h>
#include <MultiresolutionSurface/fastfft.h>
#include <MultiresolutionSurface/fftAlgo.h>
#include <MultiresolutionSurface/FastInverse.h>
#include <ByteOrder/ByteSwapping.h>
#include <fftw3.h>
#include <SmoothingFunction/CardinalBSpline.h>

bool getBlurMap(fftw_complex* blurmap, fftw_complex* rbfCenterValues, double* rbfCenterPositions, fftw_complex* rbf1D, int numCenters, int numFreq, int gridSize, int gridSizeBig, double alpha, int interpFuncExtent, int rbfWidth, float min[3], float max[3], SmoothingFunction* smoothingFunction, double* min2, double* max2);

// here we're breaking getBlurmap into four parts.  In order of call: getCenterFreq, getRBFFreq, multiplyWithTensor, inverseFFT
bool getCenterFreq(fftw_complex* centerFrequencies, fftw_complex* rbfCenterValues, double* rbfCenterPositions, int numCenters, int numFreq, int gridSize, double alpha, int interpFuncExtent, int rbfWidth, float min[3], float max[3], SmoothingFunction* smoothingFunctionm);
bool getRBFFreq(fftw_complex* rbf1D, int gridSize, fftw_complex* RBFFrequencies);
bool multiplyWithTensor(fftw_complex* centerFrequencies, fftw_complex* RBFFrequencies, int gridSize, int numFreq);
bool inverseFFT(fftw_complex* blurmap, fftw_complex* centerFrequencies, int numFreq, int gridSize);
bool initializeFFTW(int size);
#endif
