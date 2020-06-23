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
#include <MultiresolutionSurface/fftAlgo.h>

extern double getTime();
extern bool computeFFT1d(fftw_complex* dest,fftw_complex* src, int dim, int type, unsigned int flags);
extern bool compute3DFFT(fftw_complex* dest,fftw_complex* src, int dim1, int dim2, int dim3, int type, unsigned int flags);
extern bool multiplyWith1DTensor(fftw_complex* dest, fftw_complex* input1, fftw_complex* input2, int length, int length2);
extern bool getNormalizedPositions(double* rbfCenterPositions, double* xk, double* yk, double* zk, int rbfWidth, int gridSize, int numCenters, int interpFuncExtent, int numFreq, double alpha, double* min, double* max);
bool getBlurMapWithFullFFT(fftw_complex* blurmap, fftw_complex* rbfCenterValues, double* rbfCenterPositions, fftw_complex* rbf1D, int numCenters, int numFreq, int gridSizeSmall, int gridSize, double alpha, int interpFuncExtent, int rbfWidth, double* min, double* max)
{
	if(!blurmap || !rbfCenterValues || !rbfCenterPositions || !rbf1D || (numCenters<1) || (numFreq<1) || (gridSize<1) || (alpha<1) || (interpFuncExtent<1))
	{
		return false;
	}
	fftw_complex* centerFrequencies = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * gridSize*gridSize*gridSize);
	{
		int i;
		for(i=0; i<gridSize*gridSize*gridSize; i++)
		{
			centerFrequencies[i][0] = centerFrequencies[i][1] = 0;
		}
	}
	double fftAlgoTime;
	// get the frequencies of the centers.
	{
		printf("Obtaining the normalized positions\n");
		// shift values to [-0.5 .. 0.5)
		// There needs to be a gap to allow for convolution
		//    If small extent rbf functions, then we only need to 0 -pad with that size,
		//    else, we may need to 0-pad with cuberoot(3)* N.
		double* xk = new double[numCenters];
		double* yk = new double[numCenters];
		double* zk = new double[numCenters];
		if(!getNormalizedPositions(rbfCenterPositions, xk, yk, zk, rbfWidth, gridSizeSmall, numCenters, interpFuncExtent, numFreq, alpha, min, max))
		{
			printf("Could not obtain the normalized positions\n");
			delete [] xk;
			delete [] yk;
			delete [] zk;
			return false;
		}
		fftw_complex* spaceValsGrid = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * gridSize*gridSize*gridSize);
		// fill up the grid to perform full FFT
		{
			int i;
			for(i=0; i<gridSize*gridSize*gridSize; i++)
			{
				spaceValsGrid[i][0] = spaceValsGrid[i][1] = 0;
			}
			int minX=1000;
			int maxX=-1000;
			int minY=1000;
			int maxY=-1000;
			int minZ=1000;
			int maxZ=-1000;
			for(i=0; i<numCenters; i++)
			{
				int xpos = (int)((xk[i]+0.5) * gridSize);
				int ypos = (int)((yk[i]+0.5) * gridSize);
				int zpos = (int)((zk[i]+0.5) * gridSize);
				if(minX > xpos)
				{
					minX = xpos;
				}
				if(minY > ypos)
				{
					minY = ypos;
				}
				if(minZ > zpos)
				{
					minZ = zpos;
				}
				if(maxX < xpos)
				{
					maxX = xpos;
				}
				if(maxY < ypos)
				{
					maxY = ypos;
				}
				if(maxZ < zpos)
				{
					maxZ = zpos;
				}
				int pos = zpos*gridSize*gridSize + ypos*gridSize + xpos;
				spaceValsGrid[pos][0] = rbfCenterValues[i][0];
				spaceValsGrid[pos][1] = rbfCenterValues[i][1];
			}
			printf("The min max pos was [%d %d %d] - [%d %d %d]\n", minX, minY, minZ, maxX, maxY, maxZ);
		}
		fftAlgoTime = getTime();
		{
			double freqTime = getTime();
			if(!compute3DFFT(centerFrequencies, spaceValsGrid, gridSize, gridSize, gridSize, FFTW_FORWARD, FFTW_ESTIMATE))
			{
				printf("Could not perform the full FFT\n");
				delete [] xk;
				delete [] yk;
				delete [] zk;
				fftw_free(centerFrequencies);
				fftw_free(spaceValsGrid);
				return false;
			}
			printf("The time taken to do full FFT  was %f\n", getTime()-freqTime);
		}
		delete [] xk;
		delete [] yk;
		delete [] zk;
		fftw_free(spaceValsGrid);
	}
	fftw_complex* RBFFrequencies = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * gridSize);
	if(!computeFFT1d(RBFFrequencies, rbf1D, gridSize, FFTW_FORWARD, FFTW_ESTIMATE))
	{
		printf("Could not compute the RBF frequencies\n");
		fftw_free(centerFrequencies);
		fftw_free(RBFFrequencies);
		return false;
	}
	fftw_complex* frequencies = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * gridSize*gridSize*gridSize);
	if(!multiplyWith1DTensor(frequencies, centerFrequencies, RBFFrequencies, gridSize, gridSize))
	{
		printf("Could not multiply the two frequency maps\n");
		fftw_free(frequencies);
		fftw_free(centerFrequencies);
		fftw_free(RBFFrequencies);
		return 0;
	}
	fftw_free(centerFrequencies);
	fftw_free(RBFFrequencies);
	if(!compute3DFFT(blurmap, frequencies, gridSize, gridSize, gridSize, FFTW_BACKWARD, FFTW_ESTIMATE))
	{
		printf("Could not compute the blurmap\n");
		return 0;
	}
	fftw_free(frequencies);
	printf("The Full FFT algo took %lf seconds\n", getTime() -fftAlgoTime);
	return true;
}
