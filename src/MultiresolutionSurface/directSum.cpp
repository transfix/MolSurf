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

extern double getTime();
extern bool computeFFT1d(fftw_complex* dest,fftw_complex* src, int dim, int type, unsigned int flags);
extern bool computeFFT3d(fftw_complex* dest,fftw_complex* src, int dim1, int dim2, int dim3, int type, unsigned int flags);
extern bool multiplyWith1DTensor(fftw_complex* dest, fftw_complex* input1, fftw_complex* input2, int length, int length2);
extern bool getNormalizedPositions(double* rbfCenterPositions, double* xk, double* yk, double* zk, int rbfWidth, int gridSize, int numCenters, int interpFuncExtent, int numFreq, double alpha, double* min, double* max);
bool getBlurMapWithSlowAlgo(fftw_complex* blurmap, fftw_complex* rbfCenterValues, double* rbfCenterPositions, fftw_complex* rbf1D, int numCenters, int numFreq, int gridSize, double alpha, int interpFuncExtent, int rbfWidth, double* min, double* max)
{
	if(!blurmap || !rbfCenterValues || !rbfCenterPositions || !rbf1D || (numCenters<1) || (numFreq<1) || (gridSize<1) || (alpha<1) || (interpFuncExtent<1))
	{
		return false;
	}
	printf("Obtaining the normalized positions\n");
	// shift values to [-0.5 .. 0.5)
	// There needs to be a gap to allow for convolution
	//    If small extent rbf functions, then we only need to 0 -pad with that size,
	//    else, we may need to 0-pad with cuberoot(3)* N.
	double* xk = new double[numCenters];
	double* yk = new double[numCenters];
	double* zk = new double[numCenters];
	if(!getNormalizedPositions(rbfCenterPositions, xk, yk, zk, rbfWidth, gridSize, numCenters, interpFuncExtent, numFreq, alpha, min, max))
	{
		printf("Could not obtain the normalized positions\n");
		delete [] xk;
		delete [] yk;
		delete [] zk;
		return false;
	}
	int i;
	for(i=0; i<gridSize*gridSize*gridSize; i++)
	{
		blurmap[i][0] = blurmap[i][1] = 0;
	}
	double slowAlgoTime = getTime();
	for(i=0; i<numCenters; i++)
	{
		int xpos = (int)((xk[i]+0.5) * gridSize);
		int ypos = (int)((yk[i]+0.5) * gridSize);
		int zpos = (int)((zk[i]+0.5) * gridSize);
		int xStart = xpos-rbfWidth+1;
		int yStart = ypos-rbfWidth+1;
		int zStart = zpos-rbfWidth+1;
		int xEnd = xpos+rbfWidth-1;
		int yEnd = ypos+rbfWidth-1;
		int zEnd = zpos+rbfWidth-1;
		if(xStart < 0)
		{
			xStart = 0;
		}
		if(yStart < 0)
		{
			yStart = 0;
		}
		if(zStart < 0)
		{
			zStart = 0;
		}
		if(xEnd >= gridSize)
		{
			xEnd = gridSize-1;
		}
		if(yEnd >= gridSize)
		{
			yEnd = gridSize-1;
		}
		if(zEnd >= gridSize)
		{
			zEnd = gridSize-1;
		}
		int c1, c2, c3;
		for(c1=zStart; c1<=zEnd; c1++)
		{
			double zGaussian = rbf1D[abs(zpos-c1)][0];
			for(c2=yStart; c2<=yEnd; c2++)
			{
				double zyGaussian = zGaussian*rbf1D[abs(ypos-c2)][0];
				for(c3=xStart; c3<=xEnd; c3++)
				{
					blurmap[c1* gridSize* gridSize + c2* gridSize + c3][0] += zyGaussian*rbf1D[abs(xpos-c3)][0];
				}
			}
		}
	}
	printf("The slow algo takes %lf seconds\n", getTime() - slowAlgoTime);
	return true;
}
