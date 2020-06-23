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
#include <MultiresolutionSurface/blurring.h>

extern bool getBlurMapWithSlowAlgo(fftw_complex* blurmap, fftw_complex* rbfCenterValues, double* rbfCenterPositions, fftw_complex* rbf1D, int numCenters, int numFreq, int gridSize, double alpha, int interpFuncExtent, int rbfWidth);


// Get the current time in seconds as a double value
double getTime2()
{
#ifdef _WIN32
	time_t ltime;
	_timeb tstruct;
	time(&ltime);
	_ftime(&tstruct);
	return (double)(ltime + 1e-3*(tstruct.millitm));
#else
	struct timeval t;
	gettimeofday(&t, NULL);
	return (double)(t.tv_sec + 1e-6*t.tv_usec);
#endif
}


/* Computes FFT using the FFTW3 code. Please call fftw_whatever_they_want_to_do with the */
/* appropriate sizes before calling this. Call this with preferably FFTW_MEASURE flag. */
bool computeFFT1d(fftw_complex* dest,fftw_complex* src, int dim, int type, unsigned int flags)
{
	if(!dest || !src || (dim<1))
	{
		if(!dest)
		{
			cout<<"no dest"<<endl;
		}
		if(!src)
		{
			cout<<"no src"<<endl;
		}
		if(dim<1)
		{
			cout<<"dim less than 1"<<endl;
		}
		return false;
	}
	cout<<"dim: "<<dim<<endl;
	cout<<"flag: "<<flags<<endl;
	cout<<"type: "<<type<<endl;
	fftw_plan p;
	p = fftw_plan_dft_1d(dim, src, dest, type, flags);
	fftw_execute(p);
	fftw_destroy_plan(p);
	return true;
}

bool computeFFT2d(fftw_complex* dest,fftw_complex* src, int dim1, int dim2, int type, unsigned int flags)
{
	if(!dest || !src || dim1<1 || dim2<1)
	{
		return false;
	}
	fftw_plan p;
	p = fftw_plan_dft_2d(dim1, dim2, src, dest, type, flags);
	fftw_execute(p);
	fftw_destroy_plan(p);
	return true;
}

bool compute3DFFT(fftw_complex* dest,fftw_complex* src, int dim1, int dim2, int dim3, int type, unsigned int flags)
{
	if(!dest || !src || dim1<1 || dim2<1 || dim3<1)
	{
		return false;
	}
	fftw_plan p;
	p = fftw_plan_dft_3d(dim1, dim2, dim3, src, dest, type, flags);
	fftw_execute(p);
	fftw_destroy_plan(p);
	return true;
}

/* We should call this function for all sizes of ffts that we want to compute later */
/* It is a good idea to call this somewhere in the main function as soon as we know what */
/* sizes are relevant. */
bool initializeFFTW(int size)
{
	if(size<1)
	{
		return false;
	}
	fftw_complex* data = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*size*size*size);
	{
		int i;
		for(i=0; i<size*size*size; i++)
		{
			data[i][0] = (double)rand()/((double)(RAND_MAX));
			data[i][1] = (double)rand()/((double)(RAND_MAX));
		}
	}
	if(!compute3DFFT(data,data, size, size, size, FFTW_BACKWARD, FFTW_MEASURE))
	{
		return false;
	}
	if(!compute3DFFT(data,data, size, size, size, FFTW_FORWARD, FFTW_MEASURE))
	{
		return false;
	}
	fftw_free(data);
	return true;
}

/* Do some sanity check on the inputs. There are many many small things to take care of. */
int checkInputs(int lowResGridSize, int interpFuncExtent, int numFreq, int gridSize, char* file1, char* file2)
{
	int v = (int)(ceil(pow(numFreq, 1.0/3.0)));
	int vby2 = v/2;
	if(v* v* v != numFreq)
	{
		printf("numFreq should be a cube of an even integer\n");
		return 0;
	}
	if(v < 2)
	{
		printf("numFreq should be a cube of an even integer\n");
		return 0;
	}
	if(vby2+vby2 != v)
	{
		printf("numFreq should be a cube of an even integer\n");
		return 0;
	}
	if(interpFuncExtent < 1)
	{
		printf("interpFuncExtent is too small\n");
		return 0;
	}
	// check for M to be a cube and its cube root to be an even number!
	return 1;
}

bool fftshift3D(fftw_complex* data, int length)
{
	if(!data || (length<1))
	{
		return false;
	}
	fftw_complex* temp = (fftw_complex*)malloc(sizeof(fftw_complex) * length * length * length);
	int c;
	for(c=0; c< length * length * length; c++)
	{
		temp[c][0] = data[c][0];
		temp[c][1] = data[c][1];
	}
	int i, j, k;
	for(i=0; i<length; i++)
	{
		for(j=0; j<length; j++)
		{
			for(k=0; k<length; k++)
			{
				int shift_i = (i+length/2)%length;
				int shift_j = (j+length/2)%length;
				int shift_k = (k+length/2)%length;
				data[shift_i* length* length + shift_j* length + shift_k][0] = temp[i*length*length + j*length + k][0];
				data[shift_i* length* length + shift_j* length + shift_k][1] = temp[i*length*length + j*length + k][1];
			}
		}
	}
	fftw_free(temp);
	return true;
}

bool multiplyWithLarger1DTensor(fftw_complex* frequencies, int numFreq, fftw_complex* largerbf, int gridSize)
{
	if(!frequencies || !largerbf || (numFreq<1) || (gridSize<1))
	{
		return false;
	}
	int c = 0;
	int i, j, k;
	fftw_complex* rbf = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*numFreq);
	{
		int i;
		for(i=0; i<numFreq; i++)
		{
			if(i<numFreq/2)
			{
				rbf[i][0] = largerbf[i][0];
				rbf[i][1] = largerbf[i][1];
			}
			else
			{
				rbf[i][0] = largerbf[gridSize-numFreq+i][0];
				rbf[i][1] = largerbf[gridSize-numFreq+i][1];
			}
		}
	}
	for(j=0; j<numFreq; j++)
	{
		for(k=0; k<numFreq; k++)
		{
			for(i=0; i<numFreq; i++)
			{
				double val_re = frequencies[c][0];
				double val_im = frequencies[c][1];
				double t1_re = val_re * rbf[i][0] - val_im * rbf[i][1];
				double t1_im = val_re * rbf[i][1] + val_im * rbf[i][0];
				double t2_re = t1_re * rbf[j][0] - t1_im * rbf[j][1];
				double t2_im = t1_re * rbf[j][1] + t1_im * rbf[j][0];
				double t3_re = t2_re * rbf[k][0] - t2_im * rbf[k][1];
				double t3_im = t2_re * rbf[k][1] + t2_im * rbf[k][0];
				frequencies[c][0] = t3_re;
				frequencies[c][1] = t3_im;
				c++;
			}
		}
	}
	fftw_free(rbf);
	return true;
}

// careful with the shifts
bool multiplyWith1DTensor(fftw_complex* dest, fftw_complex* input1, fftw_complex* input2, int length, int length2)
{
	if(!dest || !input1 || !input2 || (length<1) || (length2<1))
	{
		return false;
	}
	fftw_complex* shiftedRBFFrequencies = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*length);
	{
		int i;
		for(i=0; i<length; i++)
		{
			shiftedRBFFrequencies[i][0] = shiftedRBFFrequencies[i][1] = 0;
		}
		for(i=0; i<length/2; i++)
		{
			shiftedRBFFrequencies[i][0] = input2[i][0];
			shiftedRBFFrequencies[i][1] = input2[i][1];
		}
		for(i=length/2; i<length; i++)
		{
			shiftedRBFFrequencies[i][0] = input2[i+length2-length][0];
			shiftedRBFFrequencies[i][1] = input2[i+length2-length][1];
		}
	}
	{
		int i;
		for(i=0; i<length; i++)
		{
			shiftedRBFFrequencies[i][0] = input2[i][0];
			shiftedRBFFrequencies[i][1] = input2[i][1];
		}
	}
	int c = 0;
	int i, j, k;
	double sum_re = 0, sum_im = 0;
	for(j=0; j<length; j++)
	{
		for(k=0; k<length; k++)
		{
			for(i=0; i<length; i++)
			{
				double val_re = input1[c][0];
				double val_im = input1[c][1];
				sum_re += val_re;
				sum_im += val_im;
				double t1_re = val_re * shiftedRBFFrequencies[i][0] - val_im * shiftedRBFFrequencies[i][1];
				double t1_im = val_re * shiftedRBFFrequencies[i][1] + val_im * shiftedRBFFrequencies[i][0];
				double t2_re = t1_re * shiftedRBFFrequencies[j][0] - t1_im * shiftedRBFFrequencies[j][1];
				double t2_im = t1_re * shiftedRBFFrequencies[j][1] + t1_im * shiftedRBFFrequencies[j][0];
				double t3_re = t2_re * shiftedRBFFrequencies[k][0] - t2_im * shiftedRBFFrequencies[k][1];
				double t3_im = t2_re * shiftedRBFFrequencies[k][1] + t2_im * shiftedRBFFrequencies[k][0];
				dest[c][0] = t3_re;
				dest[c][1] = t3_im;
				c++;
			}
		}
	}
	fftw_free(shiftedRBFFrequencies);
	return true;
}

bool zeroPad(fftw_complex* dest, int destSize, fftw_complex* input, int inputSize)
{
	if(!dest || (destSize<1) || !input || (inputSize<1))
	{
		return false;
	}
	int i,j, k;
	int c = 0;
	for(i=0; i<destSize*destSize*destSize; i++)
	{
		dest[i][0] = dest[i][1] = 0;
	}
	for(i=0; i<inputSize; i++)
	{
		for(j=0; j<inputSize; j++)
		{
			for(k=0; k<inputSize; k++)
			{
				int i_index = i;
				int j_index = j;
				int k_index = k;
				if(i>=inputSize/2)
				{
					i_index = i+destSize-inputSize;
				}
				if(j>=inputSize/2)
				{
					j_index = j+destSize-inputSize;
				}
				if(k>=inputSize/2)
				{
					k_index = k+destSize-inputSize;
				}
				dest[i_index* destSize* destSize + j_index* destSize + k_index][0] = input[i*inputSize*inputSize + j*inputSize + k][0];
				dest[i_index* destSize* destSize + j_index* destSize + k_index][1] = input[i*inputSize*inputSize + j*inputSize + k][1];
			}
		}
	}
	return true;
}
bool getNormalizedPositions(double* rbfCenterPositions, double* xk, double* yk, double* zk, int rbfWidth, int gridSize, int numCenters, int interpFuncExtent, int numFreq,
							double alpha, double* min, double* max)
{
	if(!rbfCenterPositions || !xk || !yk || !zk || (rbfWidth<1) || (gridSize<1) || (numCenters<1))
	{
		return false;
	}
	if(rbfWidth > 2*gridSize)
	{
		rbfWidth = gridSize; // we dont care outside this interval
	}
	// get the largest side
	double minX, maxX, minY, maxY, minZ, maxZ;
	minX = maxX = rbfCenterPositions[0*0+0];
	minY = maxY = rbfCenterPositions[0*0+1];
	minZ = maxZ = rbfCenterPositions[0*0+2];
	minX = min[0];
	minY = min[1];
	minZ = min[2];
	maxX = max[0];
	maxY = max[1];
	maxZ = max[2];
	int i;
	double largestWidth = maxX - minX;
	if(maxY-minY > largestWidth)
	{
		largestWidth = maxY-minY;
	}
	if(maxZ-minZ > largestWidth)
	{
		largestWidth = maxZ-minZ;
	}
	double scale =1;
	double scale_all[3] = {1,1,1};
	{
		double n = alpha * numFreq;
		double m = interpFuncExtent;
		double y = (n * gridSize) / (n + 2.0 *m);
		double x = m * y / n;
		scale = y / (y + 2 * x + 2 * rbfWidth) / largestWidth;
		scale_all[0] = scale * largestWidth / (maxX - minX);
		scale_all[1] = scale * largestWidth / (maxY - minY);
		scale_all[2] = scale * largestWidth / (maxZ - minZ);
	}
	{
		double xTrans = - minX - (maxX-minX)/2.0;
		double yTrans = - minY - (maxY-minY)/2.0;
		double zTrans = - minZ - (maxZ-minZ)/2.0;
		double scaleFactor = scale;
		for(i=0; i<numCenters; i++)
		{
			xk[i] = rbfCenterPositions[i*3+0] + xTrans; // from -w/2 to +w/2
			xk[i] = xk[i] * scale_all[0]; // within -0.5 to 0.5, including appropriate 0 padding
			yk[i] = rbfCenterPositions[i*3+1] + yTrans; // from -w/2 to +w/2
			yk[i] = yk[i] * scale_all[1]; // within -0.5 to 0.5, including appropriate 0 padding
			zk[i] = rbfCenterPositions[i*3+2] + zTrans; // from -w/2 to +w/2
			zk[i] = zk[i] * scale_all[2]; // within -0.5 to 0.5, including appropriate 0 padding
		}
	}
	{
		double scale = 1;
		for(i=0; i<numCenters; i++)
		{
			xk[i] = xk[i] * scale; // within -0.5 to 0.5, including appropriate 0 padding
			yk[i] = yk[i] * scale; // within -0.5 to 0.5, including appropriate 0 padding
			zk[i] = zk[i] * scale; // within -0.5 to 0.5, including appropriate 0 padding
		}
	}
	{
		int i;
		for(i=0; i<numCenters; i++)
		{
			double t = xk[i] + 0.5;
			t = t*gridSize;
			t = (int)t;
			t = t/(double)gridSize;
			xk[i] = t-0.5;
			t = yk[i] + 0.5;
			t = t*gridSize;
			t = (int)t;
			t = t/(double)gridSize;
			yk[i] = t-0.5;
			t = zk[i] + 0.5;
			t = t*gridSize;
			t = (int)t;
			t = t/(double)gridSize;
			zk[i] = t-0.5;
		}
	}
	return true;
}

// here we're breaking getBlurmap into four parts. In order of call: getCenterFreq, getRBFFreq, multiplyWithTensor, inverseFFT
bool getCenterFreq(fftw_complex* centerFrequencies, fftw_complex* rbfCenterValues, double* rbfCenterPositions, int numCenters, int numFreq, int gridSize, double alpha, int interpFuncExtent, int rbfWidth, float min[3], float max[3], SmoothingFunction* smoothingFunction)
{
	int intAlpha=(int)alpha;
	// get the frequencies of the centers.
	{
		double* xk = new double[numCenters];
		double* yk = new double[numCenters];
		double* zk = new double[numCenters];
		double dmin[3] = {min[0], min[1], min[2]};
		double dmax[3] = {max[0], max[1], max[2]};
		if(!getNormalizedPositions(rbfCenterPositions, xk, yk, zk, rbfWidth, gridSize, numCenters, interpFuncExtent, numFreq, alpha, dmin, dmax))
		{
			return false;
		}
		if(!getCenterFrequencies(numCenters, alpha, interpFuncExtent, numFreq, gridSize, xk, yk, zk, rbfCenterValues, centerFrequencies, 0, smoothingFunction))
		{
			return false;
		}
		if(!fftshift3D(centerFrequencies, numFreq))
		{
			return false;
		}
		delete [] xk;
		delete [] yk;
		delete [] zk;
	}
	return true;
}

bool getRBFFreq(fftw_complex* rbf1D, int gridSize, fftw_complex* RBFFrequencies)
{
	if(!computeFFT1d(RBFFrequencies, rbf1D, gridSize, FFTW_FORWARD, FFTW_MEASURE))
	{
		cout<<"no good rbffreq out: ";
		for(int j=0; j<gridSize; j++)
		{
			cout<<RBFFrequencies[j][0]<<" ";
		}
		cout<<endl;
		return false;
	}
	return true;
}

bool multiplyWithTensor(fftw_complex* centerFrequencies, fftw_complex* RBFFrequencies, int gridSize, int numFreq)
{
	if(!multiplyWithLarger1DTensor(centerFrequencies, numFreq, RBFFrequencies, gridSize))
	{
		return false;
	}
	return true;
}

bool inverseFFT(fftw_complex* blurmap, fftw_complex* centerFrequencies, int numFreq, int gridSize)
{
	if(!zeroPad(blurmap, gridSize, centerFrequencies, numFreq))
	{
		return false;
	}
	if(!compute3DFFT(blurmap, blurmap, gridSize, gridSize, gridSize, FFTW_BACKWARD, FFTW_MEASURE))
	{
		return false;
	}
	return true;
}

bool getBlurMap(fftw_complex* blurmap, fftw_complex* rbfCenterValues, double* rbfCenterPositions, fftw_complex* rbf1D, int numCenters, int numFreq, int gridSize, int gridSizeBig, double alpha, int interpFuncExtent, int rbfWidth, float min[3], float max[3], SmoothingFunction* smoothingFunction)
{
	if(!blurmap || !rbfCenterValues || !rbfCenterPositions || !rbf1D || (numCenters<1) || (numFreq<1) || (gridSize<1) || (alpha<1) || (interpFuncExtent<1))
	{
		return false;
	}
	int intAlpha=(int)alpha;
	fftw_complex* centerFrequencies = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * numFreq*numFreq*numFreq);
	double forwardTime = 0;
	double ndftTime;
	// get the frequencies of the centers.
	{
		double* xk = new double[numCenters];
		double* yk = new double[numCenters];
		double* zk = new double[numCenters];
		double dmin[3] = {min[0], min[1], min[2]};
		double dmax[3] = {max[0], max[1], max[2]};
		if(!getNormalizedPositions(rbfCenterPositions, xk, yk, zk, rbfWidth, gridSize, numCenters, interpFuncExtent, numFreq, alpha, dmin, dmax))
		{
			return false;
		}
		double freqTime = getTime2();
		forwardTime = getTime2();
		ndftTime = getTime2();
		if(!getCenterFrequencies(numCenters, alpha, interpFuncExtent, numFreq, gridSize, xk, yk, zk, rbfCenterValues, centerFrequencies, 0, smoothingFunction))
		{
			return false;
		}
		if(!fftshift3D(centerFrequencies, numFreq))
		{
			return false;
		}
		delete [] xk;
		delete [] yk;
		delete [] zk;
	}
	fftw_complex* RBFFrequencies = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * gridSizeBig);
	if(!computeFFT1d(RBFFrequencies, rbf1D, gridSizeBig, FFTW_FORWARD, FFTW_MEASURE))
	{
		return false;
	}
	// perform multiplication with 1D tensor
	if(!multiplyWithLarger1DTensor(centerFrequencies, numFreq, RBFFrequencies, gridSizeBig))
	{
		return false;
	}
	fftw_free(RBFFrequencies);
	printf("Forward step took %lf \n", getTime2() - forwardTime);
	if(!zeroPad(blurmap, gridSizeBig, centerFrequencies, numFreq))
	{
		return false;
	}
	fftw_free(centerFrequencies);
	if(!compute3DFFT(blurmap, blurmap, gridSizeBig, gridSizeBig, gridSizeBig, FFTW_BACKWARD, FFTW_MEASURE))
	{
		return false;
	}
	printf("The NDFT took %lf seconds\n", getTime2() - ndftTime);
	return true;
}
