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

void getGaussian(double b, int grid_size, double error, double width, double r, fftw_complex* rbf, int& rbfWidth)
{
	double ratio = (double) grid_size / (double) width;
	double r0 = r * r;
	int c=0;
	while(1)
	{
		double dist = c * c / (ratio * ratio);
		double expVal = -b * dist /r/r + b;
		if(exp(expVal) < error)
		{
			c++;
			break;
		}
		c++;
	}
	/* some extra stuff; can be found at blurring.cpp */
	for(int g=0; g<grid_size; g++)
	{
		rbf[g][0] = rbf[g][1] = 0;
	}
	int num = c + 1;
	if(num <= grid_size / 2)
	{
		rbfWidth = num;
	}
	else
	{
		rbfWidth = grid_size /2;
	}
	for(c=0; c<rbfWidth; c++)
	{
		double dist = c * c / (ratio * ratio);
		double expVal = -b * dist / r /r + b;
		rbf[c][0] = exp(expVal) / exp(b);
		rbf[c][1] = 0;
	}
	int i;
	for(i=grid_size - 1; i>= grid_size - 1 - rbfWidth +1; i--)
	{
		rbf[i][0] = rbf[grid_size -1 - i + 1][0];
		rbf[i][1] = rbf[grid_size -1 - i + 1][1];
	}
}


void getGaussian2(int grid_size, double error, double width, double r, double b, fftw_complex* rbf, int& rbfWidth)
{
	double ratio = (double) width/ (double) grid_size;
	double r0 = r * r;
	int c=0;
	while(1)
	{
		double dist = c * c * ratio * ratio;
		double expVal = -b * dist /r/r + b;
		if(exp(expVal) < error)
		{
			c++;
			break;
		}
		c++;
	}
	/* some extra stuff; can be found at blurring.cpp */
	for(c=0; c<grid_size; c++)
	{
		rbf[c][0] = rbf[c][1] = 0;
	}
	int num = c;
	for(c=0; c<num && c<grid_size /2; c++)
	{
		double dist = c * c * (ratio * ratio);
		double expVal = -b * dist / r /r + b;
		rbf[c][0] = exp(expVal) / exp(b);
		rbf[c][1] = 0;
	}
	if(num +1 <= grid_size / 2)
	{
		rbfWidth = num +1;
	}
	else
	{
		rbfWidth = grid_size /2;
	}
	int i;
	for(i=grid_size - 1; i>= grid_size - 1 - rbfWidth +1; i--)
	{
		rbf[i][0] = rbf[grid_size -1 - i +1][0];
		rbf[i][1] = rbf[grid_size -1 - i +1][1];
	}
}
