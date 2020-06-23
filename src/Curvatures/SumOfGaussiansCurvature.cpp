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
#include <Curvatures/SumOfGaussiansCurvature.h>

SumOfGaussiansCurvature::SumOfGaussiansCurvature(int numberOfGaussians, double* gaussianCenters, int numberOfGridDivisions, double maxFunctionError, double blobbiness, int numberOfPoints, float* points, double* HandK, double* normals, double* k1Vector, double* k2Vector)
	: Curvature(numberOfPoints, points, HandK, normals, k1Vector, k2Vector)
{
	m_NumberOfGaussians = numberOfGaussians;
	m_GaussianCenters = gaussianCenters;
	m_NumberOfGridDivisions = numberOfGridDivisions;
	m_MaxFunctionError = maxFunctionError;
	m_Blobbiness = blobbiness;
	m_CurvaturesGridVoxels = 0;
	m_Min[0] = 0;
	m_Min[1] = 0;
	m_Min[2] = 0;
	m_Max[0] = 1;
	m_Max[1] = 1;
	m_Max[2] = 1;
}

SumOfGaussiansCurvature::~SumOfGaussiansCurvature()
{
	if(m_CurvaturesGridVoxels)
	{
		delete []m_CurvaturesGridVoxels;
		m_CurvaturesGridVoxels = 0;
	}
}

void SumOfGaussiansCurvature::getMinMax()
{
	int i;
	for(i=0; i<m_NumberOfPoints; i++)
	{
		int j;
		for(j=0; j<3; j++)
		{
			if(m_Points[3*i+j] < m_Min[j])
			{
				m_Min[j] = m_Points[3*i+j];
			}
			if(m_Points[3*i+j] > m_Max[j])
			{
				m_Max[j] = m_Points[3*i+j];
			}
		}
	}
}

void SumOfGaussiansCurvature::getIndices(double x, double y, double z, int* gridx, int* gridy, int* gridz)
{
	*gridx = (x-m_Min[0]) / (m_Max[0] - m_Min[0]) * m_NumberOfGridDivisions;
	if(*gridx < 0)
	{
		*gridx = 0;
	}
	if(*gridx > m_NumberOfGridDivisions-1)
	{
		*gridx = m_NumberOfGridDivisions-1;
	}
	*gridy = (y-m_Min[1]) / (m_Max[1] - m_Min[1]) * m_NumberOfGridDivisions;
	if(*gridy < 0)
	{
		*gridy = 0;
	}
	if(*gridy > m_NumberOfGridDivisions-1)
	{
		*gridy = m_NumberOfGridDivisions-1;
	}
	*gridz = (z-m_Min[2]) / (m_Max[2] - m_Min[2]) * m_NumberOfGridDivisions;
	if(*gridz < 0)
	{
		*gridz = 0;
	}
	if(*gridz > m_NumberOfGridDivisions-1)
	{
		*gridz = m_NumberOfGridDivisions-1;
	}
}

void SumOfGaussiansCurvature::populateGrid(double overlapExtent)
{
	unsigned int i;
	for(i=0; i < m_NumberOfGaussians; i++)
	{
		double x = m_GaussianCenters[i*4+0];
		double y = m_GaussianCenters[i*4+1];
		double z = m_GaussianCenters[i*4+2];
		double r = m_GaussianCenters[i*4+3];
		int mingridx, mingridy, mingridz;
		int maxgridx, maxgridy, maxgridz;
		getIndices(x-overlapExtent, y-overlapExtent, z-overlapExtent, &mingridx, &mingridy, &mingridz);
		getIndices(x+overlapExtent, y+overlapExtent, z+overlapExtent, &maxgridx, &maxgridy, &maxgridz);
		int j, k, l;
		for(j=mingridz; j<=maxgridz; j++)
		{
			for(k=mingridy; k<=maxgridy; k++)
			{
				for(l=mingridx; l<=maxgridx; l++)
				{
					int index = j*m_NumberOfGridDivisions*m_NumberOfGridDivisions + k*m_NumberOfGridDivisions + l;
					Tuple t(x, y, z, r);
					m_CurvaturesGridVoxels[index].addKernel(t);
				}
			}
		}
	}
}

void SumOfGaussiansCurvature::createGrid()
{
	if(m_CurvaturesGridVoxels)
	{
		delete []m_CurvaturesGridVoxels;
		m_CurvaturesGridVoxels = 0;
	}
	int i;
	for(i=0; i<3; i++)
	{
		m_Min[i] = 1000000000.0;
		m_Max[i] = -1000000000.0;
	}
	getMinMax();
	m_CurvaturesGridVoxels = new CurvaturesGridVoxel[m_NumberOfGridDivisions*m_NumberOfGridDivisions*m_NumberOfGridDivisions];
	double maxAtomRadius = 4.0;
	double overlapExtent = (log(m_MaxFunctionError) + m_Blobbiness) * maxAtomRadius / m_Blobbiness;
	populateGrid(overlapExtent);
}

bool SumOfGaussiansCurvature::initialize()
{
	createGrid();
	m_Initialized = true;
	return true;
}

inline void SumOfGaussiansCurvature::evalCurvatureDueToOneAtom(double* phiX, double* phiY, double* phiZ, double* phiXX,
		double* phiYY, double* phiZZ,  double* phiXY, double* phiXZ,
		double* phiYZ,  double xc, double yc, double zc, double rad,
		double x, double y, double z)
{
	double r =	(x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc);
	double expval = m_Blobbiness*r/(rad*rad) - m_Blobbiness;
	double phi = exp(expval);
	double f = 2*m_Blobbiness/(rad*rad);
	if(phi > 1e-20)
	{
		*phiX += f*(x-xc)*phi;
		*phiY += f*(y-yc)*phi;
		*phiZ += f*(z-zc)*phi;
		*phiXY += f*f*(x-xc)*(y-yc)*phi;
		*phiXZ += f*f*(x-xc)*(z-zc)*phi;
		*phiYZ += f*f*(y-yc)*(z-zc)*phi;
		*phiXX += f*f*(x-xc)*(x-xc)*phi + f*phi;
		*phiYY += f*f*(y-yc)*(y-yc)*phi + f*phi;
		*phiZZ += f*f*(z-zc)*(z-zc)*phi + f*phi;
	}
}

inline int SumOfGaussiansCurvature::getIndex(double x, double y, double z)
{
	int gridx, gridy, gridz;
	gridx = (x-m_Min[0]) / (m_Max[0] - m_Min[0]) * m_NumberOfGridDivisions;
	if(gridx < 0)
	{
		gridx = 0;
	}
	if(gridx > m_NumberOfGridDivisions-1)
	{
		gridx = m_NumberOfGridDivisions-1;
	}
	gridy = (y-m_Min[1]) / (m_Max[1] - m_Min[1]) * m_NumberOfGridDivisions;
	if(gridy < 0)
	{
		gridy = 0;
	}
	if(gridy > m_NumberOfGridDivisions-1)
	{
		gridy = m_NumberOfGridDivisions-1;
	}
	gridz = (z-m_Min[2]) / (m_Max[2] - m_Min[2]) * m_NumberOfGridDivisions;
	if(gridz < 0)
	{
		gridz = 0;
	}
	if(gridz > m_NumberOfGridDivisions-1)
	{
		gridz = m_NumberOfGridDivisions-1;
	}
	return gridz*m_NumberOfGridDivisions*m_NumberOfGridDivisions + gridy*m_NumberOfGridDivisions + gridx;
}

inline void SumOfGaussiansCurvature::evalCurvature(double* phiX,  double* phiY,  double* phiZ,
		double* phiXX, double* phiYY, double* phiZZ,
		double* phiXY, double* phiXZ, double* phiYZ,
		double x,      double y,      double z)
{
	int index = getIndex(x, y, z);
	int numAtoms = m_CurvaturesGridVoxels[index].getNumberOfKernels();
	if(numAtoms < 1)
	{
		return;
	}
	int k;
	for(k=0; k<numAtoms; k++)
	{
		Tuple a = m_CurvaturesGridVoxels[index].getKernel(k);
		double ax = a[0];
		double ay = a[1];
		double az = a[2];
		double rad = a[3];
		evalCurvatureDueToOneAtom(phiX, phiY, phiZ, phiXX, phiYY, phiZZ, phiXY, phiXZ, phiYZ,
								  ax, ay, az, rad, x, y, z);
	}
}
