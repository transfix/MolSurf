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
#include <Blurmaps/CurvaturesGrid.h>

using PDBParser::Atom;
using PDBParser::GroupOfAtoms;


CurvaturesGrid::CurvaturesGrid()
{
	grid = 0;
}

CurvaturesGrid::~CurvaturesGrid()
{
	if(grid)
	{
		delete []grid;
		grid = 0;
	}
}

void CurvaturesGrid::getMinMax(GroupOfAtoms* molecule, int numberOfPoints, float* points)
{
	int i;
	for(i=0; i<numberOfPoints; i++)
	{
		int j;
		for(j=0; j<3; j++)
		{
			if(points[3*i+j] < m_Min[j])
			{
				m_Min[j] = points[3*i+j];
			}
			if(points[3*i+j] > m_Max[j])
			{
				m_Max[j] = points[3*i+j];
			}
		}
	}
}

void CurvaturesGrid::getIndices(double x, double y, double z, int* gridx, int* gridy, int* gridz)
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

int CurvaturesGrid::getIndex(double x, double y, double z)
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

void CurvaturesGrid::fillAtomsRecursively(GroupOfAtoms* molecule, double overlapExtent, int radiusType, unsigned int level)
{
	unsigned int i;
	for(i=0; i < molecule->m_SubGroups.size(); i++)
	{
		GroupOfAtoms* mol = molecule->m_SubGroups.at(i);
		fillAtomsRecursively(mol, overlapExtent, radiusType, level);
	}
	for(i=0; i < molecule->m_Atoms.size(); i++)
	{
		Atom* at = molecule->m_Atoms.at(i);
		double x = at->m_Position[0];
		double y = at->m_Position[1];
		double z = at->m_Position[2];
		int mingridx, mingridy, mingridz;
		int maxgridx, maxgridy, maxgridz;
		getIndices(x-overlapExtent, y-overlapExtent, z-overlapExtent, &mingridx, &mingridy, &mingridz);
		getIndices(x+overlapExtent, y+overlapExtent, z+overlapExtent, &maxgridx, &maxgridy, &maxgridz);
		int j, k, l;
		float rad = at->getRadius();
		for(j=mingridz; j<=maxgridz; j++)
		{
			for(k=mingridy; k<=maxgridy; k++)
			{
				for(l=mingridx; l<=maxgridx; l++)
				{
					int index = j*m_NumberOfGridDivisions*m_NumberOfGridDivisions + k*m_NumberOfGridDivisions + l;
					Tuple t(x, y, z, rad);
					grid[index].addAtom(t);
				}
			}
		}
	}
}

int getNearestInt(float x)
{
	int ret = (int)x;
	if((ret < x+0.001) && (ret > x-0.001))
	{
		return ret;
	}
	if((ret+1 < x+0.001) && (ret+1 > x-0.001))
	{
		return ret+1;
	}
	if((ret-1 < x+0.001) && (ret-1 > x-0.001))
	{
		return ret-1;
	}
	return ret; // not needed...
}

void CurvaturesGrid::createGrid(GroupOfAtoms* molecule, int numberOfPoints, float* points,
								double blobbyness, int numberOfGridDivisions, double maxFunctionError, int radiusType, unsigned int level)
{
	delete grid;
	grid = 0;
	m_NumberOfGridDivisions = numberOfGridDivisions;
	m_Blobbyness = blobbyness;
	m_NumberOfPoints = numberOfPoints;
	int i;
	for(i=0; i<3; i++)
	{
		m_Min[i] = 1000000000.0;
		m_Max[i] = -1000000000.0;
	}
	if(numberOfPoints < 1)
	{
		return;
	}
	fprintf(stderr, "Getting extents\n");
	getMinMax(molecule, numberOfPoints, points);
	grid = new CurvaturesGridVoxel[numberOfGridDivisions*numberOfGridDivisions*numberOfGridDivisions];
	fprintf(stderr, "Constructing grid\n");
	//fillGrid(molecule, numberOfPoints, points, maxFunctionError );
	double maxAtomRadius = 4.0;
	double overlapExtent = (log(maxFunctionError) + m_Blobbyness) * maxAtomRadius / m_Blobbyness;
	fillAtomsRecursively(molecule, overlapExtent, radiusType, level);
}

bool CurvaturesGrid::getCurvatures(double* HandK, double* normals, double* k1Vec, double* k2Vec, GroupOfAtoms* molecule, int numberOfPoints, float* points,
								   double blobbyness, int numberOfGridDivisions, double maxFunctionError,
								   int radiusType, unsigned int level)
{
	createGrid(molecule, numberOfPoints, points, blobbyness, numberOfGridDivisions, maxFunctionError, radiusType, level);
	int i;
	fprintf(stderr, "Getting curvatures\n");
	// for each grid g
	//		for each point p in g
	//			get contribution from each atom in g to p
	//			add contribution to HandK at correct index
	//		endfor
	// endfor
	for(i=0; i<numberOfPoints; i++)
	{
		HandK[i*2+0] = HandK[i*2+1] = 0;
		double x = points[3*i+0];
		double y = points[3*i+1];
		double z = points[3*i+2];
		int index = getIndex(x, y, z);
		//int numAtoms = grid[index].m_Atoms.count();
		int numAtoms = grid[index].getNumberOfAtoms();
		if(numAtoms < 1)
		{
			continue;
		}
		double phiX = 0;
		double phiY = 0;
		double phiZ = 0;
		double phiXX = 0;
		double phiYY = 0;
		double phiZZ = 0;
		double phiXY = 0;
		double phiXZ = 0;
		double phiYZ = 0;
		int k;
		for(k=0; k<numAtoms; k++)
		{
			//Tuple *a = grid[index].m_Atoms.at(k);
			Tuple a = grid[index].getAtom(k);
			double ax = a[0];
			double ay = a[1];
			double az = a[2];
			double rad = a[3];
			if(rad < 0.9)
			{
				//	printf("Rubbish\n");
			}
			evalCurvature(&phiX, &phiY, &phiZ, &phiXX, &phiYY, &phiZZ, &phiXY, &phiXZ, &phiYZ,
						  ax, ay, az, rad, x, y, z, i);
		}
		if((phiX*phiX + phiY*phiY + phiZ*phiZ) < 1e-30)
		{
			//	printf("This is too low, should not happen\n");
			continue;
		}
		// mean curvature
		HandK[i*2+0] = (phiX*phiX*(phiYY+phiZZ) +
						phiY*phiY*(phiXX+phiZZ) +
						phiZ*phiZ*(phiXX+phiYY) -
						2*phiX*phiY*phiXY -
						2*phiY*phiZ*phiYZ-
						2*phiX*phiZ*phiXZ) /
					   (2*pow(phiX*phiX + phiY*phiY + phiZ*phiZ, 1.5));
		// gaussian curvature
		HandK[i*2+1] = (phiZ*phiZ*(phiXX*phiYY-phiXY*phiXY) +
						phiY*phiY*(phiXX*phiZZ-phiXZ*phiXZ) +
						phiX*phiX*(phiYY*phiZZ-phiYZ*phiYZ) +
						2*phiX*phiY*(phiXZ*phiYZ - phiXY*phiZZ) +
						2*phiX*phiZ*(phiXY*phiYZ - phiXZ*phiYY) +
						2*phiY*phiZ*(phiXY*phiXZ - phiYZ*phiXX)) /
					   (pow(phiX*phiX + phiY*phiY + phiZ*phiZ, 2.0));
		double g = sqrt(phiX*phiX + phiY*phiY + phiZ*phiZ);
		// find k1, k2:
		double H = HandK[i*2+0], K = HandK[i*2+1];
		double H2minusK = H*H -K;
		if(H2minusK < 0)
		{
			H2minusK = 0;
		}
		double k1 = H - sqrt(H2minusK);
		double k2 = H + sqrt(H2minusK);
		// find normal:
		normals[i*3+0] = phiX / g;
		normals[i*3+1] = phiY / g;
		normals[i*3+2] = phiZ / g;
		double tk1Vec[3];
		double tk2Vec[3];
		tk1Vec[0] = tk1Vec[1] = tk1Vec[2] = tk2Vec[0] = tk2Vec[1] = tk2Vec[2] = 0;
		k1Vec[i*3+0] = k1Vec[i*3+1] = k1Vec[i*3+2] = 0;
		k2Vec[i*3+0] = k2Vec[i*3+1] = k2Vec[i*3+2] = 0;
		if(fabs(k1-k2) > 1e-10)
		{
			// construct curvature matrix. Find eigenvectors corresponding to k1 and k2.
			double c11, c12, c13;
			double c21, c22, c23;
			double c31, c32, c33;
			double g3 = g*g*g;
			c11 = (phiXX*g - phiX*(phiX*phiXX + phiY*phiXY + phiZ*phiXZ)) / g3;
			c12 = (phiXY*g - phiX*(phiX*phiXY + phiY*phiYY + phiZ*phiYZ)) / g3;
			c13 = (phiXZ*g - phiX*(phiX*phiXZ + phiY*phiYZ + phiZ*phiZZ)) / g3;
			c21 = (phiXY*g - phiY*(phiX*phiXX + phiY*phiXY + phiZ*phiXZ)) / g3;
			c22 = (phiYY*g - phiY*(phiX*phiXY + phiY*phiYY + phiZ*phiYZ)) / g3;
			c23 = (phiYZ*g - phiY*(phiX*phiXZ + phiY*phiYZ + phiZ*phiZZ)) / g3;
			c31 = (phiXZ*g - phiZ*(phiX*phiXX + phiY*phiXY + phiZ*phiXZ)) / g3;
			c32 = (phiYZ*g - phiZ*(phiX*phiXY + phiY*phiYY + phiZ*phiYZ)) / g3;
			c33 = (phiZZ*g - phiZ*(phiX*phiXZ + phiY*phiYZ + phiZ*phiZZ)) / g3;
			if(CCVOpenGLMath::LinearAlgebra::solveEigenSystem(c11, c12, c13, c21, c22, c23, c31, c32, c33, tk1Vec, tk2Vec, k1, k2))
			{
				double d;
				d = sqrt(tk1Vec[0]*tk1Vec[0] + tk1Vec[1]*tk1Vec[1] + tk1Vec[2]*tk1Vec[2]);
				tk1Vec[0] /= d;
				tk1Vec[1] /= d;
				tk1Vec[2] /= d;
				d = sqrt(tk2Vec[0]*tk2Vec[0] + tk2Vec[1]*tk2Vec[1] + tk2Vec[2]*tk2Vec[2]);
				tk2Vec[0] /= d;
				tk2Vec[1] /= d;
				tk2Vec[2] /= d;
				k1Vec[i*3+0] = tk1Vec[0];
				k1Vec[i*3+1] = tk1Vec[1];
				k1Vec[i*3+2] = tk1Vec[2];
				k2Vec[i*3+0] = tk2Vec[0];
				k2Vec[i*3+1] = tk2Vec[1];
				k2Vec[i*3+2] = tk2Vec[2];
			}
		}
		if(((i+1) % 20) == 0 || (i+1) == numberOfPoints)
		{
			printf("Evaluating functions %2.2f%% done (%012d)\r", 100.0*(i+1)/(float)(numberOfPoints), i+1);
		}
	}
	printf("\n");
	return true;
}

inline void CurvaturesGrid::evalCurvature(double* phiX, double* phiY, double* phiZ, double* phiXX,
		double* phiYY, double* phiZZ,  double* phiXY, double* phiXZ,
		double* phiYZ,  double xc, double yc, double zc, double rad,
		double x, double y, double z, int vertexNumber)
{
	double r =	(x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc);
	double expval = m_Blobbyness*r/rad - m_Blobbyness;
	double phi = exp(expval);
	double f = 2*m_Blobbyness/rad;
	if(phi < 1e-30)
	{
		phi = 0;
	}
	if((fabs(f)>.00001) && phi > 1e-20)
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
	/*	if( vertexNumber == 329 )
		{
			FILE* fp = fopen("whatever.txt","a");
			fprintf( fp, "phi = %lf, *phiX = %lf, point=[%lf %lf %lf], center=[%lf %lf %lf]\n", phi, (*phiX), x, y, z, xc, yc, zc );
			fclose( fp);
		}*/
}

bool CurvaturesGrid::write(int numberOfPoints, double* HandK, double* normals, double* k1Vec, double* k2Vec, const char* filename)
{
	if(!HandK || !normals || !k1Vec || !k2Vec || (numberOfPoints < 1))
	{
		return false;
	}
	FILE* fp = fopen(filename, "w");
	if(!fp)
	{
		return false;
	}
	fprintf(fp, "#\n");
	fprintf(fp, "# Curvature file exported from TexMol.\n");
	fprintf(fp, "# CVC, The University of Texas at Austin.\n");
	fprintf(fp, "#\n");
	fprintf(fp, "# Format:\n");
	fprintf(fp, "# First see the statement: Number: number of points\n");
	fprintf(fp, "# Then, search for one of the following: Mean, Gaussian, Normal, K1Vec, K2Vec\n");
	fprintf(fp, "# Any of the above should be followed by the corresponding values with 1, 1, 3, 3, 3 in each line\n");
	fprintf(fp, "#\n");
	fprintf(fp, "#(All comments start with '#')\n");
	fprintf(fp, "#\n");
	fprintf(fp, "Number: %d\n", numberOfPoints);
	int i;
	fprintf(fp, "Mean\n");
	for(i=0; i<numberOfPoints; i++)
	{
		fprintf(fp, "%lf\n", HandK[i*2+0]);
	}
	fprintf(fp, "Gaussian\n");
	for(i=0; i<numberOfPoints; i++)
	{
		fprintf(fp, "%lf\n", HandK[i*2+1]);
	}
	fprintf(fp, "Normal\n");
	for(i=0; i<numberOfPoints; i++)
	{
		fprintf(fp, "%lf %lf %lf\n", normals[i*3+0], normals[i*3+1], normals[i*3+2]);
	}
	fprintf(fp, "K1Vec\n");
	for(i=0; i<numberOfPoints; i++)
	{
		fprintf(fp, "%lf %lf %lf\n", k1Vec[i*3+0], k1Vec[i*3+1], k1Vec[i*3+2]);
	}
	fprintf(fp, "K2Vec\n");
	for(i=0; i<numberOfPoints; i++)
	{
		fprintf(fp, "%lf %lf %lf\n", k2Vec[i*3+0], k2Vec[i*3+1], k2Vec[i*3+2]);
	}
	fclose(fp);
	return true;
}
