/*
  Copyright 2011 The University of Texas at Austin

	Authors: Muhibur Rasheed <muhibur@ices.utexas.edu>
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
#include <InterfaceStats/ComputeLeastSquarePlane.h>
#include <cstdlib>
#include <cmath>
#include <iostream>

extern "C" {void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda,
                double* w, double* work, int* lwork, int* info );
	   }
#define N 3
#define LDA N

ComputeLeastSquarePlane::ComputeLeastSquarePlane()
{

}

ComputeLeastSquarePlane::~ComputeLeastSquarePlane()
{

}

bool ComputeLeastSquarePlane::compute(vector < DPG::Point* > & pts, InterfaceStats* is)
{
	intStat = is;

	if(computePlane(pts))
	{
		computeStats(pts);
		return true;
	}
	
	return false;
}

bool ComputeLeastSquarePlane::computePlane(vector < DPG::Point* > & pts)
{
	int numPoints = pts.size();
	double centroid[3];

	// Compute the centroid
	
	centroid[0] = 0.0, centroid[1] = 0.0, centroid[2] = 0.0;

	for(int i=0; i<numPoints; i++)
	{
		centroid[0] += pts[i]->getX();
		centroid[1] += pts[i]->getY();
		centroid[2] += pts[i]->getZ();
	}

	centroid[0] /= (double)numPoints;
	centroid[1] /= (double)numPoints;
	centroid[2] /= (double)numPoints;

	// compute the symmetric matrix A

	double a00 = 0.0, a01 = 0.0, a02 = 0.0, a10 = 0.0, a11 = 0.0, a12 = 0.0, a20 = 0.0, a21 = 0.0, a22 = 0.0;

	for(int i=0; i<numPoints; i++)
	{
		double xc, yc, zc;
		xc = pts[i]->getX() - centroid[0];
		yc = pts[i]->getY() - centroid[1];
		zc = pts[i]->getZ() - centroid[2];

		a00 += xc*xc;
		a01 += xc*yc;
		a02 += xc*zc;
		a11 += yc*yc;		
		a12 += yc*zc;
		a22 += zc*zc;
	}

	a10 = a01;
	a20 = a02;
	a21 = a12;

	// set up and execute LAPACK's eigenvalue solver

	int n = N, lda = LDA, info, lwork;
	double wkopt;
	double* work;
	double w[N];

	double a[LDA*N] = {	a00, 	0.00, 	0.00,
			   	a10, 	a11, 	0.0,
				a20, 	a21, 	a22 };

	lwork = -1;
	dsyev_("Vectors", "Upper", &n, a, &lda, w, &wkopt, &lwork, &info);
	lwork = (int)wkopt;
	work = (double*)malloc(lwork*sizeof(double));
	dsyev_("Vectors", "Upper", &n, a, &lda, w, work, &lwork, &info);

	
	if(info > 0)
	{
		cout << "The algorithm failed to compute eigenvalues.\n" <<endl;
		return false;
	}
	else 
	{   
		// w[0] is the minimal eignevalue, a[0-2] is corresponding unit eigen vector. This is the least squares solution

		eigenValues[0] = w[0];
		eigenValues[1] = w[1];
		eigenValues[2] = w[2];

		plane[0] = a[0];
		plane[1] = a[1];
		plane[2] = a[2];
		plane[3] = -(a[0]*centroid[0] + a[1]*centroid[1] + a[2]*centroid[2]);

		return true;
	}
}

void ComputeLeastSquarePlane::computeStats(vector < DPG::Point* > & pts)
{
	double dist;		// distance from a point to the plane
	double rmsd = 0;

	double x,y,z;
	double normalizingDenumerator = sqrt(plane[0]*plane[0] + plane[1]*plane[1] + plane[2]*plane[2]);

	int numPoints = pts.size();

	for(int i=0; i<numPoints; i++)
	{
		DPG::Point* temp = pts[i];
		x = temp->getX();
		y = temp->getY();
		z = temp->getZ();

		dist = fabs((plane[0]*x + plane[1]*y + plane[2]*z + plane[3]) / normalizingDenumerator);

		rmsd += dist*dist;
	}
	
	rmsd = sqrt(rmsd/(double)numPoints);

	intStat->setPlanarityIndex(rmsd);
	intStat->setCircularityIndex(eigenValues[1]/eigenValues[2]);	// ratio of the minor/major components on the plane
}
