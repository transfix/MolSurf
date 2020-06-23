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
#include<Blurmaps/InterfaceSurface.h>
#include <DPG/PG.h>

InterfaceSurface::InterfaceSurface(string sf1, string sf2, double iw, bool ci)
{
	surfaceFileName1 = sf1;
	surfaceFileName2 = sf2;
	interfaceWidth = iw;
	coloredInput = ci;

	size_t found = sf1.find('.');
  	
	interfaceFileName.assign(sf1, 0, (int)found);
	interfaceFileName += "_merged_interface.rawnc";

	highlightedSurfaceFileName1.assign(sf1, 0, (int)found);
	highlightedSurfaceFileName1 += "_interface.rawnc";

	interIndexFileName1.assign(sf1,0, (int)found);
	interIndexFileName1 += "_interface_points.txt";

	found = sf2.find('.');
	highlightedSurfaceFileName2.assign(sf2, 0, (int)found);
	highlightedSurfaceFileName2 += "_interface.rawnc";

	interIndexFileName2.assign(sf2,0, (int)found);
	interIndexFileName2 += "_interface_points.txt";
}

InterfaceSurface::~InterfaceSurface()
{

}

bool InterfaceSurface::detect(double & area1, double & area2)
{
/* Opening input files */

	FILE *surfaceFile1, *surfaceFile2;

	surfaceFile1 = fopen(surfaceFileName1.c_str(), "rt");
	if(!surfaceFile1)
	{
		cout<<"could not open surface file "<< surfaceFileName1.c_str() << endl;
		return false;
	}

	surfaceFile2 = fopen(surfaceFileName2.c_str(), "rt");
	if(!surfaceFile2)
	{
		cout<<"could not open surface file "<< surfaceFileName2.c_str() << endl;
		fclose(surfaceFile1);
		return false;
	}

	
	FILE *interfaceFile1, *interfaceFile2;

	interfaceFile1 = fopen(interIndexFileName1.c_str(), "wt");
	if(!interfaceFile1)
	{
		cout<<"could not open surface file "<< interIndexFileName1.c_str() << endl;
		return false;
	}

	interfaceFile2 = fopen(interIndexFileName2.c_str(), "wt");
	if(!interfaceFile2)
	{
		cout<<"could not open surface file "<< interIndexFileName2.c_str() << endl;
		fclose(interfaceFile1);
		return false;
	}

/* Reading points from the first surface file and inserting into DPG */

	int retval = 0;
	int numVert1, numTriang1;
	retval = fscanf(surfaceFile1, "%d", &numVert1);
	retval = fscanf(surfaceFile1, "%d", &numTriang1);

	DPG::PG *pg1 = new DPG::PG(10.0, 1000.0, 3.0);

	double *xArray1 = new double[numVert1];
	double *yArray1 = new double[numVert1];
	double *zArray1 = new double[numVert1];
	double *nxArray1 = new double[numVert1];
	double *nyArray1 = new double[numVert1];
	double *nzArray1 = new double[numVert1];
	double *rArray1 = new double[numVert1];
	double *gArray1 = new double[numVert1];
	double *bArray1 = new double[numVert1];

	bool *interface1 = new bool[numVert1];

	for(int i=0; i<numVert1; i++)
	{
		retval = fscanf(surfaceFile1, "%lf", &xArray1[i]);
		retval = fscanf(surfaceFile1, "%lf", &yArray1[i]);
		retval = fscanf(surfaceFile1, "%lf", &zArray1[i]);
		retval = fscanf(surfaceFile1, "%lf", &nxArray1[i]);
		retval = fscanf(surfaceFile1, "%lf", &nyArray1[i]);
		retval = fscanf(surfaceFile1, "%lf", &nzArray1[i]);

		if(coloredInput)
		{
			retval = fscanf(surfaceFile1, "%lf", &rArray1[i]);
			retval = fscanf(surfaceFile1, "%lf", &gArray1[i]);
			retval = fscanf(surfaceFile1, "%lf", &bArray1[i]);
		}

		interface1[i] = false;

		DPG::Point *p = new DPG::Point(xArray1[i], yArray1[i], zArray1[i]);
		pg1->addPoint(p);
	}


/* Reading points from the second surface file and inserting into DPG */

	int numVert2, numTriang2;
	retval = fscanf(surfaceFile2, "%d", &numVert2);
	retval = fscanf(surfaceFile2, "%d", &numTriang2);

	DPG::PG *pg2 = new DPG::PG(10.0, 1000.0, 3.0);

	double *xArray2 = new double[numVert2];
	double *yArray2 = new double[numVert2];
	double *zArray2 = new double[numVert2];
	double *nxArray2 = new double[numVert2];
	double *nyArray2 = new double[numVert2];
	double *nzArray2 = new double[numVert2];
	double *rArray2 = new double[numVert2];
	double *gArray2 = new double[numVert2];
	double *bArray2 = new double[numVert2];

	bool *interface2 = new bool[numVert2];

	for(int i=0; i<numVert2; i++)
	{
		retval = fscanf(surfaceFile2, "%lf", &xArray2[i]);
		retval = fscanf(surfaceFile2, "%lf", &yArray2[i]);
		retval = fscanf(surfaceFile2, "%lf", &zArray2[i]);
		retval = fscanf(surfaceFile2, "%lf", &nxArray2[i]);
		retval = fscanf(surfaceFile2, "%lf", &nyArray2[i]);
		retval = fscanf(surfaceFile2, "%lf", &nzArray2[i]);

		if(coloredInput)
		{
			retval = fscanf(surfaceFile2, "%lf", &rArray2[i]);
			retval = fscanf(surfaceFile2, "%lf", &gArray2[i]);
			retval = fscanf(surfaceFile2, "%lf", &bArray2[i]);
		}

		interface2[i] = false;

		DPG::Point *p = new DPG::Point(xArray2[i], yArray2[i], zArray2[i]);
		pg2->addPoint(p);
	}


	cout<<"Read input surfaces"<<endl;


/* Opening output files */

	FILE *highlightedSurfaceFile1, *highlightedSurfaceFile2, *interfaceFile;

	highlightedSurfaceFile1 = fopen(highlightedSurfaceFileName1.c_str(), "wt");
	if(!highlightedSurfaceFile1)
	{
		cout<<"could not open surface file "<< highlightedSurfaceFileName1.c_str() << endl;
		return false;
	}

	highlightedSurfaceFile2 = fopen(highlightedSurfaceFileName2.c_str(), "wt");
	if(!highlightedSurfaceFile2)
	{
		cout<<"could not open surface file "<< highlightedSurfaceFileName2.c_str() << endl;
		fclose(highlightedSurfaceFile1);
		return false;
	}

/*	interfaceFile = fopen(interfaceFileName.c_str(), "rt");
	if(!interfaceFile)
	{
		cout<<"could not open surface file "<< interfaceFileName.c_str() << endl;
		fclose(highlightedSurfaceFile1);
		fclose(highlightedSurfaceFile2);
		return false;
	}*/


/* Identifying interface points on the first surface and writing to highlightedSurfaceFile1 */

	fprintf(highlightedSurfaceFile1, "%d\n", numVert1);
	fprintf(highlightedSurfaceFile1, "%d\n", numTriang1);

	for(int i=0; i<numVert1; i++)
	{
		fprintf(highlightedSurfaceFile1, "%lf %lf %lf %lf %lf %lf", xArray1[i], yArray1[i], zArray1[i], nxArray1[i], nyArray1[i], nzArray1[i]);

		DPG::Point *p = new DPG::Point(xArray1[i], yArray1[i], zArray1[i]);

		if(pg2->pointsWithinRange(p, interfaceWidth))	// point is on the interface. (re)color it
		{
			interface1[i] = true;

			fprintf(interfaceFile1, "%d\n", i);

			if(coloredInput)
			{
				double r = rArray1[i]+0.5;
				if(r>1.0) r = 1.0;
				double g = gArray1[i]+0.5;
				if(g>1.0) g = 1.0;
				double b = bArray1[i];
				if(b>1.0) b = 1.0;	

//				r = 1.0; b = 0.5; g = 0.5;

				fprintf(highlightedSurfaceFile1, " %lf %lf %lf\n", r, g, b);
			}
			else
				fprintf(highlightedSurfaceFile1, " %lf %lf %lf\n", 0.5, 0.0, 0.0);
		}
		
		else
		{
			if(coloredInput)
			{
				double r = rArray1[i];
				double g = gArray1[i];
				double b = bArray1[i];

				fprintf(highlightedSurfaceFile1, " %lf %lf %lf\n", r, g, b);
			}
			else
				fprintf(highlightedSurfaceFile1, " %lf %lf %lf\n", 0.0, 0.5, 0.5);
		}
	}

	area1 = 0;

	cout<<"Computed interface for the first molecule"<<endl;	

	for(int i=0; i<numTriang1; i++)	//copy over the triangles and compute the area
	{
		int index1, index2, index3;
		retval = fscanf(surfaceFile1, "%d", &index1);
		retval = fscanf(surfaceFile1, "%d", &index2);
		retval = fscanf(surfaceFile1, "%d", &index3);

		fprintf(highlightedSurfaceFile1, "%d %d %d\n", index1, index2, index3);

		if(interface1[index1] && interface1[index2] && interface1[index3])
		{
			double dx, dy, dz;
			double a, b, c;
			double s;
			double area;

			dx = xArray1[index2] - xArray1[index1];
			dy = yArray1[index2] - yArray1[index1];
			dz = zArray1[index2] - zArray1[index1];
			a = sqrt(dx*dx + dy*dy + dz*dz);

			dx = xArray1[index3] - xArray1[index2];
			dy = yArray1[index3] - yArray1[index2];
			dz = zArray1[index3] - zArray1[index2];
			b = sqrt(dx*dx + dy*dy + dz*dz);

			dx = xArray1[index1] - xArray1[index3];
			dy = yArray1[index1] - yArray1[index3];
			dz = zArray1[index1] - zArray1[index3];
			c = sqrt(dx*dx + dy*dy + dz*dz);

			s = (a+b+c)/2.0;

			area = sqrt(s*(s-a)*(s-b)*(s-c));

			area1 += area;
		}
	}

	cout<<"Computed interface area for the first molecule"<<endl;	


/* Identifying interface points on the second surface and writing to highlightedSurfaceFile2 */

	fprintf(highlightedSurfaceFile2, "%d\n", numVert2);
	fprintf(highlightedSurfaceFile2, "%d\n", numTriang2);

	for(int i=0; i<numVert2; i++)
	{
		fprintf(highlightedSurfaceFile2, "%lf %lf %lf %lf %lf %lf", xArray2[i], yArray2[i], zArray2[i], nxArray2[i], nyArray2[i], nzArray2[i]);

		DPG::Point *p = new DPG::Point(xArray2[i], yArray2[i], zArray2[i]);

		if(pg1->pointsWithinRange(p, interfaceWidth))	// point is on the interface. (re)color it
		{
			interface2[i] = true;

			fprintf(interfaceFile2, "%d\n", i);

			if(coloredInput)
			{
				double r = rArray2[i];
				if(r>1.0) r = 1.0;
				double g = gArray2[i]+0.5;
				if(g>1.0) g = 1.0;
				double b = bArray2[i]+0.5;
				if(b>1.0) b = 1.0;
				
//				r = 0.5; b = 0.5; g = 1.0;

				fprintf(highlightedSurfaceFile2, " %lf %lf %lf\n", r, g, b);
			}
			else
				fprintf(highlightedSurfaceFile2, " %lf %lf %lf\n", 0.0, 0.0, 0.5);
		}
		
		else
		{
			if(coloredInput)
			{
				double r = rArray2[i];
				double g = gArray2[i];
				double b = bArray2[i];

				fprintf(highlightedSurfaceFile2, " %lf %lf %lf\n", r, g, b);
			}
			else
				fprintf(highlightedSurfaceFile2, " %lf %lf %lf\n", 0.5, 0.5, 0.0);
		}
	}

	cout<<"Computed interface for the second molecule"<<endl;	

	area2 = 0;

	for(int i=0; i<numTriang2; i++)	//copy over the triangles and compute the area
	{
		int index1, index2, index3;
		retval = fscanf(surfaceFile2, "%d", &index1);
		retval = fscanf(surfaceFile2, "%d", &index2);
		retval = fscanf(surfaceFile2, "%d", &index3);

		fprintf(highlightedSurfaceFile2, "%d %d %d\n", index1, index2, index3);

		if(interface2[index1] && interface2[index2] && interface2[index3])
		{
			double dx, dy, dz;
			double a, b, c;
			double s;
			double area;

			dx = xArray2[index2] - xArray2[index1];
			dy = yArray2[index2] - yArray2[index1];
			dz = zArray2[index2] - zArray2[index1];
			a = sqrt(dx*dx + dy*dy + dz*dz);

			dx = xArray2[index3] - xArray2[index2];
			dy = yArray2[index3] - yArray2[index2];
			dz = zArray2[index3] - zArray2[index2];
			b = sqrt(dx*dx + dy*dy + dz*dz);

			dx = xArray2[index1] - xArray2[index3];
			dy = yArray2[index1] - yArray2[index3];
			dz = zArray2[index1] - zArray2[index3];
			c = sqrt(dx*dx + dy*dy + dz*dz);

			s = (a+b+c)/2.0;

			area = sqrt(s*(s-a)*(s-b)*(s-c));

			area2 += area;
		}
	}

	cout<<"Computed interface area for the second molecule"<<endl;	

/* Cleaning up */

	if(surfaceFile1) fclose(surfaceFile1);
	if(surfaceFile2) fclose(surfaceFile2);
	if(highlightedSurfaceFile1) fclose(highlightedSurfaceFile1);
	if(highlightedSurfaceFile2) fclose(highlightedSurfaceFile2);
	if(interfaceFile1) fclose(interfaceFile1);
	if(interfaceFile2) fclose(interfaceFile2);

	delete [] xArray1;
	delete [] yArray1;
	delete [] zArray1;
	delete [] nxArray1;
	delete [] nyArray1;
	delete [] nzArray1;
	delete [] rArray1;
	delete [] gArray1;
	delete [] bArray1;
	delete [] interface1;

	delete [] xArray2;
	delete [] yArray2;
	delete [] zArray2;
	delete [] nxArray2;
	delete [] nyArray2;
	delete [] nzArray2;
	delete [] rArray2;
	delete [] gArray2;
	delete [] bArray2;
	delete [] interface2;

	return true;
}


