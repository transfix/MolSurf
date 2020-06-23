/*
  Copyright 2011 The University of Texas at Austin

        Author: Muhibur Rasheed <muhib@ices.utexas.edu>
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
#include <Hausdorff/Hausdorff.h>

#include <cstdio>

bool HausdorffDistance::readFile(string fileName1, string fileName2, bool isRAWN)
{
	FILE* fpX = fopen(fileName1.c_str(), "rt");
	FILE* fpY = fopen(fileName2.c_str(), "rt");

	if(fpX!=NULL && fpY!=NULL)
	{
		int temp;
		double px, py, pz;
		double tempnx, tempny, tempnz;

		double minx = 0, miny = 0, minz = 0;

//		cout<<"Reading files "<<fileName1<<", "<<fileName2<<endl;

		if( fscanf(fpX, "%d %d", &numPointsX, &temp) != 2)
		{
			cout<<"Error reading number of vertex for file 1"<<endl;
			return false;
		}		
		if( fscanf(fpY, "%d %d", &numPointsY, &temp) != 2)
		{
			cout<<"Error reading number of vertex for file 2"<<endl;
			return false;
		}		

//		cout<<"mesh 1 has "<<numPointsX<<" points"<<endl;
//		cout<<"mesh 2 has "<<numPointsY<<" points"<<endl;

		for(int i=0; i<numPointsX; i++)
		{
			if(isRAWN)
			{ 
				if(fscanf(fpX, "%lf %lf %lf %lf %lf %lf", &px, &py, &pz, &tempnx, &tempny, &tempnz) != 6)
				{
					cout<<"Error reading vertex "<<i<<endl;
					return false;
				}
			}
			else
			{
				if(fscanf(fpX, "%lf %lf %lf", &px, &py, &pz) != 3)
				{
					cout<<"Error reading vertex "<<i<<endl;
					return false;
				}
			}

			DPG::Point* p = new DPG::Point(px, py, pz);
			pointSetX.push_back(p);

			if(minx>px) minx = px;
			if(miny>py) miny = py;
			if(minz>pz) minz = pz;
		}

//		cout<<"Read all points of mesh 1"<<endl;

		double min = minx;
		if(min>miny) min = miny;
		if(min>minz) min = minz;

		if(min < 0) min = -min;

		min *= 2;

		pgX = new DPG::PG(2, min, 2);

		for(int i=0; i<numPointsX; i++)
		{	
			pgX->addPoint(pointSetX[i]);
		}	

//		cout<<"Inserted all points of mesh 1 to DPG"<<endl;

		minx = 0, miny = 0, minz = 0;

		for(int i=0; i<numPointsY; i++)
		{
			if(isRAWN)
			{
				if(fscanf(fpY, "%lf %lf %lf %lf %lf %lf", &px, &py, &pz, &tempnx, &tempny, &tempnz) != 6)
				{
					cout<<"Error reading vertex "<<i<<endl;
					return false;
				}
			}
			else 
			{
				if(fscanf(fpY, "%lf %lf %lf", &px, &py, &pz) != 3)
				{
					cout<<"Error reading vertex "<<i<<endl;
					return false;
				}
			}

			DPG::Point* p = new DPG::Point(px, py, pz);
			pointSetY.push_back(p);

			if(minx>px) minx = px;
			if(miny>py) miny = py;
			if(minz>pz) minz = pz;
		}

//		cout<<"Read all points of mesh 2"<<endl;

		min = minx;
		if(min>miny) min = miny;
		if(min>minz) min = minz;

		if(min < 0) min = -min;

		min *= 2;

		pgY = new DPG::PG(2, min, 2);

		for(int i=0; i<numPointsY; i++)
		{	
			pgY->addPoint(pointSetY[i]);
		}

//		cout<<"Inserted all points of mesh 2 to DPG"<<endl;

		fclose(fpX);
		fclose(fpY);

		return true;
	}
	else if(fpX!=NULL)
	{
		cout<<"Failed to open file "<<fileName2<<endl;
		fclose(fpX);
		return false;
	}
	else if(fpY!=NULL)
	{
		cout<<"Failed to open file "<<fileName1<<endl;
		fclose(fpY);
		return false;
	}
	else
	{
		cout<<"Failed to open files "<<fileName1<<" and "<<fileName2<<endl;
		return false;
	}
}

bool HausdorffDistance::computeAllDistances()
{
	vector<DPG::Point*> closePoints;
	double queryDist = 2.0;
	double currentMinDist = 400.0;

//	cout<<"Computing distances from mesh 1"<<endl;

	for(int i=0; i<numPointsX; i++)
	{	
		currentMinDist = 400.0;
		queryDist = 1.0;
		closePoints.clear();

		DPG::Point* x = pointSetX[i];

		while(! pgY->pointsWithinRange(x, queryDist ) )
			queryDist += 1.0;

		closePoints = pgY->range(x, queryDist);

		for(int j=0; j<closePoints.size(); j++)
		{
			DPG::Point* y = closePoints[j];

			double dist = x->distsq(y);

			if(currentMinDist > dist) currentMinDist = dist;
		}

/*		for(int j=0; j<pointSetY.size(); j++)
		{
			DPG::Point* y = pointSetY[j];

			double dist = x->distsq(y);

			if(currentMinDist > dist) currentMinDist = dist;
		}
*/
		currentMinDist = sqrt(currentMinDist);

//		cout<<i<<" "<< currentMinDist <<endl;

		if(oneWayHausdorffDistXY < currentMinDist)
			oneWayHausdorffDistXY = currentMinDist;

		oneWayAverageDistXY += currentMinDist;
	}

	oneWayAverageDistXY /= (double)numPointsX;
	
//	cout<<"Computing distances from mesh 2"<<endl;

	for(int i=0; i<numPointsY; i++)
	{	
		currentMinDist = 400.0;
		queryDist = 1.0;
		closePoints.clear();

		DPG::Point* y = pointSetY[i];

		while(! pgX->pointsWithinRange(y, queryDist ) )
			queryDist += 1.0;

		closePoints = pgX->range(y, queryDist);

		for(int j=0; j<closePoints.size(); j++)
		{
			DPG::Point* x = closePoints[j];

			double dist = y->distsq(x);

			if(currentMinDist > dist) currentMinDist = dist;
		}

/*		for(int j=0; j<pointSetX.size(); j++)
		{
			DPG::Point* x = pointSetX[j];

			double dist = x->distsq(y);

			if(currentMinDist > dist) currentMinDist = dist;
		}*/

		currentMinDist = sqrt(currentMinDist);

		if(oneWayHausdorffDistYX < currentMinDist)
			oneWayHausdorffDistYX = currentMinDist;

		oneWayAverageDistYX += currentMinDist;
	}
	closePoints.clear();

	oneWayAverageDistYX /= (double)numPointsY;

	averageDist = (oneWayAverageDistYX + oneWayAverageDistXY)/2.0;

	hausdorffDist = oneWayHausdorffDistXY > oneWayHausdorffDistYX ? oneWayHausdorffDistXY : oneWayHausdorffDistYX;
}

bool HausdorffDistance::cleanup()
{
	pointSetX.clear();
	pointSetY.clear();
}

HausdorffDistance::HausdorffDistance(string fileName1, string fileName2, bool isRAWN)
{
	if(readFile(fileName1, fileName2, isRAWN))
		computeAllDistances();	
}

HausdorffDistance::~HausdorffDistance()
{
	cleanup();
}
