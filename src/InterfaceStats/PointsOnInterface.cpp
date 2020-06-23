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
#include <InterfaceStats/PointsOnInterface.h>
#include <iostream>
#include <fstream>

PointsOnInterface::PointsOnInterface( vector<DPG::WeightedPoint *>& pts1, vector<DPG::WeightedPoint *>& pts2 )
{
	setDefaults( );

	copyPoints( pts1, true );
	copyPoints( pts2, false );

	buildOctrees( );	
	buildPGs();
}


PointsOnInterface::~PointsOnInterface( )
{
	freeMemory( );
	if ( printStatus ) printf( "\n" );
}


void PointsOnInterface::freeMemory( void )
{
	staticOctree.clear();
	movingOctree.clear();  

	pointList1.clear();
	pointList2.clear();

	interfacePointList.clear();
}


void PointsOnInterface::setDefaults( void )
{
	staticOctreeBuilt = false;
	movingOctreeBuilt = false;

	staticOctreeRoot = -1;
	movingOctreeRoot = -1;

	curNode = maxNode = 0;
	curSubtreeRoot = maxSubtreeRoot = 0;

	minRadius = 4.0;
	minRadiusUsed = -1;
	maxLeafSize = 500;
	maxLeafSizeUsed = -1;
	epsilon = 0.5;
	interfaceWidth = 2.2;

	transMatrix.reset( );

	printStatus = true;
}


void PointsOnInterface::copyPoints(vector<DPG::WeightedPoint *>& pts, bool firstSet)
{
	int size = pts.size();
	
	for(int i=0; i<size; i++)
	{
		if(firstSet)
			pointList1.push_back(pts[i]);
		else
			pointList2.push_back(pts[i]);
	}
}


void PointsOnInterface::computeXlateForPG( )
{
	double minXYZ;
	DPG::Point* p;
	double x, y, z;

	double minX1, minY1, minZ1;
	double maxX1, maxY1, maxZ1;
	double minX2, minY2, minZ2;
	double maxX2, maxY2, maxZ2;
	
	double xDim1, yDim1, zDim1;
	double xDim2, yDim2, zDim2;

	double maxDim1, maxDim2;
	double minXYZ1, minXYZ2;

	int numPoints = pointList1.size();
   
	for ( int i = 0; i < numPoints; i++ )
	{
		p = pointList1[i];
	
		x = p->getX();
		y = p->getY();
		z = p->getZ();	

		if ( i == 0 )
		{ 
			minX1 = maxX1 = x;
			minY1 = maxY1 = y;
			minZ1 = maxZ1 = z;
		}

		if ( x < minX1 ) minX1 = x;
		else if ( x > maxX1 ) maxX1 = x;

		if ( y < minY1 ) minY1 = y;
		else if ( y > maxY1 ) maxY1 = y;

		if ( z < minZ1 ) minZ1 = z;
		else if ( z > maxZ1 ) maxZ1 = z;
     	}

	minXYZ1 = minX1;
	if(minXYZ1 > minY1) minXYZ1 = minY1;
	if(minXYZ1 > minZ1) minXYZ1 = minZ1;

	numPoints = pointList2.size();
   
	for ( int i = 0; i < numPoints; i++ )
	{
		p = pointList2[i];
	
		x = p->getX();
		y = p->getY();
		z = p->getZ();	

		if ( i == 0 )
		{ 
			minX2 = maxX2 = x;
			minY2 = maxY2 = y;
			minZ2 = maxZ2 = z;
		}

		if ( x < minX2 ) minX2 = x;
		else if ( x > maxX2 ) maxX2 = x;

		if ( y < minY2 ) minY2 = y;
		else if ( y > maxY2 ) maxY2 = y;

		if ( z < minZ2 ) minZ2 = z;
		else if ( z > maxZ2 ) maxZ2 = z;
     	}

	minXYZ2 = minX2;
	if(minXYZ2 > minY2) minXYZ2 = minY2;
	if(minXYZ2 > minZ2) minXYZ2 = minZ2;

	xDim1 = maxDim1 = maxX1 - minX1;
	yDim1 = maxY1 - minY1;
	zDim1 = maxZ1 - minZ1;

	if(maxDim1 < yDim1) maxDim1 = yDim1;
	if(maxDim1 < zDim1) maxDim1 = zDim1;

	xDim2 = maxDim2 = maxX2 - minX2;
	yDim2 = maxY2 - minY2;
	zDim2 = maxZ2 - minZ2;

	if(maxDim2 < yDim2) maxDim2 = yDim2;
	if(maxDim2 < zDim2) maxDim2 = zDim2;

	minXYZ1 -= 2*maxDim2;
	minXYZ2 -= 2*maxDim1;

	xlatePG1 = xlatePG2 = 0.0;

	if(minXYZ1 < 0) xlatePG1 = -minXYZ1;
	if(minXYZ2 < 0) xlatePG2 = -minXYZ2;
}


void PointsOnInterface::buildPGs( )
{
	computeXlateForPG();

	pgPoints1 = new DPG::PG( 3.0, xlatePG1, 3.0 );  
	pgPoints2 = new DPG::PG( 3.0, xlatePG2, 3.0 );  

	int numPoints = pointList1.size();
	DPG::Point* p;

	for ( int i = 0; i < numPoints; i++ )
	{
		p = pointList1[i];
		pgPoints1 -> addPoint(p);
	}	

	numPoints = pointList2.size();

	for ( int i = 0; i < numPoints; i++ )
	{
		p = pointList2[i];
		pgPoints2 -> addPoint(p);
	}
}


// Recursively creates the octrees

int PointsOnInterface::constructOctree( int pointsStartID, int pointsEndID, vector<DPG::WeightedPoint *>& ptsT, bool firstSet)
{
	POINTS_OCTREE_NODE pon;

	DPG::WeightedPoint *pt;
	if(firstSet) pt = pointList1[pointsStartID];
	else pt = pointList2[pointsStartID];

	double minX = pt->getX();
	double minY = pt->getY();
	double minZ = pt->getZ();
	double maxX = pt->getX();
	double maxY = pt->getY();
	double maxZ = pt->getZ();

	for ( int i = pointsStartID + 1; i <= pointsEndID; i++ )
 	{
		DPG::WeightedPoint *temp;
		if(firstSet) temp = pointList1[ i ];
		else temp = pointList2[ i ];

		double x = temp->getX(), y = temp->getY(), z = temp->getZ();

		if ( x < minX ) minX = x;
		if ( x > maxX ) maxX = x;

		if ( y < minY ) minY = y;
		if ( y > maxY ) maxY = y;

		if ( z < minZ ) minZ = z;
		if ( z > maxZ ) maxZ = z;
	} 

	double cx = pon.cx = ( minX + maxX ) / 2;
	double cy = pon.cy = ( minY + maxY ) / 2;
	double cz = pon.cz = ( minZ + maxZ ) / 2;

	DPG::Point center(cx,cy,cz);

	double d2 = center.distsq(pt);

	for ( int i = pointsStartID + 1; i <= pointsEndID; i++ )
	{
		DPG::WeightedPoint *temp;
		if(firstSet) temp = pointList1[ i ];
		else temp = pointList2[ i ];

		double d2T = center.distsq(temp);

 		if ( d2T > d2 ) d2 = d2T;
 	} 

	double cr = pon.cr = sqrt( d2 );

	
	if ( ( pointsEndID - pointsStartID +  1 > maxLeafSize ) && ( cr > minRadius ) )	//checking if splitting is needed
	{
		pon.leaf = false;

		int ptsCount[ 8 ] = { 0, 0, 0, 0, 0, 0, 0, 0 };

		// Finding bucket sizes for the splitcells

		for ( int i = pointsStartID; i <= pointsEndID; i++ )
		{
			if(firstSet) ptsT[ i ] = pointList1[ i ];
			else ptsT[ i ] = pointList2[ i ];

			DPG::WeightedPoint *temp;
			if(firstSet) temp = pointList1[ i ];
			else temp = pointList2[ i ];

			int splitCellIndex = ( zeroIfLess( temp->getZ(), cz ) << 2 ) + ( zeroIfLess( temp->getY(), cy ) << 1 ) + ( zeroIfLess( temp->getX(), cx ) );

			ptsCount[ splitCellIndex ]++;	
		}
  
		int ptsStartIndex[ 8 ];
		int ptsCurIndex[ 8 ];  

		ptsCurIndex[ 0 ] = ptsStartIndex[ 0 ] = pointsStartID;

		for ( int i = 1; i < 8; i++ )
			ptsCurIndex[ i ] = ptsStartIndex[ i ] = ptsStartIndex[ i - 1 ] + ptsCount[ i - 1 ];

		//Bucketing the atoms for the splitcells

		for ( int i = pointsStartID; i <= pointsEndID; i++ )	
		{  
			int splitCellIndex = ( zeroIfLess( ptsT[ i ]->getZ(), cz ) << 2 ) + ( zeroIfLess( ptsT[ i ]->getY(), cy ) << 1 ) + ( zeroIfLess( ptsT[ i ]->getX(), cx ) );
 
			if(firstSet)
				pointList1[ ptsCurIndex[ splitCellIndex ] ] = ptsT[ i ];
			else
				pointList2[ ptsCurIndex[ splitCellIndex ] ] = ptsT[ i ];

			ptsCurIndex[ splitCellIndex ]++;  
		}  

		// recurse for the splitcells
	  
		for ( int i = 0; i < 8; i++ ) 
		{
			if ( ptsCount[ i ] > 0 )
			{
				pon.cPtr[ i ] = constructOctree( ptsStartIndex[ i ], ptsStartIndex[ i ] + ptsCount[ i ] - 1, ptsT, firstSet);
			}

			else
			{
				pon.cPtr[ i ] = -1; 
			}
		}  
	}

	else
	{
		pon.leaf = true;
	}

	pon.pointsStartID = pointsStartID;
	pon.pointsEndID = pointsEndID;

	if(firstSet)
	{
		staticOctree.push_back(pon);
		return staticOctree.size()-1;	// returning index, so that parent can set child pointers
	}
	else
	{
		movingOctree.push_back(pon);
		return movingOctree.size()-1;	// returning index, so that parent can set child pointers
	}
}


bool PointsOnInterface::buildStaticOctree(  )
{  
	if ( printStatus ) printf( "\nbuilding static octree...\n " );

	int numPoints = pointList1.size();

	vector<DPG::WeightedPoint*> ptsT;

	for(int i=0; i<numPoints; i++)
		ptsT.push_back(NULL);

	staticOctreeRoot = constructOctree( 0, numPoints - 1, ptsT, true );
	staticOctreeBuilt = true;

	ptsT.clear();

	return true;
}


bool PointsOnInterface::buildMovingOctree(  )
{  
	if ( printStatus ) printf( "\nbuilding moving octree... \n" );

	int numPoints = pointList2.size();

	vector<DPG::WeightedPoint*> ptsT;

	for(int i=0; i<numPoints; i++)
		ptsT.push_back(NULL);

	movingOctreeRoot = constructOctree( 0, numPoints - 1, ptsT, false );
	movingOctreeBuilt = true;

	ptsT.clear();

	return true;
}


bool PointsOnInterface::buildOctrees( void )
{
	if ( ( minRadius != minRadiusUsed ) || ( maxLeafSize != maxLeafSizeUsed ) || ( !staticOctreeBuilt ) ) 				buildStaticOctree( );  
	if ( ( minRadius != minRadiusUsed ) || ( maxLeafSize != maxLeafSizeUsed ) || ( !movingOctreeBuilt ) ) 				buildMovingOctree( );
	minRadiusUsed = minRadius;
	maxLeafSizeUsed = maxLeafSize;
}


Matrix PointsOnInterface::computeInverse(Matrix transMat)
{
	double invR[3][3];
	double oldT[3];
	double newT[3] = {0,0,0};

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			invR[i][j] = transMat.get(j,i);

			newT[i] -= invR[i][j]*transMat.get(j,3);
		}
	}

	Matrix invMat;

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			invMat.set(i,j, invR[i][j]);
		}
		invMat.set(i,3, newT[i]);
	}

	for(int j=0;j<4;j++)
	{
		invMat.set(3,j, transMat.get(3,j));
	}

	return invMat;
}


void PointsOnInterface::approximateInteractions( Matrix transMat, Matrix transMatInv, int nodeS, int nodeM )
{
	double sumRad = staticOctree[ nodeS ].cr + movingOctree[ nodeM ].cr;
	double sumRad2 = sumRad * sumRad;

	Vector oldPos( movingOctree[ nodeM ].cx, movingOctree[ nodeM ].cy, movingOctree[ nodeM ].cz, 1.0 );
	Vector newPos = transMat * oldPos;
	double cxM = newPos[ 0 ], cyM = newPos[ 1 ], czM = newPos[ 2 ];

	double dx = staticOctree[ nodeS ].cx - cxM,
	dy = staticOctree[ nodeS ].cy - cyM,
	dz = staticOctree[ nodeS ].cz - czM;
	double d2 = dx * dx + dy * dy + dz * dz;

	bool farEnough = false;

	if ( d2 > sumRad2 )// && ( d2 > (sumRad2/(epsilon * epsilon)) )  ) 
		farEnough = true;

	if ( farEnough ) 
	{
		return;
	}  
	else
	{
		if ( staticOctree[ nodeS ].leaf && movingOctree[ nodeM ].leaf )	// Base case: both are leaves.
		{
			for ( int i = movingOctree[ nodeM ].pointsStartID; i <= movingOctree[ nodeM ].pointsEndID; i++ )
			{
				Vector oldPos( pointList2[ i ]->getX(), pointList2[ i ]->getY(), pointList2[ i ]->getZ(), 1.0 );
				Vector newPos = transMat * oldPos;
				double xM = newPos[ 0 ], yM = newPos[ 1 ], zM = newPos[ 2 ];

				DPG::WeightedPoint *temp = new DPG::WeightedPoint(xM, yM, zM, pointList2[ i ]->getWeight());
			  
				if(pgPoints1->pointsWithinRange(temp, interfaceWidth ))
				{
					if(!movingInterfacePoints[i])
					{
						interfacePointList.push_back(temp);
						movingInterfacePoints[i] = true;	
						area += temp->getWeight();
					}
				}
			}

			for ( int j = staticOctree[ nodeS ].pointsStartID; j <= staticOctree[ nodeS ].pointsEndID; j++ )
			{
				Vector oldPos( pointList1[ j ]->getX(), pointList1[ j ]->getY(), pointList1[ j ]->getZ(), 1.0 );
				Vector newPos = transMatInv * oldPos;
				double xM = newPos[ 0 ], yM = newPos[ 1 ], zM = newPos[ 2 ];

				DPG::WeightedPoint *temp = new DPG::WeightedPoint(xM, yM, zM, pointList1[ j ]->getWeight());
			  
				if(pgPoints2->pointsWithinRange(temp, interfaceWidth ))
				{
					if(!staticInterfacePoints[j])
					{
						interfacePointList.push_back(pointList1[ j ]);
						staticInterfacePoints[j] = true;	
					}
				}
			} 
		}

		else if ( !staticOctree[ nodeS ].leaf && !movingOctree[ nodeM ].leaf )
		{
			for ( int i = 0; i < 8; i++ )
			{
				if ( staticOctree[ nodeS ].cPtr[ i ] >= 0 )
				{
					for ( int j = 0; j < 8; j++ )
					{
						if ( movingOctree[ nodeM ].cPtr[ j ] >= 0 ) 
 						{
							approximateInteractions( transMat, transMatInv, staticOctree[ nodeS ].cPtr[ i ], movingOctree[ nodeM ].cPtr[ j ]);
						}  
					}
				}
			}
		}

		else if ( !staticOctree[ nodeS ].leaf )
		{
			for ( int i = 0; i < 8; i++ )
			{
				if ( staticOctree[ nodeS ].cPtr[ i ] >= 0 )
 				{
					approximateInteractions( transMat, transMatInv, staticOctree[ nodeS ].cPtr[ i ], nodeM);
 				}  
			}
 		}

		else 
		{
			for ( int j = 0; j < 8; j++ )
			{
				if ( movingOctree[ nodeM ].cPtr[ j ] >= 0 )
				{
					approximateInteractions( transMat, transMatInv, nodeS, movingOctree[ nodeM ].cPtr[ j ] );
				}  
			}
		}  
	}
}



bool PointsOnInterface::computeInteractions( Matrix transMat, InterfaceStats *intStat )
{
	if ( !staticOctreeBuilt || !movingOctreeBuilt ) return false;

	interfacePointList.clear();

	int size = pointList1.size();

	for(int i=0; i<size; i++)
		staticInterfacePoints.push_back(false);

	size = pointList2.size();

	for(int i=0; i<size; i++)
		movingInterfacePoints.push_back(false);
	

	Matrix transMatInv = computeInverse(transMat);

	area = 0;
 
	approximateInteractions( transMat, transMatInv, staticOctreeRoot, movingOctreeRoot ); // fills the interfacePointList

	cout<<"Identified "<<interfacePointList.size() << " interface points" <<endl;
	cout<<"Total points "<< pointList1.size() +  pointList2.size() <<endl;

/*	std::ofstream outfile;
	outfile.open("interfacePoints.xyz");
	for(int i=0;i<interfacePointList.size();i++)
		outfile<<interfacePointList[i]->getX()<<" "<<interfacePointList[i]->getY()<<" "<<interfacePointList[i]->getZ()<<endl;
	outfile.close();*/

	staticInterfacePoints.clear();
	movingInterfacePoints.clear();

	intStat->setInterfaceArea(area);

	return true;
}
