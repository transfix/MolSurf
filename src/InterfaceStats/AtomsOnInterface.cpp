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
#include <InterfaceStats/AtomsOnInterface.h>

AtomsOnInterface::AtomsOnInterface( vector<PDBParser::Atom *>& molecule1, vector<PDBParser::Atom *>& molecule2 )
{
	setDefaults( );
	initStats( );

	copyAtoms( molecule1, true );
	copyAtoms( molecule2, false );
 
	buildOctrees( );
}


AtomsOnInterface::~AtomsOnInterface( )
{
	freeMemory( );
	if ( printStatus ) printf( "\n" );
}


void AtomsOnInterface::freeMemory( void )
{
	staticAtomsOctree.clear();
	movingAtomsOctree.clear();  
	freeMem( movingAtomsSubtreeRoots ); 
	movingAtoms.clear();
	staticAtoms.clear();
	staticInterfaceAtoms.clear();
	movingInterfaceAtoms.clear();
	staticResidueIndexMap.clear();
	movingResidueIndexMap.clear();
	staticResidueIDMap.clear();
	movingResidueIDMap.clear();
	pairList.clear();
	residueMap.clear();
}


void AtomsOnInterface::setDefaults( void )
{
	staticAtomsOctreeBuilt = false;
	movingAtomsOctreeBuilt = false;
	movingAtomsSubtreeRoots = NULL;
	staticAtomsOctreeRoot = -1;
	movingAtomsOctreeRoot = -1;
	curNode = maxNode = 0;
	curSubtreeRoot = maxSubtreeRoot = 0;

	minRadius = 2.0;
	minRadiusUsed = -1;
	maxLeafSize = 10;
	maxLeafSizeUsed = -1;
	epsilon = 0.5;
	interfaceWidth = 10.0;

	priorComputationCleared = true;
	transMatrix.reset( );

	printStatus = true;

	string residueNames[20] = {"ILE", "VAL", "LEU", "PHE", "CYS", "MET", "ALA", "GLY", "THR", "SER", "TRP", "TYR", "PRO", "HIS", "GLU", "GLN", "ASP", "ASN", "LYS", "ARG"};

	for(int i=0; i<20; i++)
	{
		isPolar[i] = false;
		residueMap[residueNames[i]] = i;
	}

	isPolar[0] = true;
	isPolar[1] = true;
	isPolar[2] = true;
	isPolar[3] = true;
	isPolar[5] = true;
	isPolar[6] = true;
	isPolar[7] = true;
	isPolar[12] = true;
}


void AtomsOnInterface::initStats()
{
	for(int i = 0; i<2; i++)
	{
		numAtoms[i] = 0;
		numPolar[i] = 0;
		numNonPolar[i] = 0;
		numCore[i] = 0;
		numRim[i] = 0;
	}

	for(int i=0; i<20; i++)
		for(int j=0; j<20; j++)
			resresCont[i][j] = 0;

	localDensityIndex = 0.0;
	globalDensityIndex = 0.0;
}


void AtomsOnInterface::copyAtoms(vector<PDBParser::Atom *>& atoms, bool stc)
{
	int size = atoms.size();
	
	for(int i=0; i<size; i++)
	{
		if(stc)
			staticAtoms.push_back(atoms[i]);
		else
			movingAtoms.push_back(atoms[i]);
	}
}


// Recursively creates the octrees

int AtomsOnInterface::constructAtomsOctree( int atomsStartID, int atomsEndID, vector<PDBParser::Atom *>& atomsT, bool stc )
{
	ATOMS_OCTREE_NODE aon;

	PDBParser::Atom *atm;
	if(stc) atm = staticAtoms[ atomsStartID ];
	else atm = movingAtoms[ atomsStartID ];

	double minX = atm->getX();
	double minY = atm->getY();
	double minZ = atm->getZ();
	double maxX = atm->getX();
	double maxY = atm->getY();
	double maxZ = atm->getZ();
	double maxR = atm->getRadius();

	for ( int i = atomsStartID + 1; i <= atomsEndID; i++ )
 	{
		PDBParser::Atom *temp;
		if(stc) temp = staticAtoms[ i ];
		else temp = movingAtoms[ i ];

		double x = temp->getX(), y = temp->getY(), z = temp->getZ(), r = temp->getRadius();

		if ( x < minX ) minX = x;
		if ( x > maxX ) maxX = x;

		if ( y < minY ) minY = y;
		if ( y > maxY ) maxY = y;

		if ( z < minZ ) minZ = z;
		if ( z > maxZ ) maxZ = z;

		if ( r > maxR ) maxR = r;
	} 

	double cx = aon.cx = ( minX + maxX ) / 2;
	double cy = aon.cy = ( minY + maxY ) / 2;
	double cz = aon.cz = ( minZ + maxZ ) / 2;

	DPG::Point center(cx,cy,cz);

	double d2 = center.distsq(atm);

	for ( int i = atomsStartID + 1; i <= atomsEndID; i++ )
	{
		PDBParser::Atom *temp;
		if(stc) temp = staticAtoms[ i ];
		else temp = movingAtoms[ i ];

		double d2T = center.distsq(temp);

 		if ( d2T > d2 ) d2 = d2T;
 	} 

	double cr = aon.cr = sqrt( d2 ) + maxR;

	if ( ( atomsEndID - atomsStartID +  1 > maxLeafSize ) && ( cr > minRadius ) )		//checking if splitting is needed
	{
		aon.leaf = false;

		int atomsCount[ 8 ] = { 0, 0, 0, 0, 0, 0, 0, 0 };

		// Finding bucket sizes for the splitcells

		for ( int i = atomsStartID; i <= atomsEndID; i++ )
		{
			if(stc) atomsT[ i ] = staticAtoms[ i ];
			else atomsT[ i ] = movingAtoms[ i ];

			PDBParser::Atom *temp;
			if(stc) temp = staticAtoms[ i ];
			else temp = movingAtoms[ i ];

			int splitCellIndex = ( zeroIfLess( temp->getZ(), cz ) << 2 ) + ( zeroIfLess( temp->getY(), cy ) << 1 ) + ( zeroIfLess( temp->getX(), cx ) );

			atomsCount[ splitCellIndex ]++;	
		}

  
		int atomsStartIndex[ 8 ];
		int atomsCurIndex[ 8 ];  

		atomsCurIndex[ 0 ] = atomsStartIndex[ 0 ] = atomsStartID;

		for ( int i = 1; i < 8; i++ )
			atomsCurIndex[ i ] = atomsStartIndex[ i ] = atomsStartIndex[ i - 1 ] + atomsCount[ i - 1 ];


		//Bucketing the atoms for the splitcells

		for ( int i = atomsStartID; i <= atomsEndID; i++ )	
		{  
			int splitCellIndex = ( zeroIfLess( atomsT[ i ]->getZ(), cz ) << 2 ) + ( zeroIfLess( atomsT[ i ]->getY(), cy ) << 1 ) + ( zeroIfLess( atomsT[ i ]->getX(), cx ) );
 
			if(stc)
				staticAtoms[ atomsCurIndex[ splitCellIndex ] ] = atomsT[ i ];
			else
				movingAtoms[ atomsCurIndex[ splitCellIndex ] ] = atomsT[ i ];

			atomsCurIndex[ splitCellIndex ]++;  
		}  

		// recurse for the splitcells
	  
		for ( int i = 0; i < 8; i++ ) 
		{
			if ( atomsCount[ i ] > 0 )
			{
				aon.cPtr[ i ] = constructAtomsOctree( atomsStartIndex[ i ], atomsStartIndex[ i ] + atomsCount[ i ] - 1, atomsT, stc );
			}

			else
			{
				aon.cPtr[ i ] = -1; 
			}
		}  
	}

	else
	{
		aon.leaf = true;
	}

	aon.atomsStartID = atomsStartID;
	aon.atomsEndID = atomsEndID;

	if(stc)
	{
		staticAtomsOctree.push_back(aon);
		return staticAtomsOctree.size()-1;	// returning index, so that parent can set child pointers
	}
	else
	{
		movingAtomsOctree.push_back(aon);
		return movingAtomsOctree.size()-1;	// returning index, so that parent can set child pointers
	}
}



void AtomsOnInterface::fillMovingAtomsSubtreeRootsArray( int nodeID, int maxNodesInLevel )
{
	if ( nodeID < 0 ) 
		return;

	if ( nodeID == movingAtomsOctreeRoot ) 
		maxSubtreeRoot = 0;

	if ( movingAtomsOctree[ nodeID ].leaf || ( maxNodesInLevel >= movingAtomsSubtreeRootsSize ) ) 
	{
		movingAtomsSubtreeRoots[ maxSubtreeRoot++ ] = nodeID;
	}
	else 
	{
		for ( int i = 0; i < 8; i++ ) 
			if ( movingAtomsOctree[ nodeID ].cPtr[ i ] >= 0 )
				fillMovingAtomsSubtreeRootsArray( movingAtomsOctree[ nodeID ].cPtr[ i ], 8 * maxNodesInLevel );
	}
}


bool AtomsOnInterface::buildStaticAtomsOctree(  )
{  
	if ( printStatus ) printf( "\nbuilding moving atoms octree... " );

	int size = staticAtoms.size();

	vector<PDBParser::Atom*> atomsT;

	for(int i=0; i<size; i++)
		atomsT.push_back(staticAtoms[i]);

	staticAtomsOctreeRoot = constructAtomsOctree( 0, size - 1, atomsT, true );

	staticAtomsOctreeBuilt = true;

	atomsT.clear();

	for(int i=0; i<size; i++)
		staticInterfaceAtoms.push_back(false);

	return true;
}



bool AtomsOnInterface::buildMovingAtomsOctree(  )
{  
	if ( printStatus ) printf( "\nbuilding moving atoms octree... " );

	int size = movingAtoms.size();

	vector<PDBParser::Atom*> atomsT;

	for(int i=0; i<size; i++)
		atomsT.push_back(movingAtoms[i]);

	movingAtomsOctreeRoot = constructAtomsOctree( 0, size - 1, atomsT, false );

	movingAtomsOctreeBuilt = true;

	atomsT.clear();

	for(int i=0; i<size; i++)
		movingInterfaceAtoms.push_back(false);

//	fillMovingAtomsSubtreeRootsArray( movingAtomsOctreeRoot, 1 );

	return true;
}


bool AtomsOnInterface::buildOctrees( void )
{
	if ( ( minRadius != minRadiusUsed ) || ( maxLeafSize != maxLeafSizeUsed ) || ( !staticAtomsOctreeBuilt ) ) 				buildStaticAtomsOctree( );  
	if ( ( minRadius != minRadiusUsed ) || ( maxLeafSize != maxLeafSizeUsed ) || ( !movingAtomsOctreeBuilt ) ) 				buildMovingAtomsOctree( );
	minRadiusUsed = minRadius;
	maxLeafSizeUsed = maxLeafSize;
}


void AtomsOnInterface::identifyRimCore()
{
/*	map<int,int>::iterator resStart = movingResidueIDMap.begin();
	map<int,int>::iterator resEnd = movingResidueIDMap.end();

	map<int,int>::iterator resIt;

	for ( resIt=resStart ; resIt!= resEnd; resIt++ )
	{
		int numInterfaceAtom = (*resIt).second;

		int resID = (*resIt).first;
		int resIndex = movingResidueIndexMap[resID];
		int numAllAtoms = numResidueAtoms[resIndex];
	}	*/

	numCore[0] = 0;
	numCore[1] = 0;
	numRim[0] = 0;
	numRim[1] = 0;
}


void AtomsOnInterface::updateStats(int staticIndex, int movingIndex )
{
	PDBParser::Atom *sa = staticAtoms[staticIndex];
	PDBParser::Atom *ma = movingAtoms[movingIndex];
	bool newMovingAtomFound = false;
	bool newStaticAtomFound = false;
	bool newMovingResidueFound = false;
	bool newStaticResidueFound = false;
	bool isStaticHydrogen = false;
	bool isMovingHydrogen = false;

	int movingAtomID = ma->getAtomNum();
	int staticAtomID = sa->getAtomNum();
	int movingResidueID = ma->getResidueNum();
	int staticResidueID = sa->getResidueNum();
	int pairID = movingResidueID*10000+staticResidueID;

	string movingResidueName = string(ma->getResidueName());
	string staticResidueName = string(sa->getResidueName());	
	int staticResidueIndex = residueMap[staticResidueName];
	int movingResidueIndex = residueMap[movingResidueName];

	string movingAtomName = string(ma->getName());
	string staticAtomName = string(sa->getName());

	if(staticAtomName[0] == 'H')
		isStaticHydrogen = true;
	if(movingAtomName[0] == 'H')
		isMovingHydrogen = true;

	if(!movingInterfaceAtoms[movingIndex])	
	{
		newMovingAtomFound = true;		
		movingInterfaceAtoms[movingIndex] = true;	
	}	

	if(!staticInterfaceAtoms[staticIndex])	
	{
		newStaticAtomFound = true;		
		staticInterfaceAtoms[staticIndex] = true;	
	}


	if(movingResidueIDMap.find(movingResidueID)==movingResidueIDMap.end())
	{
		newMovingResidueFound = true;
		movingResidueIndexMap[movingResidueID] = residueMap[movingResidueName];

		if(newMovingAtomFound) 
		{
			if(isMovingHydrogen)	
				movingResidueIDMap[movingResidueID] = 0;
			else
				movingResidueIDMap[movingResidueID] = 1;
		}

		if(isPolar[movingResidueIndex])
			numPolar[1]++;
		else 
			numNonPolar[1]++;
	}
	else
	{
		newMovingResidueFound = false;

		int oldval = movingResidueIDMap[movingResidueID];
		if(newMovingAtomFound && !isMovingHydrogen) 
			movingResidueIDMap[movingResidueID] = oldval+1;
	}

	if(staticResidueIDMap.find(staticResidueID)==staticResidueIDMap.end())
	{
		newStaticResidueFound = true;
		staticResidueIndexMap[staticResidueID] = residueMap[staticResidueName];

		if(newStaticAtomFound) 
		{
			if(isStaticHydrogen)	
				staticResidueIDMap[staticResidueID] = 0;
			else
				staticResidueIDMap[staticResidueID] = 1;
		}

		if(isPolar[staticResidueIndex])
			numPolar[0]++;
		else 
			numNonPolar[0]++;
	}
	else 
	{
		newStaticResidueFound = false;

		int oldval = staticResidueIDMap[staticResidueID];
		if(newStaticAtomFound && !isStaticHydrogen) 
			staticResidueIDMap[staticResidueID] = oldval+1;
	}

	if(newStaticResidueFound || newMovingResidueFound)	//update resresCont
	{
		pairList.push_back(pairID);
		resresCont[staticResidueIndex][movingResidueIndex]++;
	}
}


void AtomsOnInterface::approximateInteractions( Matrix transMat, int nodeS, int nodeM )
{
	double sumRad = staticAtomsOctree[ nodeS ].cr + movingAtomsOctree[ nodeM ].cr + interfaceWidth;
	double sumRad2 = sumRad * sumRad;

	Vector oldPos( movingAtomsOctree[ nodeM ].cx, movingAtomsOctree[ nodeM ].cy, movingAtomsOctree[ nodeM ].cz, 1.0 );
	Vector newPos = transMat * oldPos;
	double cxM = newPos[ 0 ], cyM = newPos[ 1 ], czM = newPos[ 2 ];

	double dx = staticAtomsOctree[ nodeS ].cx - cxM,
	dy = staticAtomsOctree[ nodeS ].cy - cyM,
	dz = staticAtomsOctree[ nodeS ].cz - czM;
	double d2 = dx * dx + dy * dy + dz * dz;

	double neighborDistCutoff;
	bool farEnough = false;

	if ( d2 > sumRad2 ) 
		farEnough = true;

	if ( farEnough ) 
	{
		return;
	}  
	else
	{
		if ( staticAtomsOctree[ nodeS ].leaf && movingAtomsOctree[ nodeM ].leaf )	// Base case: both are leaves.
		{
			for ( int i = movingAtomsOctree[ nodeM ].atomsStartID; i <= movingAtomsOctree[ nodeM ].atomsEndID; i++ )
			{
				Vector oldPos( movingAtoms[ i ]->getX(), movingAtoms[ i ]->getY(), movingAtoms[ i ]->getZ(), 1.0 );
				Vector newPos = transMat * oldPos;
				double xM = newPos[ 0 ], yM = newPos[ 1 ], zM = newPos[ 2 ];
				DPG::Point* tempM = new DPG::Point(xM,yM,zM);
  
				for ( int j = staticAtomsOctree[ nodeS ].atomsStartID; j <= staticAtomsOctree[ nodeS ].atomsEndID; j++ )
				{
					dx = staticAtoms[ j ]->getX() - xM;
					dy = staticAtoms[ j ]->getY() - yM;
					dz = staticAtoms[ j ]->getZ() - zM;

					d2 = dx * dx + dy * dy + dz * dz; 

					neighborDistCutoff = interfaceWidth;

					if( d2 < neighborDistCutoff*neighborDistCutoff )	// these atoms are on the interface
					{
						if(!movingInterfaceAtoms[i])	
						{
							numAtoms[1]++;
							interfaceAtoms.push_back(staticAtoms[j]);
						}	

						if(!staticInterfaceAtoms[j])	
						{
							numAtoms[0]++;
							interfaceAtoms.push_back(tempM);	// transformed point
						}	

						updateStats(j, i);
					}
				}
 			}
		}

		else if ( !staticAtomsOctree[ nodeS ].leaf && !movingAtomsOctree[ nodeM ].leaf )
		{
			for ( int i = 0; i < 8; i++ )
			{
				if ( staticAtomsOctree[ nodeS ].cPtr[ i ] >= 0 )
				{
					for ( int j = 0; j < 8; j++ )
					{
						if ( movingAtomsOctree[ nodeM ].cPtr[ j ] >= 0 ) 
 						{
							approximateInteractions( transMat, staticAtomsOctree[ nodeS ].cPtr[ i ], movingAtomsOctree[ nodeM ].cPtr[ j ]);
						}  
					}
				}
			}
		}

		else if ( !staticAtomsOctree[ nodeS ].leaf )
		{
			for ( int i = 0; i < 8; i++ )
			{
				if ( staticAtomsOctree[ nodeS ].cPtr[ i ] >= 0 )
 				{
					approximateInteractions( transMat, staticAtomsOctree[ nodeS ].cPtr[ i ], nodeM);
 				}  
			}
 		}

		else 
		{
			for ( int j = 0; j < 8; j++ )
			{
				if ( movingAtomsOctree[ nodeM ].cPtr[ j ] >= 0 )
				{
					approximateInteractions( transMat, nodeS, movingAtomsOctree[ nodeM ].cPtr[ j ] );
				}  
			}
		}  
	}
}



bool AtomsOnInterface::computeInteractions( Matrix transMat, InterfaceStats *intStat )
{
	if ( !staticAtomsOctreeBuilt || !movingAtomsOctreeBuilt ) return false;

	int size = staticAtoms.size();

	for(int i=0; i<size; i++)
		staticInterfaceAtoms.push_back(false);

	size = movingAtoms.size();

	for(int i=0; i<size; i++)
		movingInterfaceAtoms.push_back(false);


	staticResidueIndexMap.clear();
	movingResidueIndexMap.clear();
	staticResidueIDMap.clear();
	movingResidueIDMap.clear();
	pairList.clear();
	initStats();
 
	cout<<"Before computation"<<endl;
	cout<<numAtoms[0]<<" "<< numCore[0]<< " "<< numRim[0]<< " "<< numPolar[0]<< " "<< numNonPolar[0]<< " "<<endl;

	approximateInteractions( transMat, staticAtomsOctreeRoot, movingAtomsOctreeRoot );

	identifyRimCore();

/*	std::ofstream outfile;
	outfile.open("interfaceAtoms.xyz");
	for(int i=0;i<interfaceAtoms.size();i++)
		outfile<<interfaceAtoms[i]->getX()<<" "<<interfaceAtoms[i]->getY()<<" "<<interfaceAtoms[i]->getZ()<<endl;
	outfile.close();*/

	intStat->setLocalDensityIndex(localDensityIndex);
	intStat->setGlobalDensityIndex(globalDensityIndex);

	cout<<"Matrix: \n";
	transMat.print();
	cout<<"Found "<< numAtoms[0] << " and "<< numAtoms[1] << " interface atoms"<<endl;
	cout<<"Total atoms "<< staticAtoms.size() << " and "<< movingAtoms.size()<<endl;

	for(int i=0; i<20; i++)
		for(int j=0; j<20; j++)
			intStat->setResresContact(i,j,resresCont[i][j]);

	cout<<"After computation"<<endl;
	cout<<numAtoms[0]<<" "<< numCore[0]<< " "<< numRim[0]<< " "<< numPolar[0]<< " "<< numNonPolar[0]<< " "<<endl;

	intStat->setNumAtoms(numAtoms[0],0);
	intStat->setNumCoreResidue(numCore[0],0);
	intStat->setNumRimResidue(numRim[0],0);
	intStat->setNumPolarResidue(numPolar[0],0);
	intStat->setNumNonPolarResidue(numNonPolar[0],0);

	intStat->setNumAtoms(numAtoms[1],1);
	intStat->setNumCoreResidue(numCore[1],1);
	intStat->setNumRimResidue(numRim[1],1);
	intStat->setNumPolarResidue(numPolar[1],1);
	intStat->setNumNonPolarResidue(numNonPolar[1],1);

	return lsp.compute(interfaceAtoms, intStat);
}


vector<int> AtomsOnInterface::getLigandInterfaceAtoms()
{
	vector<int> ligandInterfaceAtoms;

	for(int i=0;i<movingAtoms.size();i++)
		if(movingInterfaceAtoms[i])
			ligandInterfaceAtoms.push_back(i);

	return ligandInterfaceAtoms;
}
