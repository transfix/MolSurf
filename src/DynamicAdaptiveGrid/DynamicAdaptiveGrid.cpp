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
#include <DynamicAdaptiveGrid/DynamicAdaptiveGrid.h>
#include <iostream>
#include <fstream>

using namespace DynamicAdaptiveGrid;

DynamicAdaptiveGrids::DynamicAdaptiveGrids(double grid_spacing, double pg_spacing, int m, double solvent_rad, double iso, double mnx, double mny, double mnz, double mxx, double mxy, double mxz)
{
	gridSpacing = grid_spacing;
	solventRadius = solvent_rad;
	dpgCellSpacing = pg_spacing;
	minx = mnx, miny = mny, minz = mnz, maxx = mxx, maxy = mxy, maxz = mxz;
	isovalue = iso;
	mode = m;
	init();
}

DynamicAdaptiveGrids::~DynamicAdaptiveGrids()
{
	cleanup();
}

/*
void DynamicAdaptiveGrids::createGrid()
{
	double x,y,z;

	gridPointArray = new GridPoint**[zCells+1];

	double centerDist = gridSpacing * 0.87;	// using slightly higher than sqrt(3)/2. It is used as query distance.

	for(int i=0; i<zCells+1; i++)
	{
		gridPointArray[i] = new GridPoint*[yCells+1];

		for(int j=0; j<yCells+1; j++)
		{
			gridPointArray[i][j] = new GridPoint[xCells+1];

			for(int k=0; k<xCells+1; k++)
			{
				x = minx + k*gridSpacing;
				y = miny + j*gridSpacing;
				z = minz + i*gridSpacing;

				gridPointArray[i][j][k] = GridPoint(x,y,z);

				gridPointPG->addPoint(&gridPointArray[i][j][k]);
			}
		}
	}


	gridCellArray = new GridCell**[zCells];

	for(int i=0; i<zCells; i++)
	{
		gridCellArray[i] = new GridCell*[yCells];

		for(int j=0; j<yCells; j++)
		{
			gridCellArray[i][j] = new GridCell[xCells];

			for(int k=0; k<xCells; k++)
			{
				x = minx + k*gridSpacing + gridSpacing/2;
				y = miny + j*gridSpacing + gridSpacing/2;
				z = minz + i*gridSpacing + gridSpacing/2;

				gridCellArray[i][j][k] = GridCell(x,y,z);

				gridCellArray[i][j][k].setVertexContainer(&vertCon);

				gridCellPG->addPoint(&gridCellArray[i][j][k]);

				vector<DPG::Point *> nei = gridPointPG->range(&gridCellArray[i][j][k], centerDist);
				gridCellArray[i][j][k].setGridPoints(nei);
			}
		}
	}


	for(int i=0; i<zCells+1; i++)
	{
		for(int j=0; j<yCells+1; j++)
		{
			for(int k=0; k<xCells+1; k++)
			{	
				gridPointArray[i][j][k].setGridCells(gridCellPG->range(&gridPointArray[i][j][k], centerDist));
			}
		}
	}

}

void DynamicAdaptiveGrids::destroyGrid()
{
	for(int i=0; i<zCells+1; i++)
	{
		for(int j=0; j<yCells+1; j++)
		{
			delete [] gridPointArray[i][j];
		}
		delete [] gridPointArray[i];
	}
	delete [] gridPointArray;

	for(int i=0; i<zCells; i++)
	{
		for(int j=0; j<yCells; j++)
		{
			delete [] gridCellArray[i][j];
		}
		delete [] gridCellArray[i];
	}
	delete [] gridCellArray;
}
*/

void DynamicAdaptiveGrids::init()
{
	xDim = maxx-minx;
	yDim = maxy-miny;
	zDim = maxz-minz;

	xCells = xDim/gridSpacing+1;
	yCells = yDim/gridSpacing+1;
	zCells = zDim/gridSpacing+1;

	double maxDim = xDim;
	if(maxDim<yDim) maxDim = yDim;
	if(maxDim<zDim) maxDim = zDim;

	dpgxlate = minx;
	if(miny<dpgxlate) dpgxlate = miny;
	if(minz<dpgxlate) dpgxlate = minz;

	if(dpgxlate<0) dpgxlate = -dpgxlate;

	dpgxlate += maxDim;

	if(dpgCellSpacing < 2*solventRadius)
		dpgCellSpacing = 2*solventRadius;


	allAtomPG = new DPG::PG(dpgCellSpacing, dpgxlate, solventRadius);
	exposedAtomPG = new DPG::PG(dpgCellSpacing, dpgxlate, solventRadius);
	exposedGridCellPG = new DPG::PG(dpgCellSpacing, dpgxlate, solventRadius);
	exposedGridPointPG = new DPG::PG(dpgCellSpacing, dpgxlate, solventRadius);
	gridPointPG = new DPG::PG(dpgCellSpacing, dpgxlate, solventRadius);
	gridCellPG = new DPG::PG(dpgCellSpacing, dpgxlate, solventRadius);

	numExposedAtoms = 0;
	numExposedCells = 0;

	gridCellMap.clear();
	gridPointMap.clear();
	newlyUpdated.clear();
//	createGrid();
}

void DynamicAdaptiveGrids::cleanup()
{
	gridCells.clear();
	gridPoints.clear();

	toBeDeleted.clear();
	toBeUpdated.clear();
	toBeCreated.clear();

	gridCellMap.clear();
	gridPointMap.clear();
	newlyUpdated.clear();
//	destroyGrid();
}

void DynamicAdaptiveGrids::updateMinMax(PDBParser::Atom* a)
{
	double ax = a->getX();
	double ay = a->getY();
	double az = a->getZ();

	if(ax > maxx) maxx = ax;
	if(ay > maxy) maxy = ay;
	if(az > maxz) maxz = az;

	if(ax < minx) minx = ax;
	if(ay < miny) miny = ay;
	if(az < minz) minz = az;

	xDim = maxx-minx;
	yDim = maxy-miny;
	zDim = maxz-minz;
}


void DynamicAdaptiveGrids::findGridPointsInBall(double x, double y, double z, double d, vector<GridPoint*>& gps)
{
	double numPoints = 2.0*d/gridSpacing + 1;

	double bboxMinX = x-d;
	double bboxMinY = y-d;
	double bboxMinZ = z-d;

	DPG::Point* p = new DPG::Point(x,y,z);
	double d2 = d*d;

	
	double bboxMinGPX = gridSpacing * (double)((int)(bboxMinX/gridSpacing));	// find the min corner of the bounding box
	if(bboxMinX<0) bboxMinGPX -= gridSpacing;
	double bboxMinGPY = gridSpacing * (double)((int)(bboxMinY/gridSpacing));
	if(bboxMinY<0) bboxMinGPY -= gridSpacing;
	double bboxMinGPZ = gridSpacing * (double)((int)(bboxMinZ/gridSpacing));
	if(bboxMinZ<0) bboxMinGPZ -= gridSpacing;

	for(int i=0; i<numPoints; i++)
	{
		for(int j=0; j<numPoints; j++)
		{
			for(int k=0; k<numPoints; k++)
			{
				GridPoint* gp = new GridPoint(bboxMinGPX+gridSpacing*i, bboxMinGPY+gridSpacing*j, bboxMinGPZ+gridSpacing*k);
		
//				cout<<bboxMinGPX+gridSpacing*i<<" "<< bboxMinGPY+gridSpacing*j<<" "<< bboxMinGPZ+gridSpacing*k<<endl;

				if(p->distsq(gp) <= d2)
				{	
					gp->setPG(exposedAtomPG);
					gp->setSolventRad(1.4);
					gp->setMode(mode);
					string str = gp->toString();

					map<string,GridPoint*>::iterator it = gridPointMap.find(str);

					if(it == gridPointMap.end())		// New point
					{
						gridPointMap[str] = gp;
						gridPoints.push_back(gp);
						gps.push_back(gp);
					}	
					else
					{
						gps.push_back((*it).second);
					}
				}
			}
		}
	}
}


bool DynamicAdaptiveGrids::updateContour()
{
	int numDestroy = 0;
	int numCreate = 0;

	for(int i=0; i<toBeDeleted.size(); i++)
	{
		toBeDeleted[i]->destroyContour();
	}
	for(int i=0; i<toBeUpdated.size(); i++)
	{
		toBeUpdated[i]->destroyContour();
		toBeUpdated[i]->computeContour(isovalue);
	}
	for(int i=0; i<toBeCreated.size(); i++)
	{
		toBeCreated[i]->computeContour(isovalue);	
	}

	toBeDeleted.clear();
	toBeUpdated.clear();
	toBeCreated.clear();
}


bool DynamicAdaptiveGrids::updateCell(GridCell *gc)
{
	double x = gc->getX();
	double y = gc->getY();
	double z = gc->getZ();
	double dist = gridSpacing/2.0;

	//according to the order in the marching cube code
	double mod[8][3] = { {-dist, -dist, -dist}, {dist, -dist, -dist}, {dist, dist, -dist}, {-dist, dist, -dist}, {-dist, -dist, dist}, {dist, -dist, dist}, {dist, dist, dist}, {-dist, dist, dist} }; 

	vector<GridPoint*> gps;
	int count = 0;
	bool old = false;

	for(int i=0; i<8; i++)
	{
		double xmod = x+mod[i][0];
		double ymod = y+mod[i][1];
		double zmod = z+mod[i][2];

		if(xmod<0.001 && xmod>-0.001) xmod = 0.0;
		if(ymod<0.001 && ymod>-0.001) ymod = 0.0;
		if(zmod<0.001 && zmod>-0.001) zmod = 0.0;

		GridPoint* temp = new GridPoint(xmod, ymod, zmod);
		temp->setPG(exposedAtomPG);
		temp->setSolventRad(solventRadius);
		temp->setMode(mode);
		string str = temp->toString();

		map<string,GridPoint*>::iterator it = gridPointMap.find(str);

		if(it == gridPointMap.end())		// New point
		{
			gridPointMap[str] = temp;
			gridPoints.push_back(temp);
			gps.push_back(temp);
//			cout<< str <<endl;
		}	
		else
		{
			gps.push_back((*it).second);
		}
	}

	gc->setGridPoints(gps);
	gc->classifyCell();
}


bool DynamicAdaptiveGrids::updateCells()
{
	toBeDeleted.clear();
	toBeUpdated.clear();
	toBeCreated.clear();

	map<string,GridCell*>::iterator cellStart = newlyUpdated.begin();
	map<string,GridCell*>::iterator cellEnd = newlyUpdated.end();

	map<string,GridCell*>::iterator cellIt;

//	cout<<"Total cells marked for update: "<<newlyUpdated.size()<<endl;


	int index = 0;
	for(cellIt=cellStart; cellIt!=cellEnd; cellIt++)
	{
		GridCell* temp = (*cellIt).second;

		int oldClass = temp->getState();
		updateCell(temp);
		int newClass = temp->getState();

		int oldContourMarking = temp->getContourState();

		if(newClass == VDW_CELL || newClass == BAND_CELL)		// Need to create/update mesh.
		{		
			if( oldClass == VDW_CELL || oldClass == BAND_CELL )	//update
			{
				temp->markForUpdate();
				toBeUpdated.push_back(temp);
			}
			else							//create
			{
				temp->markForCreate();
				exposedGridCellPG->addPoint(temp);
				numExposedCells++;
				toBeCreated.push_back(temp);
			}
		}	
		else if(oldClass == VDW_CELL  || oldClass == BAND_CELL)		// Need to delete mesh.
		{	
//			cout<<"New class is "<<newClass<<endl;
			temp->markForDelete();
			exposedGridCellPG->removePoint(temp);
			numExposedCells--;
			toBeDeleted.push_back(temp);
		}
	}

//	cout<<"Cells to be deleted: "<< toBeDeleted.size()<<endl;
//	cout<<"Cells to be updated: "<< toBeUpdated.size()<<endl;
//	cout<<"Cells to be created: "<< toBeCreated.size()<<endl;

	newlyUpdated.clear();
}


bool DynamicAdaptiveGrids::markNeighborCells(GridPoint *gp)
{
	double x = gp->getX();
	double y = gp->getY();
	double z = gp->getZ();
	double dist = gridSpacing/2.0;
	double mod[2] = {-dist, dist};

	for(int i=0; i<2; i++)
	{
		for(int j=0; j<2; j++)
		{
			for(int k=0; k<2; k++)
			{
				double xmod = x+mod[i];
				double ymod = y+mod[j];
				double zmod = z+mod[k];

				if(xmod<0.001 && xmod>-0.001) xmod = 0.0;
				if(ymod<0.001 && ymod>-0.001) ymod = 0.0;
				if(zmod<0.001 && zmod>-0.001) zmod = 0.0;

				GridCell* temp = new GridCell(xmod, ymod, zmod);
				string str = temp->toString();
				
				map<string,GridCell*>::iterator it = gridCellMap.find(str);

				if(it == gridCellMap.end())		// New cell found
				{
					temp->setVertexContainer(&vertCon);
					gridCellMap[str] = temp;
					gridCells.push_back(temp);
				}
				else
				{
					temp = (*it).second;
				}			

				if(newlyUpdated.find(str) == newlyUpdated.end())	// Need to add it to the update list
				{
					newlyUpdated[str] = temp;
				}				
			}
		}
	}
}

			
bool DynamicAdaptiveGrids::addAtom(PDBParser::Atom* a)
{
//	int pendingUpdates = newlyUpdated.size();
	double rad = a->getRadius();
	double radsq = (rad-0.2)*(rad-0.2);

	updateMinMax(a);

	//Find the gridpoints inside the SAS and VDW volume of the atom. Mark them as exposed/burried and update the neighboring cells if necessary.

	vector <GridPoint *> gridPointsInSAS;
//	gridPointsInSAS = gridPointPG->range(a, rad + solventRadius);

	findGridPointsInBall(a->getX(), a->getY(), a->getZ(), rad+solventRadius, gridPointsInSAS);

	int numGridPointsInSAS = gridPointsInSAS.size();
	int numGridPointsInVDW = 0;
	int sasToVdwPoint = 0;
	int outToSasPoint = 0;
	int outToVdwPoint = 0;
	int sasToSasPoint = 0;
	int vdwToVdwPoint = 0;

	bool isExposedAtom = false;

	for(int i=0; i<numGridPointsInSAS; i++)			
	{
		GridPoint* gp = (GridPoint*)gridPointsInSAS[i];

		int oldState = gp->getState();

		int newState = SAS_POINT;

		if(a->distsq(gp) < radsq)
			newState = VDW_POINT;

		if(oldState != VDW_POINT)
			gp -> setState(newState);

//		gp -> gaussAtPointSingleAtom(a, false);
		
		// the three cases, where an update is necessary are OUT->SAS, OUT->VDW, SAS->VDW
		// update may be needed for VDW->VDW, SAS->SAS
		// OUT->OUT,  cannot happen for gp
		// VDW->SAS cannot happen during insertion of an atom


		if(newState == VDW_POINT)
		{
			numGridPointsInVDW++;

			if(oldState == OUT_POINT)
			{
				isExposedAtom = true;
				markNeighborCells(gp);
				outToVdwPoint++;
			}

			else if(oldState == SAS_POINT)		// remove point from exposed list
			{
				exposedGridPointPG -> removePoint(gp);
				markNeighborCells(gp);
				sasToVdwPoint++;
			}

			else if(oldState == VDW_POINT)
			{
				markNeighborCells(gp);		//might need to update the mesh. but not the state
				vdwToVdwPoint++;
			}
		}
		else if(newState == SAS_POINT)		
		{
			if(oldState == OUT_POINT)		// add point to exposed list
			{
				exposedGridPointPG -> addPoint(gp);
				isExposedAtom = true;
				markNeighborCells(gp);
				outToSasPoint++;
			}
			else if(oldState == SAS_POINT)
			{
				markNeighborCells(gp);		//might need to update the mesh. but not the state
				sasToSasPoint++;
			}		
		}

	}

/*	cout<<"numGridPointsInSAS "<< numGridPointsInSAS<<endl;	
	cout<<"numGridPointsInVDW "<< numGridPointsInVDW<<endl;
	cout<<"sasToVdwPoint "<< sasToVdwPoint<<endl;
	cout<<"outToVdwPoint "<< outToVdwPoint<<endl;
	cout<<"vdwToVdwPoint "<< vdwToVdwPoint<<endl;
	cout<<"sasToSasPoint "<< sasToSasPoint<<endl;
	cout<<"outToSasPoint "<< outToSasPoint<<endl;*/


	//If this atom contributed to the surface (converted a GridPoint), then mark it as exposed. Check if it burried any previously exposed atom.

	if(isExposedAtom)	
	{
		vector <DPG::Point *> neighborExposedAtoms = exposedAtomPG->range(a, rad + solventRadius);

		int numNeighbors = neighborExposedAtoms.size();

		for(int i=0; i<numNeighbors; i++)
		{
			PDBParser::Atom* na = (PDBParser::Atom*)neighborExposedAtoms[i];

			if(!exposedGridPointPG->pointsWithinRange(na, rad + solventRadius))
			{
				exposedAtomPG->removePoint(na);
				numExposedAtoms--;
				na->setBoundaryFlag(false);
			}
		}

		a -> setBoundaryFlag(true);
		exposedAtomPG -> addPoint(a);
		numExposedAtoms++;
	}

	allAtomPG->addPoint(a);

//	updateCells();
//	updateContour();

//	cout<<"numExposedCells "<< numExposedCells<<endl;
//	cout<<"numExposedAtoms "<< numExposedAtoms<<endl;

}

bool DynamicAdaptiveGrids::removeAtom(PDBParser::Atom* a)
{
	int pendingUpdates = newlyUpdated.size();
	double rad = a->getRadius();
	double radsq = (rad-0.2)*(rad-0.2);

	allAtomPG->removePoint(a);

	if(a->isBoundaryAtom())
	{
		exposedAtomPG->removePoint(a);
		numExposedAtoms--;
	}

	vector <GridPoint *> gridPointsInSAS;
//	gridPointsInSAS = gridPointPG->range(a, rad + solventRadius);

	findGridPointsInBall(a->getX(), a->getY(), a->getZ(), rad+solventRadius, gridPointsInSAS);

	int numGridPointsInSAS = gridPointsInSAS.size();

	for(int i=0; i<numGridPointsInSAS; i++)			
	{
		GridPoint* gp = (GridPoint*)gridPointsInSAS[i];
		gp -> gaussAtPointSingleAtom(a, true);

		int oldState = gp->getState();

		if(oldState == SAS_POINT)
		{
			if(!allAtomPG->pointsWithinRange(gp, 2*solventRadius))	//not inside SAS of other atoms
			{
				gp->setState(OUT_POINT);
				exposedGridPointPG->removePoint(gp);
				markNeighborCells(gp);
			}
			else	// inside the SAS of other atom
			{
				markNeighborCells(gp);			
			}
		}

		else if(oldState == VDW_POINT)
		{
			if(!allAtomPG->pointsWithinRange(gp, solventRadius))  //  it is no longer VDW
			{
				vector <DPG::Point *> newlyExposedAtom = allAtomPG->range(gp, 2*solventRadius);
				int numNewlyExposedAtom = newlyExposedAtom.size();

				for(int j=0; j<numNewlyExposedAtom; j++)
				{
					PDBParser::Atom * na = (PDBParser::Atom*)newlyExposedAtom[j];

					if(!na->isBoundaryAtom())
					{
//						cout<<"   0 "<<j<<endl;
						na->setBoundaryFlag(true);
						exposedAtomPG->addPoint(na);
						numExposedAtoms++;
					}
				}

				if(numNewlyExposedAtom>0)	//SAS
				{
					gp->setState(SAS_POINT);
					exposedGridPointPG->addPoint(gp);
					markNeighborCells(gp);
				}
				
				else				//OUT
				{
					gp->setState(OUT_POINT);
					markNeighborCells(gp);
				}
			}
		}
	}

//	updateCells();
//	updateContour();

//	cout<<"numExposedCells "<< numExposedCells<<endl;
//	cout<<"numExposedAtoms "<< numExposedAtoms<<endl;

}

bool DynamicAdaptiveGrids::moveAtom(PDBParser::Atom* from, PDBParser::Atom* to)
{
	removeAtom(from);
	addAtom(to);
}


bool DynamicAdaptiveGrids::addAtoms(vector<PDBParser::Atom*> &atms)
{
	int numAtoms = atms.size();

	for(int i=0;i<numAtoms;i++)
		addAtom(atms[i]);
//	cout<<"added atoms and updated cell values"<<endl;
	updateCells();
	updateContour();
}

bool DynamicAdaptiveGrids::removeAtoms(vector<PDBParser::Atom*> &atms)
{
	int numAtoms = atms.size();

	for(int i=0;i<numAtoms;i++)
		removeAtom(atms[i]);

//	cout<<"removed atoms and updated cell values"<<endl;
	updateCells();
	updateContour();
}


bool DynamicAdaptiveGrids::printSurface(string surfaceFileName)
{
	double cx = (maxx+minx)/2;
	double cy = (maxy+miny)/2;
	double cz = (maxz+minz)/2;

	double d = xDim;
	if(d<yDim) d = yDim;
	if(d<zDim) d = zDim;

	DPG::Point* p = new DPG::Point(cx,cy,cz);

	vector<DPG::Point *> cells = exposedGridCellPG->range(p, d);

	int size = cells.size();

// Printing the exposed cells
/*
	cout<<"Found exposed cells "<<size<<endl;

	string cellFileName = "cell" + surfaceFileName + ".xyzr";
	ofstream cellFile(cellFileName.c_str());

	for(int i=0; i<size; i++)
	{
		cellFile<< cells[i]->getX() <<" "<< cells[i]->getY() <<" "<< cells[i]->getZ() <<" 0.05"<<endl;
	}
	cellFile.close();


	cout<<"Found cells "<<gridCellMap.size()<<endl;

	string cellFileName2 = "allcell" + surfaceFileName + ".xyzr";
	ofstream cellFile2(cellFileName2.c_str());

	map<string,GridCell*>::iterator cellStart = gridCellMap.begin();
	map<string,GridCell*>::iterator cellEnd = gridCellMap.end();

	map<string,GridCell*>::iterator cellIt;
	for ( cellIt=cellStart ; cellIt != cellEnd; cellIt++ )
	{
		GridCell* cell = (*cellIt).second;
		cellFile2<< cell->getX() <<" "<< cell->getY() <<" "<< cell->getZ() <<" 0.05"<<endl;
	}
	cellFile2.close();
*/
// Printing the surface

	string surfFileName = surfaceFileName;
	ofstream surfFile(surfFileName.c_str());

	int numVerts = vertCon.size();
	int numTraingles = 0;

	for(int i=0; i<size; i++)
	{
		numTraingles += ((GridCell *)cells[i])->getNumSimplices();		
	}

	surfFile << numVerts << " "<< numTraingles << endl;		//printing the number of vertices and triangles 

	vertCon.printVerticesAndUpdateMap(&surfFile);			//prints the mesh vertices and updates a map from vertex->index in the raw file

	for(int i=0; i<size; i++)
	{
		((GridCell *)cells[i])->printMesh(&surfFile);		//prints the mesh in the cell	
	}
	surfFile.close();

// Printing the exposed gridpoints

	size = gridPointMap.size();
	cout<<"Found points "<<size<<endl;
	string vdwPointFileName = "vdw" + surfaceFileName + ".xyzr";
	string sasPointFileName = "sas" + surfaceFileName + ".xyzr";
	string outPointFileName = "out" + surfaceFileName + ".xyzr";
	ofstream vdwPointFile(vdwPointFileName.c_str());
	ofstream sasPointFile(sasPointFileName.c_str());
	ofstream outPointFile(outPointFileName.c_str());

	double vdwVal = 0;
	int vdwNum = 0;
	
	double sasVal = 0;
	int sasNum = 0;

	map<string,GridPoint*>::iterator pointStart = gridPointMap.begin();
	map<string,GridPoint*>::iterator pointEnd = gridPointMap.end();

	map<string,GridPoint*>::iterator pointIt;
	int index = 0;

	int paircount = 0;
	int triplecount = 0;
	int pairprob = 0;
	int tripleprob = 0;
	int pairtest = 0;	
	int tripletest = 0;

	for ( pointIt=pointStart ; pointIt != pointEnd; pointIt++ )
	{
		GridPoint* temp = (*pointIt).second;

		paircount += temp->getPairCount();
		triplecount += temp->getTripleCount(); 
		pairprob += temp->getPairProb(); 
		tripleprob += temp->getTripleProb();
		pairtest += temp->getPairTestCount();
		tripletest += temp->getTripleTestCount();
		
		double val = temp->getFunctionValue();

		val/= 10.0;

		if(val>=0.000)
		{
			if(temp->getState() == VDW_POINT)
			{
				vdwPointFile<< temp->getX() <<" "<< temp->getY() <<" "<< temp->getZ() <<" "<<0.05<<endl;
				vdwVal += temp->getFunctionValue();
				vdwNum++;
			}
			else if(temp->getState() == SAS_POINT)
			{
				sasPointFile<< temp->getX() <<" "<< temp->getY() <<" "<< temp->getZ() <<" "<<0.05<<endl;
				sasVal += temp->getFunctionValue();
				sasNum++;
			}
			else if(temp->getState() == OUT_POINT)
			{
				outPointFile<< temp->getX() <<" "<< temp->getY() <<" "<< temp->getZ() <<" "<<0.15<<endl;
			}
		}
		index++;
	}

	cout<<"Total points: "<<index<<endl;
	cout<<"Average value of vdw points: "<<vdwVal/(double)vdwNum<<endl;
	cout<<"Average value of sas points: "<<sasVal/(double)sasNum<<endl;

	cout<<"Pairs checked: "<<pairtest<<endl;
	cout<<"Pairs used: "<<paircount<<endl;
	cout<<"Pair test errors: "<<pairprob<<endl;

	cout<<"Triples checked: "<<tripletest<<endl;
	cout<<"Triples used: "<<triplecount<<endl;
	cout<<"Triple test errors: "<<tripleprob<<endl;

	vdwPointFile.close();
	sasPointFile.close();
	outPointFile.close();

	
}
