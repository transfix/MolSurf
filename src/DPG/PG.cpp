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
#include <DPG/PG.h>
#include <cmath>
#include <time.h>
#include <sys/time.h>

using namespace std;
using namespace DPG;

//#define TRANSLATE 2000.0

double gtod_sec = 0.0E0;

double gtod_timer()
{
   struct timeval tv;
   struct timezone Tzp;
   double sec;

   gettimeofday(&tv, &Tzp);

   if(gtod_sec == 0.0E0) 
      gtod_sec = (double)tv.tv_sec;
   sec = (double)tv.tv_sec - gtod_sec;

   return sec + 1.0E-06*(double)tv.tv_usec;
}


bool PG::doesOverlap(Ball *q, double delta) 
{
	vector<Point*> result;
	range(q, delta + q->getRadius(), result);
	int i;
	int size = result.size();
	double r;
	for(i=0;i<size;i++)
	{
		r = ((Ball*)result[i])->getRadius() + q->getRadius(); 
		if(q->distsq(result[i]) < r*r)
			return true;
	}
	return false;
}


bool PG::range(Point *q, double delta, vector<Point*>& result) 
{
	double x, y, z;
	x = q->getX() + TRANSLATE;
	y = q->getY() + TRANSLATE;
	z = q->getZ() + TRANSLATE;

	int cx, cy, cz;			// map the query point to a cell
	cx = (int) x/DIM;
	cy = (int) y/DIM;
	cz = (int) z/DIM;

	int d = (int) delta/DIM;	// query distance in terms of number of cells
	int dSq = d*d; 

	int lx, hx, ly, hy, lz, hz;	// used for defining the integer range. l for low, h for high
 
	double delsq = delta*delta;  

	lx = cx - d - 1;
	hx = cx + d + 1;
	ly = cy - d - 1;
	hy = cy + d + 1;  
	lz = cz - d - 1;
	hz = cz + d + 1;
	vector<tuple<plane*> > planesInRange = g.RR->report(lz,hz);		// searching planes in grid

	if(planesInRange.empty())
	{
		return false;
	}

	else									// some planes found in the range. 
	{
		int numPlanes = planesInRange.size();

		vector<tuple<line*> >  linesInRange;

		for(int i = 0; i < numPlanes; i++) 				// for each plane, look for lines in them within the range
		{
/*			int cellDistZ = cz - (planesInRange[i].ptr)->ID.z;
			int cellDistZSq = cellDistZ*cellDistZ;

			int reducedQueryDistYSq = dSq - cellDistZSq;		// in planes further from the point, a smaller range need to checked
			int reducedQueryDistY = sqrt(reducedQueryDistYSq);	// due to the query being spherical

			ly = cy - reducedQueryDistY - 1;			
			hy = cy + reducedQueryDistY + 1;	*/		

			linesInRange = (planesInRange[i].ptr)->RR->report(ly,hy);

			if(linesInRange.empty())
				continue;

			else
			{
    				int numLines = (int) linesInRange.size();

				vector<tuple<gridcell*> > cellsInRange;

    				for(int j = 0; j < numLines; j++) 
    				{
/*					int cellDistY = cy - (linesInRange[j].ptr)->ID.y;
      					int cellDistYSq = cellDistY*cellDistY;

					int reducedQueryDistXSq = reducedQueryDistYSq - cellDistYSq;	// in lines further from the point, a smaller range need to checked
					int reducedQueryDistX = sqrt(reducedQueryDistXSq);		// due to the query being spherical. 

					lx = cx - reducedQueryDistX - 1;			
					hx = cx + reducedQueryDistX + 1;*/	

					cellsInRange = (linesInRange[j].ptr)->RR->report(lx,hx);
        
					if(cellsInRange.empty())
						continue;

					else
					{
		    				int numCells = (int) cellsInRange.size();

						for(int k=0; k<numCells; k++)
						{
							int atomsInCell = cellsInRange[k].ptr->balls.size();

							for(int l=0; l<atomsInCell; l++)
							{
								Point *oa = cellsInRange[k].ptr->balls[l];

								if(oa->distsq(q) <= delsq)
									result.push_back(oa);
							}
						}
	
					      	cellsInRange.clear();
      					}

    				} //end for j

    				linesInRange.clear();
  			}
		} // end for i
	}
	return true;
}


vector<Point*> PG::range(Point *q, double delta) 
{
	double x, y, z;
	x = q->getX() + TRANSLATE;
	y = q->getY() + TRANSLATE;
	z = q->getZ() + TRANSLATE;

	int cx, cy, cz;			// map the query point to a cell
	cx = (int) x/DIM;
	cy = (int) y/DIM;
	cz = (int) z/DIM;

	int d = (int) delta/DIM;	// query distance in terms of number of cells
	int dSq = d*d; 

	int lx, hx, ly, hy, lz, hz;	// used for defining the integer range. l for low, h for high
 
	vector<Point*> result;		// the resulting list of points to return

	double delsq = delta*delta;  

	lx = cx - d - 1;
	hx = cx + d + 1;
	ly = cy - d - 1;
	hy = cy + d + 1;  
	lz = cz - d - 1;
	hz = cz + d + 1;
	vector<tuple<plane*> > planesInRange = g.RR->report(lz,hz);		// searching planes in grid

	if(planesInRange.empty())
	{
		return result;
	}

	else									// some planes found in the range. 
	{
		int numPlanes = planesInRange.size();

		vector<tuple<line*> >  linesInRange;

		for(int i = 0; i < numPlanes; i++) 				// for each plane, look for lines in them within the range
		{
/*			int cellDistZ = cz - (planesInRange[i].ptr)->ID.z;
			int cellDistZSq = cellDistZ*cellDistZ;

			int reducedQueryDistYSq = dSq - cellDistZSq;		// in planes further from the point, a smaller range need to checked
			int reducedQueryDistY = sqrt(reducedQueryDistYSq);	// due to the query being spherical

			ly = cy - reducedQueryDistY - 1;			
			hy = cy + reducedQueryDistY + 1;	*/		

			linesInRange = (planesInRange[i].ptr)->RR->report(ly,hy);

			if(linesInRange.empty())
				continue;

			else
			{
    				int numLines = (int) linesInRange.size();

				vector<tuple<gridcell*> > cellsInRange;

    				for(int j = 0; j < numLines; j++) 
    				{
/*					int cellDistY = cy - (linesInRange[j].ptr)->ID.y;
      					int cellDistYSq = cellDistY*cellDistY;

					int reducedQueryDistXSq = reducedQueryDistYSq - cellDistYSq;	// in lines further from the point, a smaller range need to checked
					int reducedQueryDistX = sqrt(reducedQueryDistXSq);		// due to the query being spherical. 

					lx = cx - reducedQueryDistX - 1;			
					hx = cx + reducedQueryDistX + 1;*/	

					cellsInRange = (linesInRange[j].ptr)->RR->report(lx,hx);
        
					if(cellsInRange.empty())
						continue;

					else
					{
		    				int numCells = (int) cellsInRange.size();

						for(int k=0; k<numCells; k++)
						{
							int atomsInCell = cellsInRange[k].ptr->balls.size();

							for(int l=0; l<atomsInCell; l++)
							{
								Point *oa = cellsInRange[k].ptr->balls[l];

								if(oa->distsq(q) <= delsq)
									result.push_back(oa);
							}
						}
	
					      	cellsInRange.clear();
      					}

    				} //end for j

    				linesInRange.clear();
  			}
		} // end for i
	}
	return result;
}


bool PG::pointsWithinRange(Point *q, double delta) 
{
	double x, y, z;
	x = q->getX() + TRANSLATE;
	y = q->getY() + TRANSLATE;
	z = q->getZ() + TRANSLATE;

	int cx, cy, cz;			// map the query point to a cell
	cx = (int) x/DIM;
	cy = (int) y/DIM;
	cz = (int) z/DIM;

	int d = (int) delta/DIM;	// query distance in terms of number of cells
	int dSq = d*d; 

	int lx, hx, ly, hy, lz, hz;	// used for defining the integer range. l for low, h for high
 
	double delsq = delta*delta;  

	lx = cx - d - 1;
	hx = cx + d + 1;
	ly = cy - d - 1;
	hy = cy + d + 1;  
	lz = cz - d - 1;
	hz = cz + d + 1;
	vector<tuple<plane*> > planesInRange = g.RR->report(lz,hz);		// searching planes in grid

	if(planesInRange.empty())
	{
		return false;
	}

	else									// some planes found in the range. 
	{
		int numPlanes = planesInRange.size();

		vector<tuple<line*> >  linesInRange;

		for(int i = 0; i < numPlanes; i++) 				// for each plane, look for lines in them within the range
		{
/*			int cellDistZ = cz - (planesInRange[i].ptr)->ID.z;

			int cellDistZSq = cellDistZ*cellDistZ;

			int reducedQueryDistYSq = dSq - cellDistZSq;		// in planes further from the point, a smaller range need to checked
			int reducedQueryDistY = sqrt(reducedQueryDistYSq);	// due to the query being spherical


			ly = cy - reducedQueryDistY - 1;			
			hy = cy + reducedQueryDistY + 1;	*/		

			linesInRange = (planesInRange[i].ptr)->RR->report(ly,hy);

			if(linesInRange.empty())
				continue;

			else
			{
    				int numLines = (int) linesInRange.size();

				vector<tuple<gridcell*> > cellsInRange;

    				for(int j = 0; j < numLines; j++) 
    				{
/*					int cellDistY = cy - (linesInRange[j].ptr)->ID.y;
      					int cellDistYSq = cellDistY*cellDistY;


					int reducedQueryDistXSq = reducedQueryDistYSq - cellDistYSq;	// in lines further from the point, a smaller range need to checked
					int reducedQueryDistX = sqrt(reducedQueryDistXSq);		// due to the query being spherical. 

					lx = cx - reducedQueryDistX - 1;			

					hx = cx + reducedQueryDistX + 1;*/	

					cellsInRange = (linesInRange[j].ptr)->RR->report(lx,hx);
        
					if(cellsInRange.empty())
						continue;

					else
					{
		    				int numCells = (int) cellsInRange.size();

						for(int k=0; k<numCells; k++)
						{
							int atomsInCell = cellsInRange[k].ptr->balls.size();

							for(int l=0; l<atomsInCell; l++)
							{
								Point *oa = cellsInRange[k].ptr->balls[l];

								if(oa->distsq(q) <= delsq)
									return true;
							}
						}
	
					      	cellsInRange.clear();
      					}

    				} //end for j

    				linesInRange.clear();
  			}
		} // end for i
	}
	return false;
}


bool PG::lightRange(Point *q, double delta, vector<Point*>& result) 
{
	double x, y, z;
	x = q->getX() + TRANSLATE;
	y = q->getY() + TRANSLATE;
	z = q->getZ() + TRANSLATE;

	int cx, cy, cz;			// map the query point to a cell
	cx = (int) x/DIM;
	cy = (int) y/DIM;
	cz = (int) z/DIM;

	int d = (int) delta/DIM;	// query distance in terms of number of cells
	int dSq = d*d; 

	int lx, hx, ly, hy, lz, hz;	// used for defining the integer range. l for low, h for high

	double delsq = delta*delta;  

	lx = cx - d - 1;
	hx = cx + d + 1;
	ly = cy - d - 1;
	hy = cy + d + 1;  
	lz = cz - d - 1;
	hz = cz + d + 1;
	vector<tuple<plane*> > planesInRange = g.RR->report(lz,hz);		// searching planes in grid

	if(planesInRange.empty())
	{
		return false;
	}

	else									// some planes found in the range. 
	{
		int numPlanes = planesInRange.size();

		vector<tuple<line*> >  linesInRange;

		for(int i = 0; i < numPlanes; i++) 				// for each plane, look for lines in them within the range
		{
/*			int cellDistZ = cz - (planesInRange[i].ptr)->ID.z;
			int cellDistZSq = cellDistZ*cellDistZ;

			int reducedQueryDistYSq = dSq - cellDistZSq;		// in planes further from the point, a smaller range need to checked
			int reducedQueryDistY = sqrt(reducedQueryDistYSq);	// due to the query being spherical

			ly = cy - reducedQueryDistY - 1;			
			hy = cy + reducedQueryDistY + 1;	*/		

			linesInRange = (planesInRange[i].ptr)->RR->report(ly,hy);

			if(linesInRange.empty())
				continue;

			else
			{
    				int numLines = (int) linesInRange.size();

				vector<tuple<gridcell*> > cellsInRange;

    				for(int j = 0; j < numLines; j++) 
    				{
/*					int cellDistY = cy - (linesInRange[j].ptr)->ID.y;
      					int cellDistYSq = cellDistY*cellDistY;

					int reducedQueryDistXSq = reducedQueryDistYSq - cellDistYSq;	// in lines further from the point, a smaller range need to checked
					int reducedQueryDistX = sqrt(reducedQueryDistXSq);		// due to the query being spherical. 

					lx = cx - reducedQueryDistX - 1;			
					hx = cx + reducedQueryDistX + 1;*/	

					cellsInRange = (linesInRange[j].ptr)->RR->report(lx,hx);
        
					if(cellsInRange.empty())
						continue;

					else
					{
		    				int numCells = (int) cellsInRange.size();

						for(int k=0; k<numCells; k++)
						{
							int atomsInCell = cellsInRange[k].ptr->balls.size();

							for(int l=0; l<atomsInCell; l++)
							{
								Point *oa = cellsInRange[k].ptr->balls[l];

								result.push_back(oa);
							}
						}
	
					      	cellsInRange.clear();
      					}

    				} //end for j

    				linesInRange.clear();
  			}
		} // end for i
	}
	return true;
}


vector<Point*> PG::lightRange(Point *q, double delta) 
{
	double x, y, z;
	x = q->getX() + TRANSLATE;
	y = q->getY() + TRANSLATE;
	z = q->getZ() + TRANSLATE;

	int cx, cy, cz;			// map the query point to a cell
	cx = (int) x/DIM;
	cy = (int) y/DIM;
	cz = (int) z/DIM;

	vector<Point*> result;

	int d = (int) delta/DIM;	// query distance in terms of number of cells
	int dSq = d*d; 

	int lx, hx, ly, hy, lz, hz;	// used for defining the integer range. l for low, h for high

	double delsq = delta*delta;  

	lx = cx - d - 1;
	hx = cx + d + 1;
	ly = cy - d - 1;
	hy = cy + d + 1;  
	lz = cz - d - 1;
	hz = cz + d + 1;
	vector<tuple<plane*> > planesInRange = g.RR->report(lz,hz);		// searching planes in grid

	if(planesInRange.empty())
	{
		return result;
	}

	else									// some planes found in the range. 
	{
		int numPlanes = planesInRange.size();

		vector<tuple<line*> >  linesInRange;

		for(int i = 0; i < numPlanes; i++) 				// for each plane, look for lines in them within the range
		{
/*			int cellDistZ = cz - (planesInRange[i].ptr)->ID.z;
			int cellDistZSq = cellDistZ*cellDistZ;

			int reducedQueryDistYSq = dSq - cellDistZSq;		// in planes further from the point, a smaller range need to checked
			int reducedQueryDistY = sqrt(reducedQueryDistYSq);	// due to the query being spherical


			ly = cy - reducedQueryDistY - 1;			
			hy = cy + reducedQueryDistY + 1;	*/		

			linesInRange = (planesInRange[i].ptr)->RR->report(ly,hy);

			if(linesInRange.empty())
				continue;

			else
			{
    				int numLines = (int) linesInRange.size();

				vector<tuple<gridcell*> > cellsInRange;

    				for(int j = 0; j < numLines; j++) 
    				{
/*					int cellDistY = cy - (linesInRange[j].ptr)->ID.y;
      					int cellDistYSq = cellDistY*cellDistY;

					int reducedQueryDistXSq = reducedQueryDistYSq - cellDistYSq;	// in lines further from the point, a smaller range need to checked

					int reducedQueryDistX = sqrt(reducedQueryDistXSq);		// due to the query being spherical. 

					lx = cx - reducedQueryDistX - 1;			
					hx = cx + reducedQueryDistX + 1;*/	

					cellsInRange = (linesInRange[j].ptr)->RR->report(lx,hx);
        
					if(cellsInRange.empty())
						continue;

					else
					{
		    				int numCells = (int) cellsInRange.size();

						for(int k=0; k<numCells; k++)
						{
							int atomsInCell = cellsInRange[k].ptr->balls.size();

							for(int l=0; l<atomsInCell; l++)
							{
								Point *oa = cellsInRange[k].ptr->balls[l];

								result.push_back(oa);
							}
						}
	
					      	cellsInRange.clear();
      					}

    				} //end for j

    				linesInRange.clear();
  			}
		} // end for i
	}
	return result;
}


/*adds a ball to the collection of balls*/
void PG::addPoint(Point *a) 
{
	double x, y, z;
	x = a->getX() + TRANSLATE;
	y = a->getY() + TRANSLATE;
	z = a->getZ() + TRANSLATE;

	int cx, cy, cz;			// map the query point to a cell
	cx = (int) x/DIM;
	cy = (int) y/DIM;
	cz = (int) z/DIM;

	vector<tuple<plane*> > planesInRange = g.RR->report(cz,cz);	//identify plane
  
	if(planesInRange.empty()) 
	{
		plane* p = new plane(cz);
		tuple<plane*> tp(cz, p);
		planesInRange.push_back(tp);
		g.RR->insert(cz, p);
	}
  
	vector<tuple<line*> > linesInRange = planesInRange[0].ptr->RR->report(cy,cy);	//identify line
  
	if(linesInRange.empty()) 
	{
		line* l = new line(cy, cz); 
		tuple<line*> tl(cy, l);
		linesInRange.push_back(tl);
		planesInRange[0].ptr->RR->insert(cy, l);
	}
  
	vector<tuple<gridcell*> > cellsInRange = linesInRange[0].ptr->RR->report(cx,cx);	//identify cell
  
	if(cellsInRange.empty()) 
	{
		gridcell* c = new gridcell(cx, cy, cz);
		tuple<gridcell*> tc(cx, c);
		cellsInRange.push_back(tc);
		linesInRange[0].ptr->RR->insert(cx, c);
	}

	(cellsInRange[0].ptr->balls).push_back(a);					// add point

	return;
}


/*removes the given ball from the collection of balls*/
void PG::removePoint(Point *a) 
{
	double x, y, z;
	x = a->getX() + TRANSLATE;
	y = a->getY() + TRANSLATE;
	z = a->getZ() + TRANSLATE;

	int cx, cy, cz;			// map the query point to a cell
	cx = (int) x/DIM;
	cy = (int) y/DIM;
	cz = (int) z/DIM;
  
  	vector<tuple<plane*> > planesInRange = g.RR->report(cz,cz);

	if(planesInRange.empty()) 
	{
		cout<<"Atom does not exist"<<endl;
		return;
	}

  	vector<tuple<line*> > linesInRange = planesInRange[0].ptr->RR->report(cy,cy);
  	if(linesInRange.empty()) 
	{
		cout<<"Atom does not exist"<<endl;
		return;
	}

	vector<tuple<gridcell*> > cellsInRange = linesInRange[0].ptr->RR->report(cx,cx); 
  	if(cellsInRange.empty()) 
	{
    		cout<<"Atom does not exist"<<endl;
		return;
	}

	int index=-1;

	// cout<<"Given Point = "<<a->getX()<<" "<<a->getY()<<" "<<a->getZ()<<" "<<endl;

	for(int i=0; i<(int)(cellsInRange[0].ptr->balls).size(); i++) 
	{
		if(a->getX() == cellsInRange[0].ptr->balls[i]->getX() && a->getY() == cellsInRange[0].ptr->balls[i]->getY() && a->getZ() == cellsInRange[0].ptr->balls[i]->getZ()) //found
		{
			index = i;
			break;
		}
	}

	if(index==-1)
	{
//		cout<<"ball not found"<<endl;
	}
	else
	{
	  	cellsInRange[0].ptr->balls[index] = cellsInRange[0].ptr->balls[cellsInRange[0].ptr->balls.size()-1];
		cellsInRange[0].ptr->balls.pop_back();
		//cout<<"removed "<<endl;
	}

	return;
}
