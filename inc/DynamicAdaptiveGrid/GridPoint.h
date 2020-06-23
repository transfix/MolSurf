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
#ifndef DAG_GRIDPOINT_H
#define DAG_GRIDPOINT_H

#include <DPG/PG.h>
#include <DynamicAdaptiveGrid/GridCell.h>
#include <PDBParser/Atom.h>

namespace DynamicAdaptiveGrid
{
	enum { 	VDW_POINT,
		SAS_POINT, 
		OUT_POINT, 
		UNDEFINED_POINT
	};

	class GridCell;

	class GridPoint: public DPG::Point
	{
		private:
			DPG::PG *exposedAtomPG;
			int defaultMode;
			double solventRad;
			static int gridPointIDGenerator;
			int gridPointID;
			bool hasVal;
			bool hasVDW;
			bool hasSES;
			bool hasGauss;
			int paircount, triplecount, pairprob, tripleprob, tripletestcount, pairtestcount;
			void generateGridPointID();
			double getPointToCircleDistance(DPG::Ball* a1, DPG::Ball* a2, double x, double y, double z, double solventRad);
			bool getSDFForIntersectingSpheres(DPG::Ball* a, DPG::Ball* b, double x, double y, double z, double solventRad, double &sdf);
			bool getSDFForTripleIntersect(DPG::Ball* a, DPG::Ball* b, DPG::Ball* c, double x, double y, double z, double solventRad, double &sdf);

		protected:
			double functionValue;
			int state;
			vector<GridCell *> gridCells;
			
		public:
			GridPoint();
			GridPoint(const GridPoint &gp);
			GridPoint(double a, double b, double c);
			GridPoint(double a, double b, double c, double F);
			GridPoint(double a, double b, double c, double F, int s);

			int getState();
			double getFunctionValue();
			double getFunctionValue(int mode);	// 0 = vdw, 1 = ses, 2 = gauss
			int getGridPointID();
			vector<GridCell *> getGridCells();
			GridCell * getGridCell(int i);
			int getNumGridCells();
			bool hasValue();

			void setState(int s);
			void setMode(int m){defaultMode = m;}
			void markAsNoVal(){hasVal = false; hasVDW = false; hasSES = false; hasGauss = false;}
			void setPG(DPG::PG *ePG){exposedAtomPG = ePG;}
			void setSolventRad(double rad){solventRad = rad;}
			void setFunctionValue(double f);
			void addToFunctionValue(double f);
			void setGridCells(vector<Point *> gcs);

			void sesAtPoint();
			void vdwAtPoint();
			void gaussAtPoint();
			void gaussAtPointSingleAtom(PDBParser::Atom* at, bool sub);

			void setPosition(double a, double b, double c);

			int getPairCount(){return paircount;}
			int getPairProb(){return pairprob;}
			int getTripleCount(){return triplecount;}
			int getTripleProb(){return tripleprob;}
			int getTripleTestCount(){return tripletestcount;}
			int getPairTestCount(){return pairtestcount;}
	};

};

#endif
