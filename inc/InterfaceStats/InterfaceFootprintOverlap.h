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
#ifndef CCV_INTERFACE_FOOTPRINT_OVERLAP_H
#define CCV_INTERFACE_FOOTPRINT_OVERLAP_H

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <Utility/utility.h>
#include <UsefulMath/Vector.h>
#include <UsefulMath/Matrix.h>
#include <DPG/PG.h>

using namespace std;
using CCVOpenGLMath::Matrix;
using CCVOpenGLMath::Vector;

struct triangle
{
	int vertex[3];
	double area;
};

class InterfaceFootprintOverlap
{
	public:
		InterfaceFootprintOverlap(string SurfaceFileName1, string SurfaceFileName2, string SurfaceFileName3, double interfaceWidth, bool verbose = false);	
		// Surface File 1 = receptor (on which the footprint is computed)
		// Surface File 2 = ligand 1
		// Surface File 3 = ligand 2
		// xform is applied to ligand 2
		// if genSurfMode = true, then it prepares new raw files where the interfaces and overlaps are highlighted. Otherwise simply computes the areas.
		// poseId = id/rank of the pose corresponding to the transformation (it is only used to name a file, use dummy value if the rank is not available)

		virtual ~InterfaceFootprintOverlap();

		bool computeInterfaces(Matrix xform);

		bool prepareStaticInterfaceSurfaceFiles();		
		bool prepareMovingInterfaceSurfaceFiles(int poseId);		

		double getAreaOfReceptorInterfaceWithStaticLigand() {return staticInterfaceComputed ? areaOfReceptorInterfaceWithLigand1 : 0.0;}
		double getAreaOfReceptorInterfaceWithMovingLigand() {return movingInterfaceComputed ? areaOfReceptorInterfaceWithLigand2 : 0.0;}
		double getAreaOfStaticLigandInterface() {return staticInterfaceComputed ? areaOfLigandInterface1 : 0.0;}
		double getAreaOfMovingLigandInterface() {return movingInterfaceComputed ? areaOfLigandInterface2 : 0.0;}
		double getAreaOfOverlap() {return movingInterfaceComputed ? areaOfOverlap : 0.0;}

	private:
		bool verbose;

		string receptorSurfaceFileName;
		string ligandSurfaceFileName1;
		string ligandSurfaceFileName2;

		string receptorHighlightedSurfaceFileName1;
		string ligandHighlightedSurfaceFileName1;

		string receptorHighlightedSurfaceFileName2;
		string ligandHighlightedSurfaceFileName2;

		string receptorFootprintSurfaceFileName;

		double interfaceWidth;
		Matrix xform;

		double areaOfReceptorInterfaceWithLigand1;
		double areaOfReceptorInterfaceWithLigand2;
		double areaOfLigandInterface1;	
		double areaOfLigandInterface2;	

		double areaOfOverlap;

		FILE *receptorSurfaceFile, *ligandSurfaceFile1, *ligandSurfaceFile2;
		FILE *receptorHighlightedSurface1, *receptorHighlightedSurface2;
		FILE *ligandHighlightedSurface1, *ligandHighlightedSurface2;
		FILE *receptorFootprintSurface;

		int numVertR, numTriangR;
		int numVertL1, numTriangL1;
		int numVertL2, numTriangL2;

		DPG::PG *pgR, *pgL1, *pgL2;

		vector<DPG::WeightedPoint*> receptorPoints;
		vector<DPG::WeightedPoint*> staticLigandPoints;
		vector<DPG::WeightedPoint*> movingLigandPoints;
		vector<DPG::WeightedPoint*> movingLigandPointsXformed;

		vector<triangle> receptorTriangles;
		vector<triangle> staticLigandTriangles;
		vector<triangle> movingLigandTriangles;

		bool* receptorPointsOnInterfaceWithLigand1;
		bool* receptorPointsOnInterfaceWithLigand2;
		bool* ligandPointsOnInterface1;
		bool* ligandPointsOnInterface2;

		bool* receptorPointsOnFootPrintOverlap;

		bool filesRead;
		bool staticInterfaceComputed;
		bool movingInterfaceComputed;
		


		void cleanup();

		bool openFiles();
		bool readFiles();

		bool initializeStaticData();
		bool detectInterfacePointsForStaticLigand();
		bool initialize();

		bool initializeMovingData();
		bool detectInterfacePointsForMovingLigand();
};


#endif
