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
#ifndef CCV_BINDING_SITE_SUMMARIZER_H
#define CCV_BINDING_SITE_SUMMARIZER_H

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
};

class BindingSiteSummarizer
{
	public:
		BindingSiteSummarizer(string receptorFName, string ligandFName, string outputFName, double interfaceWidth, bool verbose = false);	
		// Surface File 1 = receptor (on which the footprint is computed)
		// Surface File 2 = ligand
		// xforms are applied to ligand

		virtual ~BindingSiteSummarizer();

		bool computeSiteSummary(vector<Matrix>& xforms);
		bool computeSiteSummary(vector<Matrix>& xforms, vector<double>& scores);

	private:
		bool verbose;

		string receptorSurfaceFileName;
		string ligandSurfaceFileName;

		string bindingSitesFileName;

		double interfaceWidth;
		Matrix xform;

		int numXforms;

		double maxTimesOnSite;
		double minTimesOnSite;
	
		int numVertR, numTriangR;
		int numVertL, numTriangL;

		bool useScore;
		double score;

		DPG::PG *pgL;

		vector<DPG::WeightedPoint*> receptorPoints;
		vector<DPG::Point*> ligandPoints;
		vector<DPG::Point*> ligandPointsXformed;
		vector<triangle> receptorTriangles;
		vector<triangle> ligandTriangles;

		bool filesRead;


		void cleanup();
		bool initialize();
		bool detectInterfacePoints(Matrix xform);
		bool prepareBindingSiteSurfaceFile();
};


#endif
