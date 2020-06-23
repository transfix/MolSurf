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
#include<InterfaceStats/BindingSiteSummarizer.h>
#include <algorithm>


BindingSiteSummarizer::BindingSiteSummarizer(string sf1, string sf2, string sf3, double iw, bool v)
{
	verbose = v;

	receptorSurfaceFileName = sf1;
	ligandSurfaceFileName = sf2;
	bindingSitesFileName = sf3;
	interfaceWidth = iw;

	filesRead = false;
	initialize();
}

BindingSiteSummarizer::~BindingSiteSummarizer()
{
	cleanup();
}


void BindingSiteSummarizer::cleanup()
{
	receptorPoints.clear();
	ligandPoints.clear();
	ligandPointsXformed.clear();

	receptorTriangles.clear();
	ligandTriangles.clear();

	if(pgL) delete pgL;
}


bool BindingSiteSummarizer::initialize()
{
	filesRead = false;

	FILE* receptorSurfaceFile = fopen(receptorSurfaceFileName.c_str(), "rt");
	if(!receptorSurfaceFile)
	{
		cout<<"could not open surface file "<< receptorSurfaceFileName.c_str() << endl;
		return false;
	}

	FILE* ligandSurfaceFile = fopen(ligandSurfaceFileName.c_str(), "rt");
	if(!ligandSurfaceFile)
	{
		cout<<"could not open surface file "<< ligandSurfaceFileName.c_str() << endl;
		fclose(receptorSurfaceFile);
		return false;
	}


	double x, y, z, nx, ny, nz;
	int index1, index2, index3;
	double dx, dy, dz;
	double a, b, c;
	double s;
	double area;

	if(fscanf(receptorSurfaceFile, "%d %d", &numVertR, &numTriangR) != 2)
	{
		cout<<"Wrong format in receptor surface file\n";
		return false;
	}

	if(verbose) cout<<"Receptor has: "<< numVertR << " vertices and " << numTriangR << " triangles\n";


	for(int i=0; i<numVertR; i++)
	{
		if(fscanf(receptorSurfaceFile, "%lf %lf %lf %lf %lf %lf", &x, &y, &z, &nx, &ny, &nz) != 6)
		{
			cout<<"Wrong format in receptor surface file. At vertex no: "<< i <<"\n";
			cout<< "x = " << x << ", y = " << y << ", z = " << z << ", nx = " << nx << ", ny = " << ny << ", nz = " << nz << "\n";
			return false;
		}

		DPG::WeightedPoint *p = new DPG::WeightedPoint(x, y, z, 0.0);
		receptorPoints.push_back(p);
	}

	if(verbose) cout<<"Points read from receptor surface and DPG created\n";

	for(int i=0; i<numTriangR; i++)				
	{
		if(fscanf(receptorSurfaceFile, "%d %d %d", &index1, &index2, &index3) != 3)
		{
			cout<<"Wrong format in receptor surface file\n";
			return false;
		}

		triangle tr;
		tr.vertex[0] = index1;
		tr.vertex[1] = index2;
		tr.vertex[2] = index3;

		receptorTriangles.push_back(tr);		
	}

	fclose(receptorSurfaceFile);

	if(verbose) cout<<"Triangles read from receptor surface\n";


	if(fscanf(ligandSurfaceFile, "%d %d", &numVertL, &numTriangL) != 2)
	{
		cout<<"Wrong format in ligand surface file\n";
		return false;
	}

	if(verbose) cout<<"Ligand as: "<< numVertL << " vertices and " << numTriangL << " triangles\n";


	for(int i=0; i<numVertL; i++)
	{
		if(fscanf(ligandSurfaceFile, "%lf %lf %lf %lf %lf %lf", &x, &y, &z, &nx, &ny, &nz) != 6)
		{
			cout<<"Wrong format in ligand surface file. At vertex no: "<< i <<"\n";
			cout<< "x = " << x << ", y = " << y << ", z = " << z << ", nx = " << nx << ", ny = " << ny << ", nz = " << nz << "\n";
			return false;
		}

		DPG::Point *p = new DPG::Point(x, y, z);
		ligandPoints.push_back(p);
	}

	if(verbose) cout<<"Points read from ligand surface\n";

	fclose(ligandSurfaceFile);

	return true;
}


bool BindingSiteSummarizer::detectInterfacePoints(Matrix xform)
{
	pgL = new DPG::PG(10.0, 1000.0, 3.0);
	ligandPointsXformed.clear();

	double x, y, z;
	double newx, newy, newz;

	for(int i=0; i<numVertL; i++)
	{
		x = ligandPoints[i]->getX();
		y = ligandPoints[i]->getY();
		z = ligandPoints[i]->getZ();

		newx = xform.get(0,0)*x + xform.get(0,1)*y + xform.get(0,2)*z + xform.get(0,3);
		newy = xform.get(1,0)*x + xform.get(1,1)*y + xform.get(1,2)*z + xform.get(1,3);
		newz = xform.get(2,0)*x + xform.get(2,1)*y + xform.get(2,2)*z + xform.get(2,3);

		DPG::Point *p = new DPG::Point(newx, newy, newz);
		pgL->addPoint(p);
		ligandPointsXformed.push_back(p);
	}

	if(verbose) cout<<"DPG for ligand created\n";

	for(int i=0; i<numVertR; i++)
	{
		if(pgL->pointsWithinRange(receptorPoints[i], interfaceWidth))
		{
			if(!useScore)
				receptorPoints[i]->setWeight( receptorPoints[i]->getWeight() + 1 );
			else
				receptorPoints[i]->setWeight( receptorPoints[i]->getWeight() + score );
		}
	}

	delete pgL;

	if(verbose) cout<<"Receptor's interface with moving computed\n";

	return true;
}


bool BindingSiteSummarizer::computeSiteSummary(vector<Matrix>& xforms)
{
	useScore = false;
	
	numXforms = xforms.size();

	for(int i=0; i<numXforms; i++)
	{
		if(verbose) cout<<"Computing interface for pose "<<i+1<<endl;
		detectInterfacePoints(xforms[i]);
	}
	
	maxTimesOnSite = 0;
	minTimesOnSite = numXforms;
	double w;
	for(int i=0; i<numVertR; i++)
	{
		w = receptorPoints[i]->getWeight();
		if(w > maxTimesOnSite) maxTimesOnSite = w;
		else if(w < minTimesOnSite) minTimesOnSite = w;
	}

	return prepareBindingSiteSurfaceFile();
}


bool BindingSiteSummarizer::computeSiteSummary(vector<Matrix>& xforms, vector<double>& scores)
{
	useScore = true;
	
	numXforms = xforms.size();

	double maxscore = scores[0], minscore = scores[0];
	for(int i=1; i<numXforms; i++)
	{
		if(scores[i] > maxscore) maxscore = scores[i];
		if(scores[i] < minscore) minscore = scores[i];
	}

	if(maxscore>0) 
		maxTimesOnSite = numXforms*maxscore;
	else
		maxTimesOnSite = 0.0;
	
	if(minscore<0) 
		minTimesOnSite = numXforms*minscore;
	else
		minTimesOnSite = 0.0;

	for(int i=0; i<numXforms; i++)
	{
		if(verbose) cout<<"Computing interface for pose "<<i+1<<endl;
		if(verbose) cout<<"Score "<<scores[i]<<endl;
		score = scores[i];
		detectInterfacePoints(xforms[i]);
	}
	
	double w;
	for(int i=0; i<numVertR; i++)
	{
		w = receptorPoints[i]->getWeight();
		if(w > maxTimesOnSite) maxTimesOnSite = w;
		else if(w < minTimesOnSite) minTimesOnSite = w;
	}

	if(verbose)
	{
		cout<<"Max score: "<<maxTimesOnSite<<endl;
		cout<<"Min score: "<<minTimesOnSite<<endl;
	}

	return prepareBindingSiteSurfaceFile();
}


bool BindingSiteSummarizer::prepareBindingSiteSurfaceFile()
{
	FILE *bindingSitesFile = fopen(bindingSitesFileName.c_str(), "wt");
	if(!bindingSitesFile)
	{
		cout<<"could not create surface file "<< bindingSitesFileName.c_str() << endl;
		return false;
	}

	// preparing overlap surface for receptor
	fprintf(bindingSitesFile, "%d %d\n", numVertR, numTriangR);

	double mean, spanLow, spanHigh;

	if(useScore)	// use median
	{
		vector<double> weights;
		for(int i=0; i<numVertR; i++)
		{
			weights.push_back(receptorPoints[i]->getWeight());
		}
		std::sort(weights.begin(), weights.end());

		mean = weights[numVertR/2];
		
		spanLow = mean - minTimesOnSite;
		spanHigh = maxTimesOnSite - mean;
	}

	else	// use mean
	{
		spanLow = spanHigh = (maxTimesOnSite-minTimesOnSite)/2.0;
		mean = (maxTimesOnSite+minTimesOnSite)/2.0;
	}
	
	double r,g,b;
	double sub;
	double w;

	if(verbose) cout<<"Mean: "<<mean<<endl;

	for(int i=0; i<numVertR; i++)
	{
		fprintf(bindingSitesFile, "%lf %lf %lf", receptorPoints[i]->getX(), receptorPoints[i]->getY(), receptorPoints[i]->getZ());

		w = receptorPoints[i]->getWeight();

		if(w > mean)
		{
			sub = (w - mean)/spanHigh;
			fprintf(bindingSitesFile, " %lf %lf %lf\n", 1.0-sub, 1.0, 1.0-sub);		
		}
		else
		{
			sub = (mean - w)/spanLow;
			fprintf(bindingSitesFile, " %lf %lf %lf\n", 1.0, 1.0-sub, 1.0-sub);		
		}
	}

	for(int i=0; i<numTriangR; i++)
	{
		fprintf(bindingSitesFile, "%d %d %d\n", receptorTriangles[i].vertex[0], receptorTriangles[i].vertex[1], receptorTriangles[i].vertex[2]);
	}

	
	fclose(bindingSitesFile);

	return true;
}

