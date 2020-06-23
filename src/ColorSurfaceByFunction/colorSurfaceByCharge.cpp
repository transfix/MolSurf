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
#include <ColorSurfaceByFunction/colorSurfaceByCharge.h>
#include <PDBParser/FlattenGOA.h>
#include <DPG/PG.h>
#include <string>
#include <sstream>


ColorSurfaceByCharge::ColorSurfaceByCharge(string pqrf, string sf, string csf, int m, double dist, double bl)
{
	pqrFileName = pqrf;
	surfaceFileName = sf;
	coloredSurfaceFileName = csf;
	mode = m;
	distanceCutoff = dist;
	blobbyness = bl;
}

ColorSurfaceByCharge::~ColorSurfaceByCharge()
{

}

bool ColorSurfaceByCharge::colorify()
{

/* Reading the atoms and inserting into dpg */

	PDBParser::GroupOfAtoms* molecule = 0;
	GOALoader* gLoader = new GOALoader();
	molecule = gLoader->loadFile(pqrFileName.c_str());
	delete gLoader;

	if(!molecule)
	{
		cout<<"could not open pqr file "<< surfaceFileName.c_str() << endl;
		return false;
	}

	vector<PDBParser::Atom*> m_AtomList;

	PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType;
	if(!PDBParser::GroupOfAtoms::intToRadiusType(&radiusType, 0))
	{
		return false;
	}
	PDBParser::CollectionData* collectionData = 0;
	if(molecule->type == PDBParser::COLLECTION_TYPE)
	{
		collectionData = molecule->m_CollectionData;
	}
	FlattenGOA(molecule, m_AtomList, collectionData, 0, 0, 0, radiusType, PDBParser::ATOM_TYPE, false);

	DPG::PG *pg1 = new DPG::PG(10.0, 1000.0, 3.0);

	int numAtoms = m_AtomList.size();
	DPG::PG *pg = new DPG::PG(10.0, 1000.0, 3.0);

	for(int i=0; i<numAtoms; i++)
	{
		pg->addPoint(m_AtomList[i]);
	}


/* creating output file */

	FILE *coloredSurfaceFile;
	coloredSurfaceFile = fopen(coloredSurfaceFileName.c_str(), "wt");

	if(!coloredSurfaceFile)
	{
		cout<<"could not create output file "<< coloredSurfaceFileName << endl;
		return false;
	}

/* opening surface file */
	FILE *surfaceFile;

	surfaceFile = fopen(surfaceFileName.c_str(), "rt");
	if(!surfaceFile)
	{
		cout<<"could not open surface file "<< surfaceFileName.c_str() << endl;
		return false;
	}

	int numVert, numTriang;
	if(fscanf(surfaceFile, "%d", &numVert) != 1)
	{
		cout<<"could not read number of vertices "<< surfaceFileName.c_str() << endl;
		return false;
	}
	if(fscanf(surfaceFile, "%d", &numTriang) != 1)
	{
		cout<<"could not read number of triangles "<< surfaceFileName.c_str() << endl;
		return false;
	}

	fprintf(coloredSurfaceFile, "%d %d\n", numVert, numTriang);

/* Reading points from the surface file and computing the charge density values and saving them */

	double min = 100000;
	double max  = -100000;

	for(int i=0; i<numVert;i++)
	{
		double x, y, z, nx, ny, nz;

		if(fscanf(surfaceFile, "%lf %lf %lf %lf %lf %lf", &x, &y, &z, &nx, &ny, &nz ) != 6)
		{
			cout<<"could not read vertex "<< i+1 << endl;
			return false;
		}

		DPG::Point *p = new DPG::Point(x, y, z);
		vector <DPG::Point*> results;
		results = pg->range(p, distanceCutoff);

		double charge = 0.0;

		int numAtoms = results.size();
		if( numAtoms > 0)
		{
			for(int j=0; j<numAtoms; j++)
			{
				PDBParser::Atom* atom = (PDBParser::Atom*)results[j];

				double q = atom->getCharge();

				if(mode == 1) 
					charge += q;

				else if(mode == 0)
				{
					double dx = x - atom->getX();
					double dy = y - atom->getY();
					double dz = z - atom->getZ();

					double distsq = dx*dx + dy*dy + dz*dz;

					double r = atom->getRadius();

					double expval = -blobbyness + (blobbyness*distsq)/(r*r);

					charge += q*exp(expval);
				}
			}

			if(mode == 1)
			{
				charge /= (double)numAtoms;
			}
		}

		results.clear();

		double cr = 1.0, cg = 1.0, cb = 1.0;

		if(charge > max) max = charge;
		if(charge < min) min = charge;

		printf("%lf\n", charge);

		if(charge < 0.0)
		{
			cr = 1.0;
			cg = 1.0 + charge*0.5;
			cb = 1.0 + charge*0.5;
		}
		else if(charge > 0.0)
		{
			cb = 1.0;
			cr = 1.0 - charge*0.5 ;
			cg = 1.0 - charge*0.5 ;
		}

		fprintf(coloredSurfaceFile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", x, y, z, nx, ny, nz, cr, cg, cb);

	}

	printf("%lf %lf\n", max, min);

/* Reading triangles and saving them */

	for(int i=0; i<numTriang;i++)
	{
		int t1, t2, t3;

		if(fscanf(surfaceFile, "%d %d %d", &t1, &t2, &t3 ) != 3)
		{
			cout<<"could not read triangle "<< i+1 << endl;
			return false;
		}

		fprintf(coloredSurfaceFile, "%d %d %d\n", t1, t2, t3);
	}

	m_AtomList.clear();
	fclose(surfaceFile);
	fclose(coloredSurfaceFile);
	return true;
}
