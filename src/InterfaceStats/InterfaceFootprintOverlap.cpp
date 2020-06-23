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
#include<InterfaceStats/InterfaceFootprintOverlap.h>


InterfaceFootprintOverlap::InterfaceFootprintOverlap(string sf1, string sf2, string sf3, double iw, bool v)
{
	verbose = v;

	receptorSurfaceFileName = sf1;
	ligandSurfaceFileName1 = sf2;
	ligandSurfaceFileName2 = sf3;

	interfaceWidth = iw;

	filesRead = false;
	staticInterfaceComputed = false;
	movingInterfaceComputed = false;

	initialize();
}

InterfaceFootprintOverlap::~InterfaceFootprintOverlap()
{
	cleanup();
}


void InterfaceFootprintOverlap::cleanup()
{
	receptorPoints.clear();
	staticLigandPoints.clear();
	movingLigandPoints.clear();
	movingLigandPointsXformed.clear();

	receptorTriangles.clear();
	staticLigandTriangles.clear();
	movingLigandTriangles.clear();

	if(receptorPointsOnInterfaceWithLigand1) delete receptorPointsOnInterfaceWithLigand1;
	if(receptorPointsOnInterfaceWithLigand2) delete receptorPointsOnInterfaceWithLigand2;
	if(ligandPointsOnInterface1) delete ligandPointsOnInterface1;
	if(ligandPointsOnInterface2) delete ligandPointsOnInterface2;
	if(receptorPointsOnFootPrintOverlap) delete receptorPointsOnFootPrintOverlap;

	if(pgR) delete pgR;
	if(pgL1) delete pgL1;
	if(pgL2) delete pgL2;

	if(receptorSurfaceFile) fclose(receptorSurfaceFile);
	if(ligandSurfaceFile1) fclose(ligandSurfaceFile1);
	if(ligandSurfaceFile2) fclose(ligandSurfaceFile2);

	if(receptorHighlightedSurface1) fclose(receptorHighlightedSurface1);
	if(receptorHighlightedSurface2) fclose(receptorHighlightedSurface2);
	if(receptorFootprintSurface) fclose(receptorFootprintSurface);

	if(ligandHighlightedSurface1) fclose(ligandHighlightedSurface1);
	if(ligandHighlightedSurface2) fclose(ligandHighlightedSurface2);
}



bool InterfaceFootprintOverlap::openFiles()
{
	filesRead = false;

	receptorSurfaceFile = fopen(receptorSurfaceFileName.c_str(), "rt");
	if(!receptorSurfaceFile)
	{
		cout<<"could not open surface file "<< receptorSurfaceFileName.c_str() << endl;
		return false;
	}

	ligandSurfaceFile1 = fopen(ligandSurfaceFileName1.c_str(), "rt");
	if(!ligandSurfaceFile1)
	{
		cout<<"could not open surface file "<< ligandSurfaceFileName1.c_str() << endl;
		fclose(receptorSurfaceFile);
		return false;
	}

	ligandSurfaceFile2 = fopen(ligandSurfaceFileName2.c_str(), "rt");
	if(!ligandSurfaceFile2)
	{
		cout<<"could not open surface file "<< ligandSurfaceFileName2.c_str() << endl;
		fclose(receptorSurfaceFile);
		fclose(ligandSurfaceFile1);
		return false;
	}

	return true;
}


bool InterfaceFootprintOverlap::readFiles()
{
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

	receptorPointsOnInterfaceWithLigand1 = new bool[numVertR];
	receptorPointsOnInterfaceWithLigand2 = new bool[numVertR];
	receptorPointsOnFootPrintOverlap = new bool[numVertR];

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

		dx = receptorPoints[index2]->getX() - receptorPoints[index1]->getX();
		dy = receptorPoints[index2]->getY() - receptorPoints[index1]->getY();
		dz = receptorPoints[index2]->getZ() - receptorPoints[index1]->getZ();
		a = sqrt(dx*dx + dy*dy + dz*dz);

		dx = receptorPoints[index3]->getX() - receptorPoints[index2]->getX();
		dy = receptorPoints[index3]->getY() - receptorPoints[index2]->getY();
		dz = receptorPoints[index3]->getZ() - receptorPoints[index2]->getZ();
		b = sqrt(dx*dx + dy*dy + dz*dz);

		dx = receptorPoints[index1]->getX() - receptorPoints[index3]->getX();
		dy = receptorPoints[index1]->getY() - receptorPoints[index3]->getY();
		dz = receptorPoints[index1]->getZ() - receptorPoints[index3]->getZ();
		c = sqrt(dx*dx + dy*dy + dz*dz);

		s = (a+b+c)/2.0;

		area = sqrt(s*(s-a)*(s-b)*(s-c));

		triangle tr;
		tr.vertex[0] = index1;
		tr.vertex[1] = index2;
		tr.vertex[2] = index3;
		tr.area = area;

		receptorTriangles.push_back(tr);		

		receptorPoints[index1]->setWeight( receptorPoints[index1]->getWeight() + area/3.0);
		receptorPoints[index2]->setWeight( receptorPoints[index2]->getWeight() + area/3.0);
		receptorPoints[index3]->setWeight( receptorPoints[index3]->getWeight() + area/3.0);
	}

	fclose(receptorSurfaceFile);

	if(verbose) cout<<"Triangles read from receptor surface\n";



	/* Reading points from the ligand 1 surface file and inserting into DPG */

	if(fscanf(ligandSurfaceFile1, "%d %d", &numVertL1, &numTriangL1) != 2)
	{
		cout<<"Wrong format in ligand surface file\n";
		return false;
	}

	if(verbose) cout<<"Ligand 1 has: "<< numVertL1 << " vertices and " << numTriangL1 << " triangles\n";

	ligandPointsOnInterface1 = new bool[numVertL1];

	for(int i=0; i<numVertL1; i++)
	{
		if(fscanf(ligandSurfaceFile1, "%lf %lf %lf %lf %lf %lf", &x, &y, &z, &nx, &ny, &nz) != 6)
		{
			cout<<"Wrong format in ligand 1 surface file. At vertex no: "<< i <<"\n";
			cout<< "x = " << x << ", y = " << y << ", z = " << z << ", nx = " << nx << ", ny = " << ny << ", nz = " << nz << "\n";
			return false;
		}

		DPG::WeightedPoint *p = new DPG::WeightedPoint(x, y, z, 0.0);
		staticLigandPoints.push_back(p);
	}

	if(verbose) cout<<"Points read from ligand 1 surface and DPG created\n";

	for(int i=0; i<numTriangL1; i++)				
	{
		if(fscanf(ligandSurfaceFile1, "%d %d %d", &index1, &index2, &index3) != 3)
		{
			cout<<"Wrong format in ligand surface file\n";
			return false;
		}

		dx = staticLigandPoints[index2]->getX() - staticLigandPoints[index1]->getX();
		dy = staticLigandPoints[index2]->getY() - staticLigandPoints[index1]->getY();
		dz = staticLigandPoints[index2]->getZ() - staticLigandPoints[index1]->getZ();
		a = sqrt(dx*dx + dy*dy + dz*dz);

		dx = staticLigandPoints[index3]->getX() - staticLigandPoints[index2]->getX();
		dy = staticLigandPoints[index3]->getY() - staticLigandPoints[index2]->getY();
		dz = staticLigandPoints[index3]->getZ() - staticLigandPoints[index2]->getZ();
		b = sqrt(dx*dx + dy*dy + dz*dz);

		dx = staticLigandPoints[index1]->getX() - staticLigandPoints[index3]->getX();
		dy = staticLigandPoints[index1]->getY() - staticLigandPoints[index3]->getY();
		dz = staticLigandPoints[index1]->getZ() - staticLigandPoints[index3]->getZ();
		c = sqrt(dx*dx + dy*dy + dz*dz);

		s = (a+b+c)/2.0;

		area = sqrt(s*(s-a)*(s-b)*(s-c));

		triangle tr;
		tr.vertex[0] = index1;
		tr.vertex[1] = index2;
		tr.vertex[2] = index3;
		tr.area = area;

		staticLigandTriangles.push_back(tr);	

		staticLigandPoints[index1]->setWeight( staticLigandPoints[index1]->getWeight() + area/3.0);
		staticLigandPoints[index2]->setWeight( staticLigandPoints[index2]->getWeight() + area/3.0);
		staticLigandPoints[index3]->setWeight( staticLigandPoints[index3]->getWeight() + area/3.0);
	
	}

	fclose(ligandSurfaceFile1);

	if(verbose) cout<<"Triangles read from ligand 1 surface\n";




	/* Reading points from the second surface file and inserting into DPG */

	if(fscanf(ligandSurfaceFile2, "%d %d", &numVertL2, &numTriangL2) != 2)
	{
		cout<<"Wrong format in ligand 2 surface file\n";
		return false;
	}

	if(verbose) cout<<"Ligand 2 has: "<< numVertL2 << " vertices and " << numTriangL2 << " triangles\n";

	ligandPointsOnInterface2 = new bool[numVertL2];

	for(int i=0; i<numVertL2; i++)
	{
		if(fscanf(ligandSurfaceFile2, "%lf %lf %lf %lf %lf %lf", &x, &y, &z, &nx, &ny, &nz) != 6)
		{
			cout<<"Wrong format in ligand 2 surface file. At vertex no: "<< i <<"\n";
			cout<< "x = " << x << ", y = " << y << ", z = " << z << ", nx = " << nx << ", ny = " << ny << ", nz = " << nz << "\n";
			return false;
		}

		DPG::WeightedPoint *p = new DPG::WeightedPoint(x, y, z, 0.0);
		movingLigandPoints.push_back(p);
	}

	if(verbose) cout<<"Points read from ligand 2 surface and DPG created\n";

	for(int i=0; i<numTriangL2; i++)				
	{
		if(fscanf(ligandSurfaceFile2, "%d %d %d", &index1, &index2, &index3) != 3)
		{
			cout<<"Wrong format in ligand surface file 2\n";
			return false;
		}

		dx = movingLigandPoints[index2]->getX() - movingLigandPoints[index1]->getX();
		dy = movingLigandPoints[index2]->getY() - movingLigandPoints[index1]->getY();
		dz = movingLigandPoints[index2]->getZ() - movingLigandPoints[index1]->getZ();
		a = sqrt(dx*dx + dy*dy + dz*dz);

		dx = movingLigandPoints[index3]->getX() - movingLigandPoints[index2]->getX();
		dy = movingLigandPoints[index3]->getY() - movingLigandPoints[index2]->getY();
		dz = movingLigandPoints[index3]->getZ() - movingLigandPoints[index2]->getZ();
		b = sqrt(dx*dx + dy*dy + dz*dz);

		dx = movingLigandPoints[index1]->getX() - movingLigandPoints[index3]->getX();
		dy = movingLigandPoints[index1]->getY() - movingLigandPoints[index3]->getY();
		dz = movingLigandPoints[index1]->getZ() - movingLigandPoints[index3]->getZ();
		c = sqrt(dx*dx + dy*dy + dz*dz);

		s = (a+b+c)/2.0;

		area = sqrt(s*(s-a)*(s-b)*(s-c));

		triangle tr;
		tr.vertex[0] = index1;
		tr.vertex[1] = index2;
		tr.vertex[2] = index3;
		tr.area = area;

		movingLigandTriangles.push_back(tr);	

		movingLigandPoints[index1]->setWeight( movingLigandPoints[index1]->getWeight() + area/3.0);
		movingLigandPoints[index2]->setWeight( movingLigandPoints[index2]->getWeight() + area/3.0);
		movingLigandPoints[index3]->setWeight( movingLigandPoints[index3]->getWeight() + area/3.0);	
	}

	fclose(ligandSurfaceFile2);

	if(verbose) cout<<"Triangles read from ligand 2 surface\n";

	filesRead = true;

	return true;
}


bool InterfaceFootprintOverlap::initializeStaticData()
{
	pgR = new DPG::PG(10.0, 1000.0, 3.0);
	for(int i=0; i<numVertR; i++)
	{
		receptorPointsOnInterfaceWithLigand1[i] = false;
		pgR->addPoint(receptorPoints[i]);
	}

	if(verbose) cout<<"DPG for receptor created\n";

	pgL1 = new DPG::PG(10.0, 1000.0, 3.0);
	for(int i=0; i<numVertL1; i++)
	{
		ligandPointsOnInterface1[i] = false;
		pgL1->addPoint(staticLigandPoints[i]);
	}

	if(verbose) cout<<"DPG for static ligand created\n";

	return true;
}


bool InterfaceFootprintOverlap::detectInterfacePointsForStaticLigand()
{
	areaOfReceptorInterfaceWithLigand1 = 0.0;
	areaOfLigandInterface1 = 0.0;	

	// finding points on receptor surface which are on the interface with ligand 1
	for(int i=0; i<numVertR; i++)
	{
		if(pgL1->pointsWithinRange(receptorPoints[i], interfaceWidth))
		{
			receptorPointsOnInterfaceWithLigand1[i] = true;
			areaOfReceptorInterfaceWithLigand1 += receptorPoints[i]->getWeight();
		}
		else
		{
			receptorPointsOnInterfaceWithLigand1[i] = false;
		}
	}

	if(verbose) cout<<"Receptor's interface with static ligand computed\n";


	// finding points on ligand 1 surface which are on the interface with receptor
	for(int i=0; i<numVertL1; i++)
	{
		if(pgR->pointsWithinRange(staticLigandPoints[i], interfaceWidth))
		{
			ligandPointsOnInterface1[i] = true;
			areaOfLigandInterface1 += staticLigandPoints[i]->getWeight();
		}
		else
		{
			ligandPointsOnInterface1[i] = false;
		}
	}

	if(verbose) cout<<"Static Ligand's interface with receptor computed\n";

	staticInterfaceComputed = true;

	return true;
}


bool InterfaceFootprintOverlap::initialize()
{
	if(!filesRead)
	{
		if(!openFiles())
		{
			cout<<"Some of the files were not opened correctly.\n";
			return false;
		}

		if(!readFiles())
		{
			cout<<"Some of the files were not formatted correctly.\n";
			return false;
		}
	}

	if(!staticInterfaceComputed)
	{
		if(!initializeStaticData())
		{
			cout<<"Failed to create DPG for receptor and ligand 1.\n";
			return false;
		}

		if(!detectInterfacePointsForStaticLigand())
		{
			cout<<"Could not compute interface with ligand 1. Terminating.\n";
			cleanup();
			return false;
		}
	}
	return (filesRead && staticInterfaceComputed);
}


bool InterfaceFootprintOverlap::initializeMovingData()
{
	pgL2 = new DPG::PG(10.0, 1000.0, 3.0);
	movingLigandPointsXformed.clear();

	double x, y, z;
	double newx, newy, newz;

	for(int i=0; i<numVertL2; i++)
	{
		ligandPointsOnInterface2[i] = false;

		x = movingLigandPoints[i]->getX();
		y = movingLigandPoints[i]->getY();
		z = movingLigandPoints[i]->getZ();

		newx = xform.get(0,0)*x + xform.get(0,1)*y + xform.get(0,2)*z + xform.get(0,3);
		newy = xform.get(1,0)*x + xform.get(1,1)*y + xform.get(1,2)*z + xform.get(1,3);
		newz = xform.get(2,0)*x + xform.get(2,1)*y + xform.get(2,2)*z + xform.get(2,3);

		DPG::WeightedPoint *p = new DPG::WeightedPoint(newx, newy, newz, movingLigandPoints[i]->getWeight());
		pgL2->addPoint(p);
		movingLigandPointsXformed.push_back(p);
	}

	if(verbose) cout<<"DPG for moving ligand created\n";

	for(int i=0; i<numVertR; i++)
	{
		receptorPointsOnInterfaceWithLigand2[i] = false;
		receptorPointsOnFootPrintOverlap[i] = false;
	}

	return true;
}


bool InterfaceFootprintOverlap::detectInterfacePointsForMovingLigand()
{
	areaOfReceptorInterfaceWithLigand2 = 0.0;
	areaOfLigandInterface2 = 0.0;	
	areaOfOverlap = 0.0;

	// finding points on receptor surface which are on the interface with ligand 2
	for(int i=0; i<numVertR; i++)
	{
		if(pgL2->pointsWithinRange(receptorPoints[i], interfaceWidth))
		{
			receptorPointsOnInterfaceWithLigand2[i] = true;
			areaOfReceptorInterfaceWithLigand2 += receptorPoints[i]->getWeight();

			if(receptorPointsOnInterfaceWithLigand1[i])
			{
				receptorPointsOnFootPrintOverlap[i] = true;
				areaOfOverlap += receptorPoints[i]->getWeight();
			}
			else
			{
				receptorPointsOnFootPrintOverlap[i] = false;
			}
		}
		else
		{
			receptorPointsOnInterfaceWithLigand2[i] = false;
		}
	}

	if(verbose) cout<<"Receptor's interface with moving ligand computed\n";


	// finding points on ligand 2 surface which are on the interface with receptor
	for(int i=0; i<numVertL2; i++)
	{
		if(pgR->pointsWithinRange(movingLigandPointsXformed[i], interfaceWidth))
		{
			ligandPointsOnInterface2[i] = true;
			areaOfLigandInterface2 += movingLigandPointsXformed[i]->getWeight();
		}
		else
		{
			ligandPointsOnInterface2[i] = false;
		}
	}

	if(verbose) cout<<"Moving Ligand's interface with receptor computed\n";


	movingInterfaceComputed = true;

	return true;
}



bool InterfaceFootprintOverlap::computeInterfaces(Matrix mtx)
{
	movingInterfaceComputed = false;

	xform = mtx;

	if(!initialize())
	{
		cout<<"Could not initialize. Terminating.\n";
		cleanup();
		return false;
	}

	if(!initializeMovingData())
	{
		cout<<"Failed to create DPG for moving ligand.\n";
		return false;
	}

	if(!detectInterfacePointsForMovingLigand())
	{
		cout<<"Could not compute interface with ligand 2. Terminating.\n";
		cleanup();
		return false;
	}

	return true;
}



bool InterfaceFootprintOverlap::prepareStaticInterfaceSurfaceFiles()
{
	if(!staticInterfaceComputed)
	{
		cout<<"Static interface is not ready.\n";
		return false;
	}

	receptorHighlightedSurfaceFileName1 = receptorSurfaceFileName.substr(0, receptorSurfaceFileName.length() - 5) + "_interface_with_" + ligandSurfaceFileName1.substr(0, ligandSurfaceFileName1.length() - 5) +  ".rawc";
	ligandHighlightedSurfaceFileName1 = ligandSurfaceFileName1.substr(0, ligandSurfaceFileName1.length() - 5) + "_interface_with_" + receptorSurfaceFileName.substr(0, receptorSurfaceFileName.length() - 5) +  ".rawc";

	receptorHighlightedSurface1 = fopen(receptorHighlightedSurfaceFileName1.c_str(), "wt");
	if(!receptorHighlightedSurface1)
	{
		cout<<"could not create surface file "<< receptorHighlightedSurfaceFileName1.c_str() << endl;
		return false;
	}


	ligandHighlightedSurface1 = fopen(ligandHighlightedSurfaceFileName1.c_str(), "wt");
	if(!ligandHighlightedSurface1)
	{
		cout<<"could not create surface file "<< ligandHighlightedSurfaceFileName1.c_str() << endl;
		fclose(receptorHighlightedSurface1);
		return false;
	}


	// preparing interface surface for receptor
	fprintf(receptorHighlightedSurface1, "%d %d\n", numVertR, numTriangR);

	for(int i=0; i<numVertR; i++)
	{
		fprintf(receptorHighlightedSurface1, "%lf %lf %lf", receptorPoints[i]->getX(), receptorPoints[i]->getY(), receptorPoints[i]->getZ());

		if(receptorPointsOnInterfaceWithLigand1[i])
			fprintf(receptorHighlightedSurface1, " %lf %lf %lf\n", 0.67, 0.33, 0.5);
		else
			fprintf(receptorHighlightedSurface1, " %lf %lf %lf\n", 0.33, 0.67, 1.0);
	}

	for(int i=0; i<numTriangR; i++)
	{
		fprintf(receptorHighlightedSurface1, "%d %d %d\n", receptorTriangles[i].vertex[0], receptorTriangles[i].vertex[1], receptorTriangles[i].vertex[2]);
	}


	// preparing interface surface for ligand 1
	fprintf(ligandHighlightedSurface1, "%d %d\n", numVertL1, numTriangL1);

	for(int i=0; i<numVertL1; i++)
	{
		fprintf(ligandHighlightedSurface1, "%lf %lf %lf", staticLigandPoints[i]->getX(), staticLigandPoints[i]->getY(), staticLigandPoints[i]->getZ());

		if(ligandPointsOnInterface1[i])
			fprintf(ligandHighlightedSurface1, " %lf %lf %lf\n", 0.33, 0.67, 0.5);
		else
			fprintf(ligandHighlightedSurface1, " %lf %lf %lf\n", 0.67, 0.33, 1.0);
	}

	for(int i=0; i<numTriangL1; i++)
	{
		fprintf(ligandHighlightedSurface1, "%d %d %d\n", staticLigandTriangles[i].vertex[0], staticLigandTriangles[i].vertex[1], staticLigandTriangles[i].vertex[2]);
	}

	fclose(receptorHighlightedSurface1);
	fclose(ligandHighlightedSurface1);

	return true;
}



bool InterfaceFootprintOverlap::prepareMovingInterfaceSurfaceFiles(int poseId)
{
	if(!movingInterfaceComputed)
	{
		cout<<"Static interface is not ready.\n";
		return false;
	}

	stringstream ss;
	ss << poseId;
	string pose = ss.str();

	receptorHighlightedSurfaceFileName2 = receptorSurfaceFileName.substr(0, receptorSurfaceFileName.length() - 5) + "_interface_with_" + ligandSurfaceFileName2.substr(0, ligandSurfaceFileName2.length() - 5) + "_pose_" + pose + ".rawc";

	ligandHighlightedSurfaceFileName2 = ligandSurfaceFileName2.substr(0, ligandSurfaceFileName2.length() - 5) + "_pose_" + pose +  "_interface_with_" + receptorSurfaceFileName.substr(0, receptorSurfaceFileName.length() - 5) +  ".rawc";

	receptorFootprintSurfaceFileName = receptorSurfaceFileName.substr(0, receptorSurfaceFileName.length() - 5) + "_footprint_overlap_of_" + ligandSurfaceFileName1.substr(0, ligandSurfaceFileName1.length() - 5) + "_and_" + ligandSurfaceFileName2.substr(0, ligandSurfaceFileName2.length() - 5) + "_pose_" + pose +  ".rawc";

	receptorHighlightedSurface2 = fopen(receptorHighlightedSurfaceFileName2.c_str(), "wt");
	if(!receptorHighlightedSurface2)
	{
		cout<<"could not create surface file "<< receptorHighlightedSurfaceFileName2.c_str() << endl;
		return false;
	}

	ligandHighlightedSurface2 = fopen(ligandHighlightedSurfaceFileName2.c_str(), "wt");
	if(!ligandHighlightedSurface2)
	{
		cout<<"could not create surface file "<< ligandHighlightedSurfaceFileName2.c_str() << endl;
		fclose(receptorHighlightedSurface2);
		return false;
	}

	receptorFootprintSurface = fopen(receptorFootprintSurfaceFileName.c_str(), "wt");
	if(!receptorFootprintSurface)
	{
		cout<<"could not create surface file "<< receptorFootprintSurfaceFileName.c_str() << endl;
		fclose(receptorHighlightedSurface2);
		fclose(ligandHighlightedSurface2);
		return false;
	}


	// preparing interface surface for receptor
	fprintf(receptorHighlightedSurface2, "%d %d\n", numVertR, numTriangR);

	for(int i=0; i<numVertR; i++)
	{
		fprintf(receptorHighlightedSurface2, "%lf %lf %lf", receptorPoints[i]->getX(), receptorPoints[i]->getY(), receptorPoints[i]->getZ());

		if(receptorPointsOnInterfaceWithLigand2[i])
			fprintf(receptorHighlightedSurface2, " %lf %lf %lf\n", 1.0, 0.67, 0.0);		
		else
			fprintf(receptorHighlightedSurface2, " %lf %lf %lf\n", 0.33, 0.67, 1.0);
	}

	for(int i=0; i<numTriangR; i++)
	{
		fprintf(receptorHighlightedSurface2, "%d %d %d\n", receptorTriangles[i].vertex[0], receptorTriangles[i].vertex[1], receptorTriangles[i].vertex[2]);
	}


	// preparing interface surface for ligand 2
	fprintf(ligandHighlightedSurface2, "%d %d\n", numVertL2, numTriangL2);

	for(int i=0; i<numVertL2; i++)
	{
		fprintf(ligandHighlightedSurface2, "%lf %lf %lf", movingLigandPoints[i]->getX(), movingLigandPoints[i]->getY(), movingLigandPoints[i]->getZ());

		if(ligandPointsOnInterface2[i])
			fprintf(ligandHighlightedSurface2, " %lf %lf %lf\n", 0.33, 0.67, 0.5);
		else
			fprintf(ligandHighlightedSurface2, " %lf %lf %lf\n", 1.0, 0.67, 0.5);
	}

	for(int i=0; i<numTriangL2; i++)
	{
		fprintf(ligandHighlightedSurface2, "%d %d %d\n", movingLigandTriangles[i].vertex[0], movingLigandTriangles[i].vertex[1], movingLigandTriangles[i].vertex[2]);
	}


	// preparing overlap surface for receptor
	fprintf(receptorFootprintSurface, "%d %d\n", numVertR, numTriangR);

	for(int i=0; i<numVertR; i++)
	{
		fprintf(receptorFootprintSurface, "%lf %lf %lf", receptorPoints[i]->getX(), receptorPoints[i]->getY(), receptorPoints[i]->getZ());

		if(receptorPointsOnFootPrintOverlap[i])
			fprintf(receptorFootprintSurface, " %lf %lf %lf\n", 0.0, 0.33, 0.0);		
		else if(receptorPointsOnInterfaceWithLigand1[i])
			fprintf(receptorFootprintSurface, " %lf %lf %lf\n", 0.67, 0.33, 0.5);
		else if(receptorPointsOnInterfaceWithLigand2[i])
			fprintf(receptorFootprintSurface, " %lf %lf %lf\n", 1.0, 0.67, 0.0);
		else
			fprintf(receptorFootprintSurface, " %lf %lf %lf\n", 0.33, 0.67, 1.0);
	}

	for(int i=0; i<numTriangR; i++)
	{
		fprintf(receptorFootprintSurface, "%d %d %d\n", receptorTriangles[i].vertex[0], receptorTriangles[i].vertex[1], receptorTriangles[i].vertex[2]);
	}

	fclose(receptorHighlightedSurface2);
	fclose(ligandHighlightedSurface2);	
	fclose(receptorFootprintSurface);
	
	return true;
}
