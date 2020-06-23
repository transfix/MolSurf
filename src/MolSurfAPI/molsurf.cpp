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
#include <MolSurfAPI/molsurf.h>
#include <PDBParser/FlattenGOA.h>

using namespace CVCUtility;
using namespace ASPLINE;

// Internal helper functions not exposed to the interface

bool isCloseTo(double x1, double y1, double z1, double r, double x2, double y2, double z2, double* d)
{
	double d1 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
	*d = sqrt(d1);
	if(d1 < r*r)
	{
		return true;
	}
	return false;
}

#ifdef WITH_CG
void ParseCGInfo(FILE* fp, MOLECULE::GroupOfAtoms* groupOfBeads)
{
	ssize_t read;
	size_t len = 0;
	char* token = NULL;
	char* line = NULL;
	char* delim = (char*)" \t\n";
	vector<int> atomList;
	while((read = getline(&line,&len,fp)) != -1)
	{
		token = strtok(line,delim);
		if(token == NULL)
		{
			break;
		}
		if(!strcmp(token,"Bead"))
		{
			token = strtok(NULL,delim);
			int beadIndex = atoi(token);
			if((read = getline(&line,&len,fp)) != -1)
			{
				token = strtok(line,delim);
				while(token != NULL)
				{
					atomList.push_back(atoi(token));
					token = strtok(NULL,delim);
				}
				groupOfBeads->m_Atoms[beadIndex]->child = new int[atomList.size()];
				for(int i = 0; i < atomList.size(); i++)
				{
					groupOfBeads->m_Atoms[beadIndex]->child[i] = atomList[i];
				}
				groupOfBeads->m_Atoms[beadIndex]->numOfChildren = atomList.size();
				atomList.clear();
			}
		}
	}
}

void loadBornRadius(MOLECULE::GroupOfAtoms* molecule, FILE* fp)
{
	double br;
	for(int i = 0; i < molecule->numOfAtoms; i++)
	{
		if(fscanf(fp, "%lf\n", &br) != 1) error("Couldn't read Born radius!");
		molecule->m_Atoms[i]->Born_radius = br;
	}
}
#endif // WITH_CG

//Use this to read in a molecule from file
PDBParser::GroupOfAtoms* loadGOA(const string& file)
{
	PDBParser::GroupOfAtoms* molecule = GOALoader().loadFile(file);
	if(!molecule)
	{
		exit(-1);
	}
	return molecule;
}

//Write a surface to file
void saveAndDeleteGeometry(const string& file, Geometry* geometry)
{

	if(!GeometryLoader().saveFile(file.c_str(), geometry))
	{
		delete geometry;
		exit(-1); // FIXME: don't exit... rather handle errors gracefully
	}
	

	// arand: no longer Deleting geometry... things get messed up...
	//delete geometry;
}

// Implementation of the interface


// Get Area
void usageArea()
{
	cout	<< endl << "Usage: MolSurf -area <raw input file> <txt output file>"
		<< endl << "    Compute the surface area of a triangular mesh."
		<< endl ;
}

bool area(int argc, char* argv[])
{
	// First read in the arguments
	if(argc != 3)	// If the wrong number of arguments, print usage
	{
		usageArea();
		return false;
	}
	string input  = string(argv[2]);
//	FILE* output = fileWrite(argv[3]);

	Geometry* geometry = GeometryLoader().loadFile(input.c_str());
	double myArea = geometry->getArea();
	delete geometry;

	printf(" %lf ", myArea);
//	fprintf(output, "%lf\n", myArea);
//	fclose(output);
	return true;
}


double area2(char * rawnFile)
{
  Geometry* geometry = GeometryLoader().loadFile(rawnFile);

  if (geometry) {
    double myArea = geometry->getArea();
    delete geometry;
    return myArea;
  } 
  return -1.0;
}


Geometry* aSplineRaw(Geometry* geometry, int resolution)
{
	ASpline aspline = ASpline(geometry);
	aspline.Triangle_Bezier_Patch(resolution);
	return aspline.getGeometry();
}

Geometry* aSplineRawQuad(Geometry* geometry, int numOfPts, int resolution, string output1, string output2, string type)
{
	ASpline aspline = ASpline(geometry);
	aspline.Triangle_Bezier_Patch_With_Quadrature(type.c_str(), numOfPts, resolution);
	aspline.printRawnWithQuadrature(output1.c_str(),resolution,numOfPts);
	aspline.printQuadratureResults(output2.c_str());
	return aspline.getGeometry();
}

// aspline quadrature
void usageASplineQuad()
{
	cout	<< endl << "Usage: MolSurf -aSpline -quad <raw input file> <quadrature rule> "
		<< endl << "          <points per triangle> <txt output file>"
		<< endl << "    Generate quadrature points for an ASMS from a triangular mesh."
		<< endl ;
}

bool aSplineQuad(int argc, char* argv[])
{
	if(argc != 7)
	{
		usageASplineQuad();
		return false;
	}
	FILE* input	= fileRead(argv[3]);
	string type	= string(argv[4]);
	int numOfPts	= atoi(argv[5]);
	FILE* output	= fileWrite(argv[6]);
	ASpline* aspline = new ASpline(input);
	aspline->Triangle_Bezier_Patch_Nodes(type.c_str(), numOfPts);

	for(int i = 0; i < aspline->numbtris; i++)
	{
		aspline->faces[i]->print(output);
	}
	fclose(input);
	fclose(output);
	return true;
}

// aspline volume
void usageASplineVolume()
{
	cout	<< endl << "Usage: MolSurf -aspline -volume <raw input file> <txt output file>"
		<< endl << "    Compute the volume bounded by an ASMS from a triangular mesh."
		<< endl ;
}

bool aSplineVolume(int argc, char* argv[])
{
	// First read in the arguments
	if(argc != 5)	// If the wrong number of arguments, print usage
	{
		usageASplineVolume();
		return false;
	}
	FILE* input  = fileRead(argv[3]);
	FILE* output = fileWrite(argv[4]);

	ASpline* aspline =  new ASpline(input);
	fclose(input);
	aspline->Triangle_Patch_Volume();
	int myVolume = aspline->volume;

	fprintf(output, "%d\n", myVolume);
	fclose(output);
	return true;
}

#ifdef WITH_CG
// Born radius
void usageBornRadius()
{
	cout	<< endl << "Usage: MolSurf -bornRadius <radius type: vdw or pqr> <input pqr file> <surface geometry file> <output file>"
		<< endl << "    Undocumented get born radius. Good luck."
		<< endl ;
}

bool bornRadius(int argc, char* argv[])
{
	if(argc != 6)	// If the wrong number of arguments, print usage
	{
		usageBornRadius();
		return false;
	}
	string radiusString	= string(argv[2]);
	FILE* inputpqr		= fileRead(argv[3]);
	FILE* inputgeo		= fileRead(argv[4]);
	FILE* output		= fileWrite(argv[5]);

	// Get radius type
	MOLECULE::RADIUS_TYPE radiusType;
	if(radiusString.compare("vdw") == 0 || radiusString.compare("VDW") == 0)
	{
		radiusType = MOLECULE::VDW_RADIUS;
	}
	else if(radiusString.compare("pqr") == 0 || radiusString.compare("PQR") == 0)
	{
		radiusType = MOLECULE::PQR_RADIUS;
	}
	else
	{
		usageBornRadius();
		return false;
	}

	MOLECULE::GroupOfAtoms* molecule = ParsePQR(inputpqr, radiusType);
	//cout << molecule->numOfAtoms << " atoms" << endl;
	fclose(inputpqr);

	// Surface geometry file
	GEOMETRY::Surface* surface = GEOMETRY::GeometryParser().SurfaceTriangulation(inputgeo);
	fclose(inputgeo);

	// Output file
	surface->Triangle_Bezier_Patch();
	molecule->getBornRadii(surface);
	for(int i = 0; i < molecule->numOfAtoms; i++)
	{
		fprintf(output, "%f\n", molecule->m_Atoms[i]->Born_radius);
	}
	fclose(output);
	return true;
}

// Coarse grain charge
void usageCoarseGrainCharge()
{
	cout	<< endl << "Usage: MolSurf -coarseGrainCharge <radius type: vdw or pqr> <input file atomic model 1> < input file atomic model 2>"
		<< endl << "    <input file cg model 1> <input file cg model 2> <input file cg model 3> <output file>"
		<< endl << "    Undocumented coarse grain charge. Good luck."
		<< endl ;
}

bool coarseGrainCharge(int argc, char* argv[])
{
	if(argc != 6)
	{
		usageCoarseGrainCharge();
		return false;
	}
	string radiusString = string(argv[2]);
	string atomicModel  = string(argv[3]);
	FILE* atomicModel1  = fileRead(argv[3]);
	FILE* atomicModel2  = fileRead(argv[4]);
	FILE* cgModel1      = fileRead(argv[5]);
	FILE* cgModel2      = fileRead(argv[6]);
	FILE* cgModel3      = fileRead(argv[7]);
	FILE* output        = fileWrite(argv[8]);

	// Get radius type
	MOLECULE::RADIUS_TYPE radiusType;
	if(radiusString.compare("vdw") == 0 || radiusString.compare("VDW") == 0)
	{
		radiusType = MOLECULE::VDW_RADIUS;
	}
	else if(radiusString.compare("pqr") == 0 || radiusString.compare("PQR") == 0)
	{
		radiusType = MOLECULE::PQR_RADIUS;
	}
	else
	{
		usageCoarseGrainCharge();
		return false;
	}
	// Input file atomic model
	MOLECULE::GroupOfAtoms* molecule;
	if(strstr(atomicModel.c_str(), "optimal"))
	{
		molecule = ParseOptimalPQR(atomicModel1, radiusType);
	}
	else
	{
		molecule = ParsePQR(atomicModel1, radiusType);
	}
	fclose(atomicModel1);

	// Input file two
	loadBornRadius(molecule, atomicModel2);
	fclose(atomicModel2);

	// Input file cg model
	MOLECULE::GroupOfAtoms* groupOfBeads = ParsePQR(cgModel1, MOLECULE::PQR_RADIUS);
	fclose(cgModel1);

	// Input file cg model 2
	loadBornRadius(groupOfBeads, cgModel2);
	fclose(cgModel2);

	// Input file cg model 3
	ParseCGInfo(cgModel3, groupOfBeads);
	cout << molecule->numOfAtoms << " atoms, " << molecule->numOfRes << " residues, " << groupOfBeads->numOfAtoms << " cg beads" << endl;
	MOLECULE::CoarseGrain(molecule).CoarseGrainChargeGB2(groupOfBeads);

	// Output file
	savePQR(groupOfBeads, output);
	fclose(output);
	return true;
}

//Coarse grain
void usageCoarseGrain()
{
	cout	<< endl << "Usage: MolSurf -coarseGrain <radius type: pqr or vdw> <PDB/PQR input> <pqr output> <cg info output>"
	  //		<< endl << "    Improve the mesh quality using geometric flow."
		<< endl << "    If the input file name contains optimal, then we parse an optimal pqr."
		<< endl ;
}

bool coarseGrain(int argc, char* argv[])
{
	if(argc != 6)
	{
		usageCoarseGrain();
		return false;
	}
	// Get radius type
	string radiusString = string(argv[2]);
	string inputName    = string(argv[3]);
	FILE*  input        = fileRead(argv[3]);
	FILE*  pqrOutput    = fileWrite(argv[4]);
	FILE*  infoOutput   = fileWrite(argv[5]);

	MOLECULE::RADIUS_TYPE radiusType;
	if(radiusString.compare("vdw") == 0 || radiusString.compare("VDW") == 0)
	{
		radiusType = MOLECULE::VDW_RADIUS;
	}
	else if(radiusString.compare("pqr") == 0 || radiusString.compare("PQR") == 0)
	{
		radiusType = MOLECULE::PQR_RADIUS;
	}
	else
	{
		usageCoarseGrain();
		return false;
	}

	// Input file
	MOLECULE::GroupOfAtoms* molecule;
	if(strstr(inputName.c_str(), "optimal"))
	{
		molecule = ParseOptimalPQR(input, radiusType);
	}
	else
	{
		molecule = ParsePQR(input, radiusType);
	}
	fclose(input);
	cout << molecule->numOfAtoms << " atoms, " << molecule->numOfRes << "residues" << endl;

	MOLECULE::GroupOfAtoms* groupOfBeads = MOLECULE::CoarseGrain().CoarseGrainMolecule(molecule);
	cout << groupOfBeads->numOfAtoms << " beads" << endl;

	// Output file
	savePQR(groupOfBeads, pqrOutput);
	fclose(pqrOutput);

	// Output file 2
	for(int i = 0; i < groupOfBeads->numOfAtoms; i++)
	{
		fprintf(infoOutput, "Bead %d\n", i);
		for(int j = 0; j < groupOfBeads->m_Atoms[i]->numOfChildren; j++)
		{
			fprintf(infoOutput, "%d ", groupOfBeads->m_Atoms[i]->child[j]);
		}
		fprintf(infoOutput, "\n");
	}
	fclose(infoOutput);
	return true;
}


// Optimize CG Charge
void usageOptimizeCGCharge()
{
	cout	<< endl << "Usage: MolSurf -optimizeCGCharge <size> <iteration> <ARM> <Hierachical Information> < ????? >"
		<< endl << "    Undocumented optimize CG charge. Good luck."
		<< endl ;
}

bool optimizeCGCharge(int argc, char* argv[])
{
	if(argc != 7)
	{
		usageOptimizeCGCharge();
		return false;
	}
	int size        = atoi(argv[2]);
	int iteration   = atoi(argv[3]);
	char* pqrfiles  = argv[4];
	char* output    = argv[5];

	std::vector<PDBParser::Atom*> cg;

	ARM arm = ARM(&pqrfiles);
	arm.setHierarchicalInformation(output);
	arm.getOptimizedCG(cg,size,iteration,output);
	return true;
}

// Optimize CG Geometry
void usageOptimizeCGGeometry()
{
	cout	<< endl << "Usage: MolSurf -optimizeCGGeometry <size> <iterations> <pdb/pqr input file> "
		<< endl << "          [[[[UNUSED<coarse-grained pdb/pqr input file>]]]] <txt hierarchy info input file> "
		<< endl << "          <pbq/pqr output file>"
		<< endl << "    Optimize coarse-grained model (radius and centers) based on hierarchical "
		<< endl << "    information."
		<< endl ;
}

bool optimizeCGGeometry(int argc, char* argv[])
{
	if(argc != 8)
	{
		usageOptimizeCGGeometry();
		return false;
	}
	int size       = atoi(argv[2]);
	int iteration  = atoi(argv[3]);
	char* pqrfiles = argv[4];
	// argv[5] is unused!
	char* output   = argv[6];

	std::vector<PDBParser::Atom*> cg;

	ARM arm = ARM(&pqrfiles);
	arm.setHierarchicalInformation(output);
	arm.getOptimizedCG(cg, size, iteration, output);
	return true;
}


#endif

#ifndef WITH_CG

bool bornRadius(int argc, char* argv[]) {
  cout << "CG Disabled" << endl;
  return false;
}

bool coarseGrainCharge(int argc, char* argv[]) {
  cout << "CG Disabled" << endl;
  return false;
}

bool coarseGrain(int argc, char* argv[]) {
  cout << "CG Disabled" << endl;
  return false;
}

bool optimizeCGCharge(int argc, char* argv[]) {
  cout << "CG Disabled" << endl;
  return false;
}

bool optimizeCGGeometry(int argc, char* argv[]) {
  cout << "CG Disabled" << endl;
  return false;
}

#endif

// Convert to rawn using HLS
// arand: not sure if this works?
void usageConvertToRawnUsingHLS()
{
	cout	<< endl << "Usage: MolSurf -convertToRawnUsingHLS <pqr/pdb input file> <rawn output file>  <outputfile2> <size> [nIterations [weight]] "
		<< endl << "    Generate molecular surface with surface normal vectors using the higher "
		<< endl << "    order level-set algorithm.    "
		<< endl ;
}

bool convertToRawnUsingHLS(int argc, char* argv[])
{
	if(argc !=  6 && argc != 7 && argc != 8)
	{
		usageConvertToRawnUsingHLS();
		return false;
	}
	// Get dimensions
	string input   = string(argv[2]);
	string output1 = string(argv[3]);
	string output2 = string(argv[4]);
	int size       = atoi(argv[5]);
	// Get iterations
	int nIterations = 3;
	if(argc == 7 || argc == 8)
	{
		nIterations = atoi(argv[6]);
	}
	// Get weight
	double weight = 0.0001;
	if(argc == 8)
	{
		weight = atof(argv[7]);
	}

	// Computation proper
	
	PDBParser::GroupOfAtoms* molecule = loadGOA(input);

	unsigned int dim[3] = {size, size, size};
	SimpleVolumeData* sData = HLevelSet().getHigherOrderLevelSetSurface(molecule, dim, nIterations, weight);
	Geometry* geometry = GeometryLoader().loadFile(output1.c_str());
	saveAndDeleteGeometry(output2, geometry);

	return true;
}


// Decimate percentage
void usageDecimatePercentage()
{
	cout	<< endl << "Usage: MolSurf --decimatePercentage <raw input file> <raw output file>"
		<< endl << "                       <percentage> <max variation> <min angle>"
		<< endl << "     Decimate a mesh by vertex removal until fixed percentage remains."
		<< endl ;
}

bool decimatePercentage(int argc, char* argv[])
{
	if(argc != 7)
	{
		usageDecimatePercentage();
		return false;
	}
	string input     = string(argv[2]);
	string output    = string(argv[3]);
	float percentage = atof(argv[4]);
	float variation  = atof(argv[5]);
	float innerangle = atof(argv[6]);

	Decimation().DecimatePercentage(input.c_str(), output.c_str(), percentage, -1.0, variation, innerangle);
	return true;
}

// Decimate percentage
// arand - updated for arandDecimate.  Requires one less argument.
void usageDecimate()
{
	cout	<< endl << "Usage: MolSurf --decimate <raw input file> <raw output file> <iterations> <max variation> <min angle>"
		<< endl << "     Decimate a mesh by vertex removal."
		<< endl ;
}

bool decimate(int argc, char* argv[])
{
	if(argc != 7)
	{
		usageDecimate();
		return false;
	}
	string input     = string(argv[2]);
	string output    = string(argv[3]);
	int niter        = atoi(argv[4]);
	float variation  = atof(argv[5]);
	float innerangle = atof(argv[6]);

	// the NULL argument is no longer used...
	Decimation().Decimate(input.c_str(), output.c_str(), niter, NULL, variation, innerangle);
	return true;
}

// Compute Normals
void usageNormals()
{
	cout	<< endl << "Usage: MolSurf -normals -[average|loop] <raw input file> <rawn output file>"
		<< endl << "Or Usage: MolSurf -normals -invert <rawn input file> <rawn output file>"
		<< endl << "    Compute surface normals at vertices of a triangular mesh using face normal "
		<< endl << "    averaging or Loop's algorithm or invert the normal direction."
		<< endl ;
}

bool normals(int argc, char* argv[])
{
	if(argc != 5)
	{
		usageNormals();
		return false;
	}

	if(strcmp(argv[2], "-invert") == 0)
	{
		Geometry* geometry = GeometryLoader().loadFile(argv[3]);
		for(int i=0; i<geometry->m_NumTriVerts*3; i++)
		{
			geometry->m_TriVertNormals[i] *= -1.0; 
		}
		saveAndDeleteGeometry(string(argv[4]), geometry);
		return true;
	}
	else {
	bool method  = (strcmp(argv[2], "-loop") == 0);
	FILE* input  = fileRead(argv[3]);
	FILE* output = fileWrite(argv[4]);

	bool printColors = strstr(argv[4], ".rawnc");

	// Function proper
	COMPUTENORMALS::Surface* mySurface = COMPUTENORMALS::GeometryParser().SurfaceTriangulation(input);
	fclose(input);

	mySurface->Produce_Normal(method);
	fprintf(output, "%d %d\n", mySurface->numbpts, mySurface->numbtris);
	for(int i = 0; i < mySurface->numbpts; i++)
	{
	  if (printColors) {
		fprintf(output, "%f %f %f %f %f %f %f %f %f\n",
				mySurface->m_Vertices[i]->Point[0], mySurface->m_Vertices[i]->Point[1], mySurface->m_Vertices[i]->Point[2],
			mySurface->m_Vertices[i]->Normal[0],mySurface->m_Vertices[i]->Normal[1], mySurface->m_Vertices[i]->Normal[2],
			mySurface->m_Vertices[i]->Color[0],mySurface->m_Vertices[i]->Color[1], mySurface->m_Vertices[i]->Color[2]);
	  } else {
		fprintf(output, "%f %f %f %f %f %f\n",
				mySurface->m_Vertices[i]->Point[0], mySurface->m_Vertices[i]->Point[1], mySurface->m_Vertices[i]->Point[2],
				mySurface->m_Vertices[i]->Normal[0],mySurface->m_Vertices[i]->Normal[1], mySurface->m_Vertices[i]->Normal[2]);
	  }

	}
	for(int i = 0; i < mySurface->numbtris; i++)
	{
		fprintf(output, "%d %d %d\n",
				mySurface->m_Faces[i]->IndexInRAW[0], mySurface->m_Faces[i]->IndexInRAW[1], mySurface->m_Faces[i]->IndexInRAW[2]);
	}
	fclose(output);
  }
}

// Reorientate the surface
void usageReorientate()
{
	cout	<< endl << "Usage: MolSurf -reoriente <input geom file> <output geom file>"
		<< endl << "    change the orientation of the surface"
		<< endl ;
}

bool reOrientate(int argc, char* argv[])
{
	if(argc != 4)
	{
		usageReorientate();
		return false;
	}
	int a;
	Geometry* geometry = GeometryLoader().loadFile(argv[2]);
		for(int i=0; i<geometry->m_NumTris; i++)
		{	
			a =	geometry->m_Tris[3*i];
			geometry->m_Tris[3*i] = geometry->m_Tris[3*i+1];
			geometry->m_Tris[3*i+1] = a;

		}
		saveAndDeleteGeometry(string(argv[3]), geometry);
		return true;

}



#ifdef WITH_CGAL
// Call Pocket Tunnel
void usagePocketTunnel()
{
	cout	<< endl << "Usage: MolSurf -pocketTunnel <input file> <opc> <otc> <output file> <output file 2>"
		<< endl << "    Undocumented call pocket tunnel. Good luck."
		<< endl ;
}

bool pocketTunnel(int argc, char* argv[])
{
	if(argc != 7)
	{
		usagePocketTunnel();
		return false;
	}
	string inputFileName	= string(argv[2]);
	int opc 		= atoi(argv[3]);
	int otc			= atoi(argv[4]);
	string outputFileName1	= string(argv[5]);
	string outputFileName2	= string(argv[6]);

	Geometry* geometryout1 = new Geometry();
	Geometry* geometryout2 = new Geometry();
	Geometry* geometryin = GeometryLoader().loadFile(inputFileName);
	pocket_tunnel_fromsurf(geometryin, &geometryout1, &geometryout2, opc, otc);
	saveAndDeleteGeometry(outputFileName1, geometryout1);
	saveAndDeleteGeometry(outputFileName2, geometryout2);
	return true;
}
#endif //WITH_CGAL


#ifndef WITH_CGAL

bool pocketTunnel(int argc, char* argv[]) {
  cout << "PocketTunnel disable to remove CGAL dependence." << endl;
}


#endif

// Call pocket
void usagePocket()
{
	cout	<< endl << "Usage: MolSurf -pocket <pdbFileName> <distance> <size> <outputFileName>"
		<< endl << "    Undocumented call pocket. Good luck."
		<< endl ;
}

bool pocket(int argc, char* argv[])
{
	if(argc != 6)
	{
		usagePocket();
		return false;
	}
	string pdbFile	= string(argv[2]);
	float distance	= atof(argv[3]);
	int size	= atoi(argv[4]);
	string output	= string(argv[5]);
	PDBParser::GroupOfAtoms* molecule = GOALoader().loadFile(pdbFile);
	Geometry* geometry;
	Pocket().constructPockets(molecule, output.c_str(), distance, &geometry, size);
	return true;
}

// Populate SAS
void usagePopulateSAS()
{
	cout	<< endl << "Usage: MolSurf -populateSAS <pdb/pqr input file> <output file> <size>"
		<< endl << "          <probe radius> <gap between surface and skin> <radius type>"
		<< endl ;
}

bool populateSAS(int argc, char* argv[])
{
	if(argc != 8)
	{
		usagePopulateSAS();
		return false;
	}
	string input		= string(argv[2]);
	string output		= string(argv[3]);
	int size 		= atoi(argv[4]);
	double probeRadius	= atof(argv[5]);
	double floatingBandRadius = atof(argv[6]);
	int radiusInt		= atoi(argv[7]);
	PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType;
	if(!PDBParser::GroupOfAtoms::intToRadiusType(&radiusType, radiusInt))
	{
		return false;
	}

	return BlurMapsDataManager::populateSAS(input.c_str(), output.c_str(), size, size, size, probeRadius, floatingBandRadius, radiusType);
}


// Populate SAS using mesh
void usagePopulateSASUsingMesh()
{
	cout	<< endl << "Usage: MolSurf -populateSASUsingMesh <pqr input file> <rawn input file> <xyz output file>"
		<< endl << "          <probe radius> <gap between surface and skin> <minimum gap between centers of skin atoms> <radius type>"
		<< endl ;
}

/*
gap between surface and skin = see the definition of floating skin layer in the F2Dock paper
minimum gap between centers of skin atoms = used to control the density of the skin atoms sampled. This is useful to keep the number of skin atoms almost constant, no matter how dense the triangulation is.
radius type = see PDBparser::GroupOfAtoms for details
*/

bool populateSASUsingMesh(int argc, char* argv[])
{
	if(argc != 9)
	{
		usagePopulateSASUsingMesh();
		return false;
	}
	string pqr		= string(argv[2]);
	string surface		= string(argv[3]);
	string output		= string(argv[4]);
	double probeRadius	= atof(argv[5]);
	double floatingBandRadius = atof(argv[6]);
	double clashDistance = atof(argv[7]);
	int radiusInt		= atoi(argv[8]);
	PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType;
	if(!PDBParser::GroupOfAtoms::intToRadiusType(&radiusType, radiusInt))
	{
		return false;
	}

	return BlurMapsDataManager::populateSASUsingMesh(pqr, surface, output, probeRadius, floatingBandRadius, clashDistance, radiusType);
}




// Generate F2d
void usageGenerateF2d()
{
	cout	<< endl << "Usage: MolSurf -generateF2d <pqr input file> <xyz file> <output file> <int receptor=1/ligand=0> <radius type> <skin atom radius>"
		<< endl ;
}

bool generateF2d(int argc, char* argv[])
{
	if(argc < 7 || argc > 8)
	{
		usageGenerateF2d();
		return false;
	}

	string pqr		= string(argv[2]);
	string xyz		= string(argv[3]);
	string output		= string(argv[4]);

	int recep 		= atoi(argv[5]);
	bool receptor = false;
	if(recep) receptor = true;

	int radiusInt		= atoi(argv[6]);
	PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType;
	if(!PDBParser::GroupOfAtoms::intToRadiusType(&radiusType, radiusInt))
	{
		return false;
	}

	double skinAtomRad = 1.1;
	if(argc == 8)
		skinAtomRad	= atof(argv[7]);

	return BlurMapsDataManager::generateF2d(pqr, xyz, output, receptor, radiusType, skinAtomRad);
}


// Generating rmsd file for docking
void usageGetInterfaceAtoms()
{
	cout	<< endl << "Usage: MolSurf -getInterfaceAtoms <receptor pqr file> <ligand pqr file> <txt output file>"
		<< endl << "          <interface width> <string atomType = {all / backbone / any combination of {C,CA,N,O} (comma separated) }>"
		<< endl ;
}

/*
InterfaceWidth: if an atom of the ligand is within interfaceWidth of any atom of the receptor, then the ligand atom is on the interface

atomType: "all" = all interface atoms are printed in the rmsd file
          "backbone" = Only C, CA and N on the interface are printed
          comma separated list = Only interface atoms whose type is in the list is printed
*/

bool getInterfaceAtoms(int argc, char* argv[])
{
	if(argc != 7)
	{
		usageGetInterfaceAtoms();
		return false;
	}
	string pqr1		= string(argv[2]);
	string pqr2		= string(argv[3]);
	string output		= string(argv[4]);
	double interfaceWidth	= atof(argv[5]);
	string atomTypes	= string(argv[6]);
	int radiusInt		= 0;
	PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType;
	if(!PDBParser::GroupOfAtoms::intToRadiusType(&radiusType, radiusInt))
	{
		return false;
	}

	return BlurMapsDataManager::getInterfaceAtoms(pqr1, pqr2, output, interfaceWidth, atomTypes, radiusType);
}



// Generating rmsd file for docking
void usageGetInterfaceSurface()
{
	cout << endl << "Usage: MolSurf -getInterfaceSurface <receptor rawn/rawnc file> <ligand rawn/rawnc file> <interface width> <int coloredSurface (0 = rawn, 1 = rawnc)>" << endl ;
}

bool getInterfaceSurface(int argc, char* argv[])
{
	if(argc != 6)
	{
		usageGetInterfaceSurface();
		return false;
	}
	string surf1		= string(argv[2]);
	string surf2		= string(argv[3]);
	double interfaceWidth	= atof(argv[4]);
	int isColored		= atoi(argv[5]);
	bool isColoredBool = (isColored != 0);

	double area1, area2;

	return BlurMapsDataManager::getInterfaceSurface(surf1, surf2, interfaceWidth, isColoredBool, area1, area2);

	return false;
}



void usageColorByCharge()
{
	cout	<< endl << "Usage: MolSurf -colorByCharge <pqr file> <surface(RAWN) file> <output surface file(RAWNC)>"
		<< endl << "          <mode = {0 = Gaussian, 1 = average}> <distance cutoff> <blobbyness (optional)>"
		<< endl ;
}


bool colorByCharge(int argc, char* argv[])
{
	if(argc != 8 && argc != 7)
	{
		usageColorByCharge();
		return false;
	}
	string pqrf		= string(argv[2]);
	string sf		= string(argv[3]);
	string csf		= string(argv[4]);
	int mode		= atoi(argv[5]);
	double distanceCutoff	= atof(argv[6]);
	double blobbyness = 0.0;
	if(argc == 8) 
		blobbyness = atof(argv[7]);

	ColorSurfaceByCharge csbc(pqrf, sf, csf, mode, distanceCutoff, blobbyness);

	return csbc.colorify();
}



// Quality improve
void usageQualityImprove()
{
	cout	<< endl << "Usage: MolSurf -qualityImprove <raw input file> <raw output file> [parameter]"
		<< endl << "    Improve the mesh quality using geometric flow."
		<< endl ;
}

bool qualityImprove(int argc, char* argv[])
{
  double param = 0.011;
  if(argc != 4 && argc != 5)
    {
      usageQualityImprove();
      return false;
    }

  if (argc == 5) {
    param = atof(argv[4]);
  }

  string input	= string(argv[2]);
  string output	= string(argv[3]);
  LBIE::LBIE_Mesher().run(input.c_str(), output.c_str(), param);
  return true;
}

void usageRotatePDB() {
  cout	<< endl << "Usage: MolSurf -rotatePDB <pdb input file> <transoformation file> <pdb output file>"
	<< endl << "    Rotate a molecule."
	<< endl ;

}



bool rotatePDB(int argc, char * argv[]) {

  if (argc < 5) {
    usageRotatePDB();
    return false;
  }

  PDBParser::GroupOfAtoms* molecule = GOALoader().loadFile(argv[2]);
  if(!molecule) {
    cout << "Failed to load PDB." << endl;
    return false;
  }
  
  double trans[16];

  ifstream fin(argv[3]);
  for (int i=0; i<16; i++) {
    fin >> trans[i];
  }


  CCVOpenGLMath::Matrix transformation(trans[0], trans[1], trans[2], trans[3],
				       trans[4], trans[5], trans[6], trans[7],
				       trans[8], trans[9], trans[10], trans[11],
				       trans[12], trans[13], trans[14], trans[5]);
  GOALoader().saveFile(argv[4],molecule, 0, &transformation);

  return true;
}

// Call SDF
void usageSignedDistanceFunction()
{
	cout	<< endl << "Usage: MolSurf -signedDistanceFunction <raw input file> <rawiv output file> <size>"
		<< endl << "    Signed distance function using multi_sdf."
		<< endl ;
}

bool signedDistanceFunction(int argc, char* argv[])
{
	if(argc != 5)
	{
		usageSignedDistanceFunction();
		return false;
	}
	string input	= string(argv[2]);
	string output	= string(argv[3]);
	int size	= atoi(argv[4]);

	FaceVertSet3D fvs = FaceVertSet3D();
	Geom3DParser().ParseRawFile(fvs, input.c_str());
	int dim[3] = {size, size, size};
	DistanceTransform dt(fvs, dim, 2);
	dt.transform();
	dt.writeRawiv(output.c_str());
	return true;
}

// Call SDF
void usageSignedDistanceFunction2()
{
	cout	<< endl << "Usage: MolSurf -signedDistanceFunction2 <raw input file> <rawiv output file> <size x> <size y> <size z> "
		<< endl << "                <bbox x> <bbox y> <bbox z> <span x> <span y> <span z>" 
		<< endl << "    Signed distance function using multi_sdf."
		<< endl ;
}

bool signedDistanceFunction2(int argc, char* argv[])
{
	if(argc < 13)
	{
		usageSignedDistanceFunction2();
		return false;
	}
	string input	= string(argv[2]);
	string output	= string(argv[3]);

	int dim[3];
	float bbox[3];
	float spacing[3];
	dim[0] = atoi(argv[4]);
	dim[1] = atoi(argv[5]);
	dim[2] = atoi(argv[6]);
	bbox[0] = atof(argv[7]);
	bbox[1] = atof(argv[8]);
	bbox[2] = atof(argv[9]);
	spacing[0] = atof(argv[10]);	
	spacing[1] = atof(argv[11]);	
	spacing[2] = atof(argv[12]);	

	FaceVertSet3D fvs = FaceVertSet3D();
	Geom3DParser().ParseRawFile(fvs, input.c_str());

	DistanceTransform dt(fvs, dim, bbox, spacing);
	dt.transform();
	dt.writeRawiv(output.c_str());
	return true;
}

// Get surface atoms
void usageSurfaceAtoms()
{
	cout	<< endl << "Usage: MolSurf -surfaceAtoms <pdb input file> <output file>"
		<< endl << "    Compute which atoms of a molecule touch the molecular surface."
		<< endl ;
}

bool surfaceAtoms(int argc, char* argv[])
{
	if(argc != 4)
	{
		usageSurfaceAtoms();
		return false;
	}
	string input	= string(argv[2]);
	string output	= string(argv[3]);

	PDBParser::GroupOfAtoms* molecule = GOALoader().loadFile(input.c_str());
	if(!molecule)
	{
		return false;
	}
	PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType = PDBParser::GroupOfAtoms::VDW_RADIUS;
	SimpleVolumeData* sData = 0;
	sData = BlurMapsDataManager::getVolume(molecule, "", 128, 128, 128, PDBParser::GroupOfAtoms::ELECTRON_DENSITY, false, -1.0, PDBParser::GroupOfAtoms::ATOM, 0, 0, radiusType, PDBParser::ATOM_TYPE);
	if(!sData)
	{
		return false;
	}
	Geometry* largestComponent = SimpleVolumeDataIsocontourer::getIsocontour(sData, 1.0);
	
	// find surface atoms of molecule as those close to the largest component
	vector<PDBParser::Atom*> atomList;
	PDBParser::CollectionData* collectionData = 0;
	if(molecule->type == PDBParser::COLLECTION_TYPE)
	{
		collectionData = molecule->m_CollectionData;
	}
	FlattenGOA(molecule, atomList, collectionData, 0, 0, 0, radiusType, PDBParser::ATOM_TYPE, false);
	bool* boundaryAtom = new bool[atomList.size()];

	for(int i=0; i<atomList.size(); i++)
	{
		boundaryAtom[i] = false;
		for(int v=0; v<largestComponent->m_NumTriVerts; v++)
		{

			double dist = 100;
			if(isCloseTo(atomList[i]->m_Position[0],
					   atomList[i]->m_Position[1],
					   atomList[i]->m_Position[2],
					   atomList[i]->getRadius() + 1.0,
					   largestComponent->m_TriVerts[v*3+0],
					   largestComponent->m_TriVerts[v*3+1],
					   largestComponent->m_TriVerts[v*3+2],
					   &dist
					  ))

			{
				boundaryAtom[i] = true;
				break;
			}
		}
	}
	delete molecule;
	delete largestComponent;

	int numBdyAtoms = 0;
	// the boundary atoms file takes in surface area exposed for the given atom.

	FILE* fp	= fileWrite(input  + ".bdy");
	FILE* fpOut	= fileWrite(output + "_interior.pdb");
	FILE* fpSkin	= fileWrite(output + "_skin.pdb");
	FILE* fpIds	= fileWrite(output + "_ids.txt");

	for(int i=0; i<atomList.size(); i++)
	{
		if(boundaryAtom[i])
		{
			fprintf(fp, "100000.0\n");
			fprintf(fpIds, "%d\n", i+1);
			fprintf(fpSkin, "ATOM  11111  CA  GLU A1111    %8.3f%8.3f%8.3f%6.2f                \n", atomList[i]->m_Position[0], atomList[i]->m_Position[1], atomList[i]->m_Position[2], atomList[i]->getRadius());
			numBdyAtoms++;
		}
		else
		{
			fprintf(fp, "0.0\n");
			fprintf(fpOut, "ATOM  11111  N   GLU A1111    %8.3f%8.3f%8.3f%6.2f                \n", atomList[i]->m_Position[0], atomList[i]->m_Position[1], atomList[i]->m_Position[2], atomList[i]->getRadius());
		}
	}
	fclose(fp);
	fclose(fpOut);
	fclose(fpSkin);
	fclose(fpIds);
	
	// write docking input file
	fp = fileWrite(output.c_str());
	fprintf(fp, "%lu\n", atomList.size());
	for(int i=0; i<atomList.size(); i++)
	{
		if(boundaryAtom[i])
		{
			continue;
		}
		PDBParser::Atom* at = atomList[i];
		fprintf(fp, "I     %12.5lf %12.5lf %12.5lf %12.5f\n", at->m_Position[0], at->m_Position[1], at->m_Position[2], at->getCharge());
	}
	for(int i=0; i<atomList.size(); i++)
	{
		if(boundaryAtom[i])
		{
			PDBParser::Atom* at = atomList[i];
			fprintf(fp, "E     %12.5lf %12.5lf %12.5lf %12.5f\n", at->m_Position[0], at->m_Position[1], at->m_Position[2], at->getCharge());
		}
	}
	fclose(fp);

	int n = atomList.size();
	for(int i=0; i<n; i++)
	{
		delete atomList[i];
	}
	atomList.clear();
	if(boundaryAtom)
	{
		delete [] boundaryAtom;
		boundaryAtom = 0;
	}
	return true;
}

// Get surface from PDB
void usageSurfaceFromPDB()
{
	cout	<< endl << "Usage: MolSurf -surfaceFromPDB <pdb file> <surface file> <isovalue> <size>"
		<< endl << "    <blur?> <blobbiness> <radius type> <level> <optimize radii?>"
		<< endl << "    Undocumented get surface from PDB. Good luck."
		<< endl ;
}

bool surfaceFromPDB(int argc, char* argv[])
{
	if(argc != 11)
	{
		usageSurfaceFromPDB();
		return false;
	}
	string pdbFileName	= string(argv[2]);
	string surfaceFileName	= string(argv[3]);
	double isovalue 	= atof(argv[4]);
	int size 		= atoi(argv[5]);
	bool blurTheVolume 	= (strcmp(argv[ 6], "true") == 0);
	bool optimizeRadii 	= (strcmp(argv[10], "true") == 0);
	double blobbiness 	= atof(argv[7]);
	int radiusInt		= atoi(argv[8]);
	int level 		= atoi(argv[9]);
	PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType;
	if(!PDBParser::GroupOfAtoms::intToRadiusType(&radiusType, radiusInt))
	{
		return false;
	}
	Geometry* geometry = 0;
	if(blurTheVolume)
	{
		// create the volume
		PDBParser::GroupOfAtoms* molecule = GOALoader().loadFile(pdbFileName);
		if(!molecule)
		{
			return false;
		}
		SimpleVolumeData* sData = 0;
		if(optimizeRadii)
		{
			sData = BlurMapsDataManager::getMultiLevelBlurredVolume(molecule, "", size, size, size, PDBParser::GroupOfAtoms::ELECTRON_DENSITY, false, blobbiness, PDBParser::GroupOfAtoms::ATOM, 0, 0, radiusType, level);
		}
		else
		{
			sData = BlurMapsDataManager::getVolume(molecule, "",  size, size, size, PDBParser::GroupOfAtoms::ELECTRON_DENSITY, false, blobbiness, PDBParser::GroupOfAtoms::ATOM, 0, 0, radiusType, level);
		}
		delete molecule;
		if(!sData)
		{
			return false;
		}

		// create the isosurface
		geometry = SimpleVolumeDataIsocontourer::getIsocontour(sData, isovalue);
		if(geometry == 0)
		{
			delete sData;
			return false;
		}
		delete sData;
	}
	else
	{
		PDBParser::GroupOfAtoms* molecule = GOALoader().loadFile(pdbFileName);
		if(!molecule)
		{
			return false;
		}
		geometry = MolecularSurface::getMolecularSurface(molecule, size);
	}
	// save the isocontour
	saveAndDeleteGeometry(surfaceFileName, geometry);
	return true;
}



Geometry* surfaceUsingAdaptiveGrid(PDBParser::GroupOfAtoms* molecule, int size)
{
  return surfaceUsingAdaptiveGrid(molecule, size, 0.0);
}


Geometry* surfaceUsingAdaptiveGrid(PDBParser::GroupOfAtoms* molecule, int size, double extraRadius)
{
  //PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType;
  Geometry* geometry = MolecularSurface::getMolecularSurface(molecule, size, extraRadius);
  return geometry;
}
	



Geometry* surfaceUsingHLS(PDBParser::GroupOfAtoms* molecule, int size, int iterations, double weight)
{
	unsigned int dim[3] = {size, size, size};
	HLevelSet hls = HLevelSet();
	SimpleVolumeData* sData = hls.getHigherOrderLevelSetSurface(molecule, dim, iterations, weight);
	Geometry* geometry = 0;
	printf("Isovalue: %f \n", IsoValueofHLS);
	geometry = SimpleVolumeDataIsocontourer::getIsocontour(sData, (IsoValueofHLS+0.6));
	hls.computeGeometryNormal(geometry, dim);
	hls.fixOrientation(geometry); 	// make the normals point outwards
	//geometry->AllocateTriVertColors();
	geometry->ClearTriVertColors();  // arand, this is kind of a hack

	delete sData;
	return geometry; // arand, return value added...
}


// Generate a surface geometry using gaussian blur
Geometry* surfaceUsingGaussianBlur(PDBParser::GroupOfAtoms* molecule, double isovalue, int size, double blobbiness, int radiusInt, int level, bool optimizeRadii)
{
	PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType;
	if(!PDBParser::GroupOfAtoms::intToRadiusType(&radiusType, radiusInt))
	{
	  return false;
	}

	// sData
	SimpleVolumeData* sData = 0;
	if(optimizeRadii)
	{
		sData = BlurMapsDataManager::getMultiLevelBlurredVolume(molecule, "", size, size, size, PDBParser::GroupOfAtoms::ELECTRON_DENSITY, false, blobbiness, PDBParser::GroupOfAtoms::ATOM, 0, 0, radiusType, level);
	}
	else
	{
		sData = BlurMapsDataManager::getVolume(molecule, "", size, size, size, PDBParser::GroupOfAtoms::ELECTRON_DENSITY, false, blobbiness, PDBParser::GroupOfAtoms::ATOM, 0, 0, radiusType,level);
	}

	if(!sData)
	{
		return false;
	}
	//delete molecule;  // arand: don't delete or TexMol crashes...

	// create the isosurface
	Geometry* geometry = SimpleVolumeDataIsocontourer::getIsocontour(sData, isovalue);
	delete sData;

	return geometry;
}


// Generate a surface geometry using gaussian blur
Geometry* surfaceUsingGaussianBlur2(PDBParser::GroupOfAtoms* molecule, double isovalue, double resolution, double blobbiness, int radiusInt, int level)
{
	PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType;
	if(!PDBParser::GroupOfAtoms::intToRadiusType(&radiusType, radiusInt))
	{
	  return false;
	}

	// sData
	SimpleVolumeData* sData = 0;

	sData = BlurMapsDataManager::getVolumeFixedResolution(molecule, "", resolution, PDBParser::GroupOfAtoms::ELECTRON_DENSITY, false, blobbiness, PDBParser::GroupOfAtoms::ATOM, 0, 0, radiusType,level);


	if(!sData)
	{
		return false;
	}
	//delete molecule;  // arand: don't delete or TexMol crashes...

	// create the isosurface
	Geometry* geometry = SimpleVolumeDataIsocontourer::getIsocontour(sData, isovalue);
	delete sData;

	return geometry;
}





// Get surface using LBIE 
void usageSurfaceUsingLBIE()
{
	cout	<< endl << "Usage: MolSurf -surfaceUsingLBIE <meshtype> <isovalue> <rawiv input file> <raw output file>"
		<< endl << "    Generate tri/quad/tet/hex mesh from volumetric data using the LBIE library."
		<< endl ;
}

bool surfaceUsingLBIE(int argc, char* argv[])
{
	if(argc != 6)
	{
		usageSurfaceUsingLBIE();
		return false;
	}

	int isovalue	= atof(argv[3]);	
	string input	= string(argv[4]);
	string output	= string(argv[5]);

	MESHTYPE meshType = INVALIDMESH;
	if(strstr(argv[2], "tri"  )) meshType = TRI;
	if(strstr(argv[2], "tet"  )) meshType = TET;
	if(strstr(argv[2], "tri2" )) meshType = TRI2;
	if(strstr(argv[2], "tet2" )) meshType = TET2;
	if(strstr(argv[2], "quad" )) meshType = QUAD;
	if(strstr(argv[2], "hex"  )) meshType = HEX;
	if(strstr(argv[2], "nurbs")) meshType = NURBS;
	if(meshType == INVALIDMESH)
	{
		cerr << "WARNING:  Invalid meshtype, assuming tri." << endl;
		meshType = TRI;
	}

	LBIE::LBIE_Mesher lbie_mesher;
	lbie_mesher.setMeshType((MESHTYPE) meshType);
	lbie_mesher.run(input.c_str(),output.c_str(),-isovalue);
	return true;
}

// Get volume
void usageVolume()
{
	cout	<< endl << "Usage: MolSurf -volume <raw input file> <txt output file>"
		<< endl << "    Compute the volume bounded by a triangular surface mesh."
		<< endl ;
}

bool volume(int argc, char* argv[])
{
	// First read in the arguments
	if(argc != 4)	// If the wrong number of arguments, print usage
	{
		usageVolume();
		return false;
	}
	string input	= string(argv[2]);
	FILE* output	= fileWrite(argv[3]);

	// read in the isocontour
	Geometry* geometry = GeometryLoader().loadFile(input.c_str());
	double myVolume = geometry->getVolume();
	delete geometry;

	// calculate surface area and print it out
	printf(" %lf ", myVolume);
	fprintf(output, "%lf\n", myVolume);
	fclose(output);
	return true;
}

// Get volume using HLS
void usageVolumeUsingHLS()
{
	cout	<< endl << "Usage: MolSurf -volumeUsingHLS <pdb/pqr input file> <rawiv output file> <size> [nIterations weight]"
		<< endl << "    Generate volumetric level-set function."
		<< endl ;
}

bool volumeUsingHLS(int argc, char* argv[])
{
	if(argc != 7 && argc != 5)
	{
		usageVolumeUsingHLS();
		false;
	}
	string input	= string(argv[2]);
	string output	= string(argv[3]);
	int size 	= atoi(argv[4]);
	int nIterations	= 3;
	double weight	= 0.0001;
	if(argc == 7)
	{
		nIterations = atoi(argv[5]);
		weight	= atof(argv[6]);
	}
	PDBParser::GroupOfAtoms* molecule = GOALoader().loadFile(input.c_str());
	unsigned int dim[3] = {size, size, size};
	SimpleVolumeData* sData = HLevelSet().getHigherOrderLevelSetSurface(molecule, dim, nIterations, weight);
	VolumeLoader().saveFile(output.c_str(), sData);
	return true;
}




// Get volume of hydrophobicity

void usageGetVolumeWithHydrophobicity()
{
	cout	<< endl << "Usage: MolSurf -getVolumeWithHydrophobicity <pdb/pqr input file> <rawiv output file> <dim1> <dim2> <dim3> <blobbiness> [<bbxmin> <bbymin> <bbzmin> <bbxmax> <bbymax> <bbzmax> ]"
		<< endl ;
}

bool getVolumeWithHydrophobicity(int argc, char* argv[])
{
	if(argc != 8 && argc != 14)
	{
		usageGetVolumeWithHydrophobicity();
		return false;
	}
	string input	= string(argv[2]);
	string output	= string(argv[3]);
	int dim1 	= atoi(argv[4]);
	int dim2 	= atoi(argv[5]);
	int dim3 	= atoi(argv[6]);
	double blobby = atof(argv[7]);



	float bbmin[3], bbmax[3];
	if(argc ==14)
	{
		bbmin[0] = atof(argv[8]);
		bbmin[1] = atof(argv[9]);
		bbmin[2] = atof(argv[10]);

		bbmax[0] = atof(argv[11]);
		bbmax[1] = atof(argv[12]);
		bbmax[2] = atof(argv[13]);
	}

	SimpleVolumeData* sData = 0;

	PDBParser::GroupOfAtoms* molecule = GOALoader().loadFile(input.c_str());

	PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType;
	if(!PDBParser::GroupOfAtoms::intToRadiusType(&radiusType, 0))
	{
		return false;
	}
	
	if(argc !=14)
		sData = BlurMapsDataManager::getVolume(molecule, output.c_str(), dim1, dim2, dim3, PDBParser::GroupOfAtoms::PER_ATOM_HYDROPHOBICITY, false, blobby, PDBParser::GroupOfAtoms::ATOM, NULL, 0, radiusType, PDBParser::ATOM_TYPE);
	else
		sData = BlurMapsDataManager::getVolume(molecule, output.c_str(), dim1, dim2, dim3, bbmin, bbmax, PDBParser::GroupOfAtoms::PER_ATOM_HYDROPHOBICITY, false, blobby, PDBParser::GroupOfAtoms::ATOM, NULL, 0, radiusType, PDBParser::ATOM_TYPE);
	if(sData == 0)
	{
		return false;
	}

	return true;
}

// Get volume from electron
void usageGetVolumeWithElectron()
{
	cout	<< endl << "Usage: MolSurf -getVolumeWithElectron <pdb/pqr input file> <rawiv output file> <dim1> <dim2> <dim3> <blobbiness>"
		<< endl ;
}

bool getVolumeWithElectron(int argc, char* argv[])
{
	if(argc != 8)
	{
		usageGetVolumeWithElectron();
		return false;
	}
	string input	= string(argv[2]);
	string output	= string(argv[3]);
	int dim1 	= atoi(argv[4]);
	int dim2 	= atoi(argv[5]);
	int dim3 	= atoi(argv[6]);
	double blobby = atof(argv[7]);

	SimpleVolumeData* sData = 0;

	PDBParser::GroupOfAtoms* molecule = GOALoader().loadFile(input.c_str());

	PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType;
	if(!PDBParser::GroupOfAtoms::intToRadiusType(&radiusType, 0))
	{
		return false;
	}
	
	sData = BlurMapsDataManager::getVolume(molecule, output.c_str(), dim1, dim2, dim3, PDBParser::GroupOfAtoms::ELECTRON_DENSITY, false, blobby, PDBParser::GroupOfAtoms::ATOM, NULL, 0, radiusType, PDBParser::ATOM_TYPE);

	if(sData == 0)
	{
		return false;
	}

	return true;
}



// Get volume from charge density
void usageGetVolumeWithChargeDensity()
{
	cout	<< endl << "Usage: MolSurf -getVolumeWithChargeDensity <pqr input file> <rawiv output file> <dim1> <dim2> <dim3> <blobbiness> [<bbxmin> <bbymin> <bbzmin> <bbxmax> <bbymax> <bbzmax> ]"
		<< endl ;
}

bool getVolumeWithChargeDensity(int argc, char* argv[])
{
	if(argc != 8 && argc != 14)
	{
		usageGetVolumeWithChargeDensity();
		return false;
	}
	string input	= string(argv[2]);
	string output	= string(argv[3]);
	int dim1 	= atoi(argv[4]);
	int dim2 	= atoi(argv[5]);
	int dim3 	= atoi(argv[6]);
	double blobby = atof(argv[7]);
	
	float bbmin[3], bbmax[3];
	if(argc ==14)
	{
		bbmin[0] = atof(argv[8]);
		bbmin[1] = atof(argv[9]);
		bbmin[2] = atof(argv[10]);

		bbmax[0] = atof(argv[11]);
		bbmax[1] = atof(argv[12]);
		bbmax[2] = atof(argv[13]);
	}
	SimpleVolumeData* sData = 0;

	PDBParser::GroupOfAtoms* molecule = GOALoader().loadFile(input.c_str());

	PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType;
	if(!PDBParser::GroupOfAtoms::intToRadiusType(&radiusType, 0))
	{
		return false;
	}
	
	if(argc !=14)
		sData = BlurMapsDataManager::getVolume(molecule, output.c_str(), dim1, dim2, dim3, PDBParser::GroupOfAtoms::ELECTROSTATIC_POTENTIAL, false, blobby, PDBParser::GroupOfAtoms::ATOM, NULL, 0, radiusType, PDBParser::ATOM_TYPE);
	else
		sData = BlurMapsDataManager::getVolume(molecule, output.c_str(), dim1, dim2, dim3, bbmin, bbmax, PDBParser::GroupOfAtoms::ELECTROSTATIC_POTENTIAL, false, blobby, PDBParser::GroupOfAtoms::ATOM, NULL, 0, radiusType, PDBParser::ATOM_TYPE);
		

	if(sData == 0)
	{
		return false;
	}

	return true;
}

