#include <Utility/utility.h>
#include <Blurmaps/MolecularSurface.h>
#include <libCG/OptimizerLib/Optimizer/optimizer.h>
#include <libCG/ARM/ARM.h>
#include <libCG/CoarseGrain/geometry.h>
#include <libCG/CoarseGrain/surface.h>
#include <PDBParser/FlattenGOA.h>

#include <fstream>

using namespace PDBParser;
using namespace GEOMETRY;

GroupOfAtoms*  ARM::read_pdb(const char* fname)
{
	GOALoader* gLoader = new GOALoader();
	PDBParser::GroupOfAtoms* molecule = gLoader->loadFile(fname);
	delete gLoader;
	if(!molecule)
	{
		return 0;
	}
	return molecule;
}

void ARM::getAtomList(std::vector<PDBParser::Atom*>& atomList, const char* fname)
{
	GroupOfAtoms* mol=NULL;
	mol = read_pdb(fname);
	printf("after read %s\n", fname);
	CollectionData* collectionData = 0;
	if(mol->type == COLLECTION_TYPE)
	{
		collectionData = mol->m_CollectionData;
	}
	GroupOfAtoms::RADIUS_TYPE radiusType = GroupOfAtoms::VDW_RADIUS;
	PDBParser::FlattenGOA(mol, atomList, collectionData, 0, 0, 0, radiusType,
									ATOM_TYPE, false);
	unsigned int noatoms = atomList.size();
	printf("number of atoms:%d\n",noatoms);
}

ARM::ARM(char* pqrfiles[ARMLevel])
{
	getAtomList(arm[0],pqrfiles[0]);
	getAtomList(arm[1],pqrfiles[1]);
}

ARM::~ARM()
{
}

void ARM::setHierarchicalInformation(char* hinfofile)
{
  // arand: changed to c++ stype input
  //        since getline not supported on mac
  // WARNING: UNTESTED!!!!
  //FILE* fp = fopen(hinfofile,"r");

	ifstream fin(hinfofile);

	char* token = NULL;
	char* line = NULL;
	char* delim = (char*)" \t\n";
	string mystr;
	//while((read = getline(&line,&len,fp)) != -1)
	while(getline(fin,mystr))
	{
	  line = (char *)mystr.c_str();
		token = strtok(line,delim);
		if(token == NULL)
		{
			break;
		}
		if(!strcmp(token,"Bead"))
		{
			token = strtok(NULL,delim);
			int index = atoi(token);
			//printf("index:%d\n",index);
			BeadStructure bead;
			bead.index = index;
			//if((read = getline(&line,&len,fp)) != -1)
			if(getline(fin,mystr))
			{
			  line = (char *) mystr.c_str();
				token = strtok(line,delim);
				int numofchilds = 0;
				while(token != NULL)
				{
					bead.childs[numofchilds] = atoi(token);
					numofchilds++;
					token = strtok(NULL,delim);
				}
				bead.numchilds = numofchilds;
			}
			armHierarchy.push_back(bead);
		}
	}
}

void copyGeometryToSurface(Geometry* geo, Surface* surface)
{
	double maxx, maxy, maxz, minx, miny, minz;
	maxx = -100000.0;
	maxy = -100000.0;
	maxz = -100000.0;
	minx = 100000.0;
	miny = 100000.0;
	minz = 100000.0;
	surface->numbpts = geo->m_NumTriVerts;
	surface->numbtris = geo->m_NumTris;
	for(int i=0; i<geo->m_NumTriVerts; i++)
	{
		float x = geo->m_TriVerts[i*3];
		float y = geo->m_TriVerts[i*3+1];
		float z = geo->m_TriVerts[i*3+2];
		Vertex* point = new Vertex();
		point->Point[0] = x;
		point->Point[1] = y;
		point->Point[2] = z;
		if(x < minx)
		{
			minx = x;
		}
		if(y < miny)
		{
			miny = y;
		}
		if(z < minz)
		{
			minz = z;
		}
		if(x > maxx)
		{
			maxx = x;
		}
		if(y > maxy)
		{
			maxy = y;
		}
		if(z > maxz)
		{
			maxz = z;
		}
		if(geo->m_bTriSmoothNormalsReady)
		{
			float nx = geo->m_TriVertNormals[i*3];
			float ny = geo->m_TriVertNormals[i*3+1];
			float nz = geo->m_TriVertNormals[i*3+2];
			float w = sqrt(nx*nx + ny*ny + nz*nz);
			point->Normal[0] = nx/w;
			point->Normal[1] = ny/w;
			point->Normal[2] = nz/w;
		}
		else
		{
			point->Normal[0] = 0.0;
			point->Normal[1] = 0.0;
			point->Normal[2] = 0.0;
		}
		surface->addPoint(point);
	}
	for(int i=0; i<geo->m_NumTris; i++)
	{
		int i1 = geo->m_Tris[i*3];
		int i2 = geo->m_Tris[i*3+1];
		int i3 = geo->m_Tris[i*3+2];
		Face* facet = new Face();
		facet->Index[0] = i1;
		facet->Index[1] = i2;
		facet->Index[2] = i3;
		facet->AdjTri[0] = -1;
		facet->AdjTri[1] = -1;
		facet->AdjTri[2] = -1;
		facet->whichnml[0] = -1;
		facet->whichnml[1] = -1;
		facet->whichnml[2] = -1;
		surface->addFacet(facet);
	}
	const double layer = 1.0;
	surface->min[0] = minx-layer;
	surface->min[1] = miny-layer;
	surface->min[2] = minz-layer;
	surface->max[0] = maxx+layer;
	surface->max[1] = maxy+layer;
	surface->max[2] = maxz+layer;
	for(int i = 0; i < 3; i++)
	{
		surface->center[i] = (surface->min[i] + surface->max[i])*0.5;
	}
}

void ARM::getAtomListOfBead(int index, std::vector<PDBParser::Atom*>& atomList)
{
	for(int i=0; i<armHierarchy[index].numchilds; i++)
	{
		int id = armHierarchy[index].childs[i];
		PDBParser::Atom* atom = arm[0][id];
		atomList.push_back(atom);
	}
}

void ARM::getOptimizedCG(std::vector<PDBParser::Atom*>& cg, int size, int iteration, char* outfilename)
{
	std::vector<PDBParser::Atom*> atomList_bead;
	float totalcharge = 0.0;
	float totalcharge2 = 0.0;
	for(int i=0; i< armHierarchy.size(); i++)
	{
		Geometry* geometry = getSurfaceOfBead(i,size);
		Surface* surface = new Surface();
		copyGeometryToSurface(geometry,surface);
		printf("surface has %lu vertices.\n",surface->m_Vertices.size());
		int id = armHierarchy[i].index;
		printf("Bead index is:%d\n",id);
		PDBParser::Atom* bead = arm[1].at(id);
		getAtomListOfBead(i,atomList_bead);
		Optimizer* optm = new Optimizer(atomList_bead, bead, surface);
		optm->optimizeCenterAndRadius(iteration);
		printf("before Optimized: %f %f %f %f\n", arm[1][id]->m_Position[0], arm[1][id]->m_Position[1],arm[1][id]->m_Position[2], arm[1][id]->getRadius());
		bead->setPosition(optm->adata[0].center[0], optm->adata[0].center[1], optm->adata[0].center[2]);
		bead->setRadius(optm->adata[0].radius);
		printf("Optimized: %f %f %f %f\n", arm[1][id]->m_Position[0], arm[1][id]->m_Position[1],arm[1][id]->m_Position[2], arm[1][id]->getRadius());
		float blob =-2.3442;
		blob = -2.0;
		blob = -1.6;
		float charge2 = 0.0;
		float charge = 0.0;
		for(int k=0; k< atomList_bead.size(); k++)
		{

			charge += atomList_bead[k]->getCharge()*exp(blob*(((atomList_bead[k]->m_Position[0]-arm[1][id]->m_Position[0])*(atomList_bead[k]->m_Position[0]-arm[1][id]->m_Position[0]) + (atomList_bead[k]->m_Position[1]-arm[1][id]->m_Position[1])*(atomList_bead[k]->m_Position[1]-arm[1][id]->m_Position[1]) +
					  (atomList_bead[k]->m_Position[2]-arm[1][id]->m_Position[2])*(atomList_bead[k]->m_Position[2]-arm[1][id]->m_Position[2])) /(atomList_bead[k]->getRadius()*atomList_bead[k]->getRadius()) -1.0));
			charge2+= atomList_bead[k]->getCharge();
		}
		if(fabs(charge-charge2)>= 1.0)
		{
			cout <<"Big difference !"<<endl;
		}
		totalcharge += charge2;
		arm[1][id]->setCharge(charge);
		totalcharge2 += charge;
		atomList_bead.clear();
		delete geometry;
		delete surface;
		delete optm;
	}
	cout<<"Total charge: "<<totalcharge <<endl;
	cout<<"Coarse Grained model charge is : "<<totalcharge2 <<endl;
	float temp= (totalcharge-totalcharge2)/armHierarchy.size();
	for(int i=0; i< armHierarchy.size(); i++)
		//for(int i=0;i<5; i++)  for test
	{
		int id = armHierarchy[i].index;
		float tmpcharge = temp + arm[1][id]->getCharge();
		arm[1][id]->setCharge(tmpcharge);
	}
	printAtomList(arm[1],outfilename);
}

void ARM::printAtomList(std::vector<PDBParser::Atom*>& atomList, char* outfilename)
{
	FILE* fp = fopen(outfilename,"w");
	for(int i=0; i<atomList.size(); i++)
	{
		fprintf(fp, "ATOM%7d  C   %3s %c%4d    %8.3f%8.3f%8.3f%8.4f%7.4f\n", i+1,
				atomList[i]->getResidueName(),
				atomList[i]->getChain(),
				atomList[i]->getResidueNum(),
				atomList[i]->m_Position[0],
				atomList[i]->m_Position[1],
				atomList[i]->m_Position[2],
				atomList[i]->getCharge(),
				atomList[i]->getRadius());
	}
	fclose(fp);
}

Geometry* ARM::getSurfaceOfBead(int index, int size)
{
	PDBParser::GroupOfAtoms* molecule = PDBParser::GroupOfAtoms::getEmptyHierarchy();
	PDBParser::GroupOfAtoms* curResidue = molecule->m_SubGroups[0]->m_SubGroups[0]->m_SubGroups[0]->m_SubGroups[0];
	for(int i=0; i<armHierarchy[index].numchilds; i++)
	{
		int id = armHierarchy[index].childs[i];
		PDBParser::Atom* atom = arm[0][id];
		curResidue->addAtom(
			atom->getAtomNum(),
			atom->getElementIndex(),
			atom->getResidueNum(),
			(char*)atom->getName(),
			atom->getResidueId(),
			atom->getChain(),
			atom->m_Position[0],
			atom->m_Position[1],
			atom->m_Position[2],
			atom->getRadius(),
			atom->getCharge(),
			false);
	}
	Geometry* geometry = MolecularSurface::getMolecularSurface(molecule, size);
	curResidue = 0;
	delete molecule;
	return geometry;
}
