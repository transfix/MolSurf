#include <string>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <Utility/utility.h>
#include <PDBParser/Atom.h>
#include <PDBParser/CollectionData.h>
#include <PDBParser/ElementInformation.h>
#include <PDBParser/GroupOfAtoms.h>
#include <GOAFileTypes/GOALoader.h>
#include <PDBParser/FlattenGOA.h>
#include <PDBParser/ResidueInformation.h>

using namespace std;

namespace PDBParser
{
	class GroupOfAtoms;
	class Atom;
};

int main(int argc, char** argv)
{
	if(argc !=4)
	{
		cout<<"Usage: computeRMSD ref_pdb xformed_pdb atomlist\n";
		return -1;
	}

	int radiusInt = 0;
	PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType;
	if(!PDBParser::GroupOfAtoms::intToRadiusType(&radiusType, radiusInt))
	{
		return false;
	}

	PDBParser::GroupOfAtoms* molecule1 = 0;
	GOALoader* gLoader1 = new GOALoader();
	molecule1 = gLoader1->loadFile(argv[1]);
	delete gLoader1;
	if(!molecule1)
	{
		return 0;
	}

	PDBParser::GroupOfAtoms* molecule2 = 0;
	GOALoader* gLoader2 = new GOALoader();
	molecule2 = gLoader2->loadFile(argv[2]);
	delete gLoader2;
	if(!molecule2)
	{
		return 0;
	}

	vector<PDBParser::Atom*> m_AtomList1;
	vector<PDBParser::Atom*> m_AtomList2;

	PDBParser::CollectionData* collectionData1 = 0;
	if(molecule1->type == PDBParser::COLLECTION_TYPE)
	{
		collectionData1 = molecule1->m_CollectionData;
	}
	FlattenGOA(molecule1, m_AtomList1, collectionData1, 0, 0, 0, radiusType, PDBParser::ATOM_TYPE, false);

	PDBParser::CollectionData* collectionData2 = 0;
	if(molecule2->type == PDBParser::COLLECTION_TYPE)
	{
		collectionData2 = molecule2->m_CollectionData;
	}
	FlattenGOA(molecule2, m_AtomList2, collectionData2, 0, 0, 0, radiusType, PDBParser::ATOM_TYPE, false);

	vector<int> atomids;

	FILE* fp = fopen(argv[3], "rt");
	int numAtoms;
	int numRMSAtoms;
	double SD = 0.0;
	double distsq;

	double x,y,z;
	int id;
	double x1, x2, y1, y2, z1, z2;

	fscanf(fp, "%d", &numAtoms);
	fscanf(fp, "%d", &numRMSAtoms);

	cout<<numAtoms<<" "<<numRMSAtoms<<endl;

	for(int i=0; i<numRMSAtoms; i++)
	{
		fscanf(fp, "%d %lf %lf %lf", &id,&x,&y,&z);

		cout<<id<<" "<<x<<" "<<y<<" "<<z<<endl;

		PDBParser::Atom* a1 = m_AtomList1[id];
		PDBParser::Atom* a2 = m_AtomList2[id];

		x1 = a1->getX();
		y1 = a1->getY();
		z1 = a1->getZ();

		x2 = a2->getX();
		y2 = a2->getY();
		z2 = a2->getZ();
	
		distsq = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
		SD += distsq;
	}

	double MSD = SD/(double)numRMSAtoms;

	double RMSD = sqrt(MSD);

	cout<<argv[1]<<"  "<<argv[2]<<"  "<<numRMSAtoms<<"  "<<RMSD<<endl;
}
