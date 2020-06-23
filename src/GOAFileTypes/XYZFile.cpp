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
#include <GOAFileTypes/XYZFile.h>
#include <math.h>
#include <PDBParser/Atom.h>
#include <PDBParser/CollectionData.h>
#include <PDBParser/GroupOfAtoms.h>
#include <PDBParser/FlattenGOA.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>

XYZFile XYZFile::ms_XYZFileRepresentative;

XYZFile::XYZFile()
{
}

XYZFile::~XYZFile()
{
}

PDBParser::GroupOfAtoms* XYZFile::loadFile(const string& fileName, bool deleteWater)
{
	FILE* fp=fopen(fileName.c_str(), "r");
	// check to make sure the file exists
	if(!fp)
	{
		printf("Error: could not open file %s\n", fileName.c_str());
		return 0;
	}
	PDBParser::GroupOfAtoms* molecule = PDBParser::GroupOfAtoms::getEmptyHierarchy();
	PDBParser::GroupOfAtoms* firstResidue = molecule->m_SubGroups[0]->m_SubGroups[0]->m_SubGroups[0]->m_SubGroups[0];
	int numAtoms = 0;
	while(true)
	{
		double x=0, y=0, z=0, radius=1;
		if(fscanf(fp, "%lf %lf %lf\n", &x, &y, &z) != 3)
		{
			break;
		}
		PDBParser::Atom* atom = new PDBParser::Atom();
		atom->m_Position[0] = x;
		atom->m_Position[1] = y;
		atom->m_Position[2] = z;
		atom->setRadius(radius);
		firstResidue->addAtom(
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
		delete atom;
		numAtoms++;
	}
	double xmin=0, ymin=0, zmin=0, xmax=1, ymax=1, zmax=1;
	molecule->getMinMax(& xmin, & ymin, & zmin, & xmax, & ymax, & zmax);
	molecule->buildHierarchicalInformation((xmin + xmax)/ 2.0, (ymin + ymax) /2.0, (zmin + zmax) /2.0);
	if(numAtoms == 0)
	{
		delete molecule;
		molecule = 0;
	}
	fclose(fp);
	return molecule;
}

bool XYZFile::checkType(const string& fileName)
{
	return false;
}

bool XYZFile::saveFile(PDBParser::GroupOfAtoms* molecule, const string& fileName, unsigned int level, CCVOpenGLMath::Matrix* transformation)
{
	FILE* fp = fopen(fileName.c_str(), "w");
	// failed to open the file? bail.
	if(!fp)
	{
		printf("Could not open file %s\n", fileName.c_str());
		return false;
	}
	// make 2 loops. First loop, get number of atoms. second loop, write them out.
	vector<PDBParser::Atom*> m_AtomList;
	PDBParser::CollectionData* collectionData = 0;
	if(molecule->type == PDBParser::COLLECTION_TYPE)
	{
		collectionData = molecule->m_CollectionData;
	}
	FlattenGOA(molecule, m_AtomList, collectionData, 0, 0, 0, PDBParser::GroupOfAtoms::VDW_RADIUS, level, false);
	int n = m_AtomList.size();
	{
		int i;
		for(i=0; i<n; i++)
		{
			fprintf(fp, "%f %f %f\n", m_AtomList[i]->m_Position[0], m_AtomList[i]->m_Position[1], m_AtomList[i]->m_Position[2]);
		}
	}
	// clean up
	fclose(fp);
	return true;
}

GOAFileType* XYZFile::getRepresentative()
{
	return &ms_XYZFileRepresentative;
}
