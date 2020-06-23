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
#include <CardinalBSpline/CardinalBSpline.h>
#include <PDBParser/FlattenGOA.h>

CardinalBSpline::CardinalBSpline()
{
}

CardinalBSpline::~CardinalBSpline()
{
}

bool CardinalBSpline::computeFunction(std::vector<PDBParser::Atom*> atomList, float* data, unsigned int* dim, float* minExtent, float* maxExtent, int order)
{
	float position[3];
	if(!data || !dim || !minExtent || !maxExtent)
	{
		return false;
	}
	for(unsigned int i = 0; i < dim[0]*dim[1]*dim[2]; i++)
	{
		data[i]=0.0f;
	}
	float dinteval[3];
	for(int i=0; i<3; i++)
	{
		dinteval[i]= (maxExtent[i]-minExtent[i])/(dim[i]-1);
	}
	int size = atomList.size();
	for(int i = 0; i < size; i++)
	{
		PDBParser::Atom* atom = atomList[i];
		if(!atom)
		{
			continue;
		}
		position[0] = atom->m_Position[0];
		position[1] = atom->m_Position[1];
		position[2] = atom->m_Position[2];
		float radius = atom->getRadius();
		float probeRadius = 0.0;
		radius += probeRadius;
		CardinalBSpline_One_Atom(position, radius, order, dim,  minExtent, dinteval, data);
	}
}

bool CardinalBSpline::getAtomListAndExtent(GroupOfAtoms* molecule, std::vector<PDBParser::Atom*> & atomList, float* minExtent, float* maxExtent)
{
	if(!molecule || !minExtent || !maxExtent)
	{
		return false;
	}
	CollectionData* collectionData = 0;
	if(molecule->type == COLLECTION_TYPE)
	{
		collectionData = molecule->m_CollectionData;
	}
	GroupOfAtoms::RADIUS_TYPE radiusType = GroupOfAtoms::VDW_RADIUS;
	FlattenGOA(molecule, atomList, collectionData, 0, 0, 0, radiusType, ATOM_TYPE, false);
	double probeRadius = 0.0;
	BlurMapsDataManager::getBoundingBox(atomList, minExtent, maxExtent, radiusType, probeRadius*2+1, 0);
	return true;
}

SimpleVolumeData* CardinalBSpline::getCardinalBSplineSurface(GroupOfAtoms* molecule, unsigned int* dim, int order)
{
	std::vector<PDBParser::Atom*> atomList;
	if(!getAtomListAndExtent(molecule, atomList, minExtent, maxExtent))
	{
		return 0;
	}
	float* data = new float[dim[0]*dim[1]*dim[2]];
	if(!computeFunction(atomList, data, dim, minExtent, maxExtent, order))
	{
		delete []data;
		return 0;
	}
	SimpleVolumeData* sData = new SimpleVolumeData(dim);
	sData->setDimensions(dim);
	sData->setNumberOfVariables(1);
	sData->setData(0,data);
	sData->setType(0,SimpleVolumeData::FLOAT);
	sData->setName(0,"CardinalBSpline");
	sData->setMinExtent(minExtent);
	sData->setMaxExtent(maxExtent);
	return sData;
}

void CardinalBSpline::CardinalBSpline_One_Atom(float* position, float radius, int order,  unsigned int* dim, float* minExtent, float* dinteval, float* data)
{
	float funct;
	float x, y, z;
	int nx, ny, nz;
	int space[3];
	nx = (int)((position[0]-minExtent[0])/dinteval[0]);
	ny = (int)((position[1]-minExtent[1])/dinteval[1]);
	nz = (int)((position[2]-minExtent[2])/dinteval[2]);
	space[0] = (int)((radius/dinteval[0]) + 1);
	space[1] = (int)((radius/dinteval[1]) + 1);
	space[2] = (int)((radius/dinteval[2]) + 1);
	int tmp=space[0];
	tmp = (space[1] > tmp)? space[1]:tmp;
	tmp = (space[2] > tmp)? space[2]:tmp;
	radius *= 0.8;
	for(int i= nx-space[0]-tmp; i< nx+space[0]+tmp; i++)
	{
		funct=0.0;
		if((i >= 0) &&(i < dim[0]))
		{
			x = minExtent[0] + i * dinteval[0] - position[0];
		}
		else
		{
			continue;
		}
		for(int j= ny-space[1]-tmp; j< ny+space[1]+tmp; j++)
		{
			if((j >= 0) &&(j < dim[1]))
			{
				y = minExtent[1] + j * dinteval[1] - position[1];
			}
			else
			{
				continue;
			}
			for(int k = nz -space[2]-tmp; k < nz +space[2]+tmp; k++)
			{
				if((k >=0) && (k < dim[2]))
				{
					z = minExtent[2] + k * dinteval[2] - position[2];
					funct = CardinalBBasis(order, x,radius)*CardinalBBasis(order, y,radius)*CardinalBBasis(order, z,radius);
					if(funct <= -0.0)
					{
						funct*= -1.0;
					}
					data[k* dim[0]*dim[1]+j* dim[0]+i] += funct;
				}
				else
				{
					continue;
				}
			}
		}
	}
}
