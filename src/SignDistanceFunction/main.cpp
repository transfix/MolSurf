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
#include <Geometry/Geometry.h>
#include <SignDistanceFunction/common.h>
#include <SimpleVolumeData/SimpleVolumeData.h>

// Include order dependency
#include <SignDistanceFunction/sdfLib.h>

using namespace SDFLibrary;

void free_memory()
{
	SDFLibrary::listnode* temp;
	SDFLibrary::listnode* currNode;
	printf("starting memory de-allocation\n");
	for(int i = 0; i < SDFLibrary::size; i++)
	{
		for(int j = 0; j < SDFLibrary::size; j++)
		{
			for(int k = 0; k < SDFLibrary::size; k++)
			{
				currNode = SDFLibrary::sdf[i][j][k].tindex;
				while(currNode != NULL)
				{
					temp = currNode;
					currNode = currNode->next;
					free(temp);
				}
			}
			free(SDFLibrary::sdf[i][j]);
		}
		free(SDFLibrary::sdf[i]);
	}
	free(SDFLibrary::sdf);
	free(SDFLibrary::values);
	if(SDFLibrary::vertices != NULL)
	{
		free(SDFLibrary::vertices);
	}
	if(SDFLibrary::surface != NULL)
	{
		free(SDFLibrary::surface);
	}
	if(SDFLibrary::normals != NULL)
	{
		free(SDFLibrary::normals);
	}
	if(SDFLibrary::distances != NULL)
	{
		free(SDFLibrary::distances);
	}
	if(SDFLibrary::queues != NULL)
	{
		free(SDFLibrary::queues);
	}
	if(SDFLibrary::bverts != NULL)
	{
		free(SDFLibrary::bverts);
	}
	printf("Memory de-allocated successfully! \n");
}

bool SDFLibrary::setParameters(int Size, int isNormalFlip, float* mins, float* maxs)
{
	//First the default values.
	SDFLibrary::init_all_vars();
	//Then, assign the actual input values.
	SDFLibrary::size = Size;
	SDFLibrary::flipNormals = isNormalFlip;
	SDFLibrary::minext[0] = mins[0];
	SDFLibrary::minext[1] = mins[1];
	SDFLibrary::minext[2] = mins[2];
	SDFLibrary::maxext[0] = maxs[0];
	SDFLibrary::maxext[1] = maxs[1];
	SDFLibrary::maxext[2] = maxs[2];
	SDFLibrary::span[0] = (maxs[0]-mins[0])/(SDFLibrary::size);
	SDFLibrary::span[1] = (maxs[1]-mins[1])/(SDFLibrary::size);
	SDFLibrary::span[2] = (maxs[2]-mins[2])/(SDFLibrary::size);

	// arand: is there any reason that the size must be a power of two?
	if((Size!=16) && (Size!=32) &&(Size!=64) && (Size!=128) && (Size!=256) &&(Size!=512) &&(Size!=1024))
	{
		printf("size is incorrect\n");
		return false;
	}
	return true;
}

float* SDFLibrary::computeSDF(int nverts, float* verts, int ntris, int* tris)
{
	int numb;
	float* sdfValues =NULL;
	float isoval;
	//Set up the volume grid
	if(!initSDF())
	{
		return 0;
	}
	//Read in the Geometry
	readGeom(nverts, verts, ntris, tris);
	//Setup the Octree
	adjustData();
	//Compute the SDF
	compute();
	//Return the SDF
	numb = (SDFLibrary::size+1)*(SDFLibrary::size+1)*(SDFLibrary::size+1);
	sdfValues = (float*)(malloc(sizeof(float)*(numb)));
	isoval = 100.0f;
	for(int i=0; i<numb; i++)
	{
		sdfValues[i] = SDFLibrary::values[i].value * SDFLibrary::values[i].signe;
	}
	free_memory();
	return (sdfValues);
}

RAWIV_header* SDFLibrary::getVolumeInfo()
{
	RAWIV_header* ret = (RAWIV_header*)(malloc(sizeof(RAWIV_header)*1));
	for(int i=0; i<3; i++)
	{
		ret->minext[i] = SDFLibrary::minext[i];
		ret->maxext[i] = SDFLibrary::maxext[i];
		ret->span[i] = SDFLibrary::span[i];
		ret->origin[i] = 0.0f;
		ret->dim[i] = SDFLibrary::size+1;
	}
	ret->ngridpts = (SDFLibrary::size+1)*(SDFLibrary::size+1)*(SDFLibrary::size+1);
	ret->ncells = (SDFLibrary::size)*(SDFLibrary::size)*(SDFLibrary::size);
	ret->size = SDFLibrary::size;
	return ret;
}

SimpleVolumeData* SDFLibrary::getSDF(Geometry* geometry, int size)
{
	if(!geometry)
	{
		return 0;
	}
	float mins[3];
	float maxs[3];
	geometry->CalculateExtents();
	mins[0] = geometry->m_Min[0];
	mins[1] = geometry->m_Min[1];
	mins[2] = geometry->m_Min[2];
	maxs[0] = geometry->m_Max[0];
	maxs[1] = geometry->m_Max[1];
	maxs[2] = geometry->m_Max[2];
	double w0 = maxs[0] - mins[0];
	double w1 = maxs[1] - mins[1];
	double w2 = maxs[2] - mins[2];
	double tenPercent0 = w0/10;
	double tenPercent1 = w1/10;
	double tenPercent2 = w2/10;
	mins[0] -= tenPercent0;
	mins[1] -= tenPercent1;
	mins[2] -= tenPercent2;
	maxs[0] += tenPercent0;
	maxs[1] += tenPercent1;
	maxs[2] += tenPercent2;
	int flipNormals = 0;
	if(!SDFLibrary::setParameters(size, flipNormals, mins, maxs))
	{
		return 0;
	}
	//Finally, call the SDF code.
	float* values = SDFLibrary::computeSDF(geometry->m_NumTriVerts, geometry->m_TriVerts, geometry->m_NumTris, (int*)(geometry->m_Tris));
	if(!values)
	{
		return 0;
	}
	float* choppedValues = new float[size*size*size];
	{
		int c=0;
		for(int i=0; i<=size; i++)
		{
			for(int j=0; j<=size; j++)
			{
				for(int k=0; k<=size; k++)
				{
					if(i!=size && j!=size && k!=size)
					{
						choppedValues[c++] = values[i*(size+1)*(size+1) + j*(size+1) + k];
					}
				}
			}
		}
	}
	delete[] values;
	//Lastly, the header info.
	SDFLibrary::RAWIV_header* volInfo = SDFLibrary::getVolumeInfo();
	if(!volInfo)
	{
		return false;
	}
	SimpleVolumeData* sData = 0;
	{
		unsigned int dims[3];
		dims[0] = size;
		dims[1] = size;
		dims[2] = size;
		sData = new SimpleVolumeData(dims); // Johns oddities
		sData->setDimensions(dims);
	}
	sData->setNumberOfVariables(1);
	sData->setData(0, choppedValues);
	sData->setType(0, SimpleVolumeData::FLOAT);
	sData->setName(0, "Sign distance function from TexMol");
	sData->setMinExtent(volInfo->minext);
	sData->setMaxExtent(volInfo->maxext);
	return sData;
}
