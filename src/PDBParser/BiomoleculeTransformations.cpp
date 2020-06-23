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
#include <PDBParser/BiomoleculeTransformations.h>

using PDBParser::BiomoleculeTransformations;

BiomoleculeTransformations::BiomoleculeTransformations()
{
}

BiomoleculeTransformations::BiomoleculeTransformations(const BiomoleculeTransformations& b)
{
	chainIds = b.chainIds;
	rotation = b.rotation;
	translation = b.translation;
}

BiomoleculeTransformations::~BiomoleculeTransformations()
{
}

bool BiomoleculeTransformations::print(FILE* fp)
{
	if(!fp)
	{
		return false;
	}
	if(chainIds.size() == 0 || rotation.size() == 0 || translation.size() == 0)
	{
		return false;
	}
	// This is my error checking for transformations.
	int numRotations = rotation.size();
	int numTranslations = translation.size();
	if((numTranslations % 3) != 0 || (numRotations % 9) != 0)
	{
		return false;
	}
	if(numRotations != numTranslations*3)
	{
		return false;
	}
	int i;
	fprintf(fp, "\tChain ids to transform: ");
	for(i=0; i<chainIds.size()-1; i++)
	{
		fprintf(fp, "%c, ", chainIds.at(i));
	}
	fprintf(fp, "%c\n", chainIds.at(chainIds.size()-1));
	for(i=0; i<numRotations/9; i++)
	{
		fprintf(fp, "\tRotation matrix:\n");
		fprintf(fp, "\t[%f %f %f]\n\t[%f %f %f]\n\t[%f %f %f]\n",
				rotation.at(i*9+0), rotation.at(i*9+1), rotation.at(i*9+2),
				rotation.at(i*9+3), rotation.at(i*9+4), rotation.at(i*9+5),
				rotation.at(i*9+6), rotation.at(i*9+7), rotation.at(i*9+8));
		fprintf(fp, "\tTranslation matrix:\n");
		fprintf(fp, "\t[%f %f %f]\n", translation.at(i*3+0), translation.at(i*3+1), translation.at(i*3+2));
	}
	return true;
}

void BiomoleculeTransformations::addChains(std::vector<char> newchainIds)
{
	int n = newchainIds.size();
	int i;
	for(i=0; i<n; i++)
	{
		chainIds.push_back(newchainIds.at(i));
	}
}

bool BiomoleculeTransformations::isTransformationPresent(char currentChainId)
{
	int n = chainIds.size();
	int i;
	for(i=0; i<n; i++)
	{
		if(chainIds.at(i) == currentChainId)
		{
			return true;
		}
	}
	return false;
}

int BiomoleculeTransformations::getNumberOfTransformations(char currentChainId)
{
	int n = chainIds.size();
	int i;
	for(i=0; i<n; i++)
	{
		if(chainIds.at(i) == currentChainId)
		{
			int s = rotation.size();
			return rotation.size()/9; // array should be a multiple of 9.
		}
	}
	return -1; // error, did not find the chainId in this transformation list
}

int BiomoleculeTransformations::getNumberOfTransformations()
{
	int s = rotation.size();
	return rotation.size()/9; // array should be a multiple of 9.
}

bool BiomoleculeTransformations::getRotations(char currentChainId, double* rotationsArray)
{
	if(!rotationsArray)
	{
		return false;
	}
	int n = chainIds.size();
	int i;
	for(i=0; i<n; i++)
	{
		if(chainIds.at(i) == currentChainId)
		{
			int j;
			for(j=0; j<rotation.size(); j++)
			{
				rotationsArray[j] = rotation.at(j);
			}
			return true;
		}
	}
	return false; // error, did not find the chainId in this transformation list
}

bool BiomoleculeTransformations::getRotations(double* rotationsArray)
{
	if(!rotationsArray)
	{
		return false;
	}
	int j;
	for(j=0; j<rotation.size(); j++)
	{
		rotationsArray[j] = rotation.at(j);
	}
	return true;
}

bool BiomoleculeTransformations::getTranslations(char currentChainId, double* translationArray)
{
	if(!translationArray)
	{
		return false;
	}
	int n = chainIds.size();
	int i;
	for(i=0; i<n; i++)
	{
		if(chainIds.at(i) == currentChainId)
		{
			int j;
			for(j=0; j<translation.size(); j++)
			{
				translationArray[j] = translation.at(j);
			}
			return true;
		}
	}
	return false; // error, did not find the chainId in this transformation list
}

bool BiomoleculeTransformations::getTranslations(double* translationArray)
{
	if(!translationArray)
	{
		return false;
	}
	int j;
	for(j=0; j<translation.size(); j++)
	{
		translationArray[j] = translation.at(j);
	}
	return true;
}
