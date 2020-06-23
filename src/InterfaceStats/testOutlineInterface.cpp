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

#include<InterfaceStats/OutlineInterface.h>
#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<vector>
#include<string>

int main(int argc, char** argv)
{
	if(argc<11)
	{
		printf("Usage: testOutlineInterface <receptor.rawc> <ligand.rawn/rawc> <interface width> <xformFile> <xformIndex> <output rawc> <rawc/rawnc [0/1]> <outlineColorR> <outlineColorG> <outlineColorB>\n");
		return -1;
	}

	double interfaceWidth = atof(argv[3]);
	FILE *xformFile = fopen(argv[4], "rt");
	int xformIndex = atoi(argv[5]);
    
    int rawcMode = atoi(argv[7]);

	OutlineInterface* oi = new OutlineInterface(std::string(argv[1]), std::string(argv[2]), interfaceWidth, rawcMode==0?false:true, true);

	if(xformFile==NULL)
	{
		printf("Could not open xform file\n");
		return -2;
	}

	int numXForm;
	fscanf(xformFile, "%d", &numXForm);

	if(xformIndex>=numXForm || xformIndex<0)
	{
		printf("Xform Index out of bounds\n");
		return -3;
	}

	Matrix mtx;
	double score, rmsd;
	for(int i=0; i<=xformIndex; i++)
	{
		for(int j=0; j<4; j++)
		{
			for(int k=0; k<4; k++)
			{
				double mtx_jk;
				fscanf(xformFile, "%lf", &mtx_jk);
				mtx.set(j,k,mtx_jk);
			}
		}

		fscanf(xformFile, "%lf", &score);
		fscanf(xformFile, "%lf", &rmsd);
	}

    vertexcolor outlineColor;
    outlineColor.r = atof(argv[8]);
    outlineColor.g = atof(argv[9]);
    outlineColor.b = atof(argv[10]);

	oi->computeOutline(mtx, std::string(argv[6]), outlineColor);

	fclose(xformFile);

	return 0;
}
