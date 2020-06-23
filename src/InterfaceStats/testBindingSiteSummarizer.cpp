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

#include<InterfaceStats/BindingSiteSummarizer.h>
#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<vector>
#include<string>

int main(int argc, char** argv)
{
	if(argc<8)
	{
		printf("Usage: testBindingSiteSummarizer <receptor.rawn> <ligand.rawn> <interface width> <xformFile> <startIndex> <endIndex> <output rawn>\n");
		return -1;
	}

	double interfaceWidth = atof(argv[3]);
	FILE *xformFile = fopen(argv[4], "rt");
	int startIndex = atoi(argv[5]);
	int endIndex = atoi(argv[6]);

	BindingSiteSummarizer* bss = new BindingSiteSummarizer(std::string(argv[1]), std::string(argv[2]), std::string(argv[7]), interfaceWidth, true);

       	std::vector<Matrix> xforms;
       	std::vector<double> scores;

	if(xformFile==NULL)
	{
		printf("Could not open xform file\n");
		return -2;
	}

	int numXForm;
	fscanf(xformFile, "%d", &numXForm);

	if(startIndex>=numXForm || startIndex<0 || endIndex>=numXForm || endIndex<0 || endIndex<startIndex)
	{
		printf("Xform Index out of bounds\n");
		return -3;
	}

	double rmsd;
	double score;
	Matrix mtx;
	for(int i=0; i<=endIndex; i++)
	{
		for(int j=0; j<4; j++)
		{
			for(int k=0; k<4; k++)
			{
				double mtx_jk;
				fscanf(xformFile, "%lf", &mtx_jk);
				//printf("%lf ", mtx_jk);
				mtx.set(j,k,mtx_jk);
			}
		}

		fscanf(xformFile, "%lf", &score);
		fscanf(xformFile, "%lf", &rmsd);

		if(i>=startIndex)
		{	 
			xforms.push_back(mtx);
			scores.push_back(score);
		}
	}

//	bss->computeSiteSummary(xforms, scores);
	bss->computeSiteSummary(xforms);

	fclose(xformFile);

	


}
