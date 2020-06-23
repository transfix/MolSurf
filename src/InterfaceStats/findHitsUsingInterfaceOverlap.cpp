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
#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<vector>
#include<string>

int main(int argc, char** argv)
{
	if(argc<9)
	{
		printf("Usage: findHitsUsingInterfaceOverlap <receptor.rawn> <ligand1.rawn> <ligand2.rawn> <overlapRatioThreshold> <interface width> <xformFile> <startIndex> <endIndex>");
		return -1;
	}

	FILE *xformFile = fopen(argv[6], "rt");
	int startIndex = atoi(argv[7]);
	int endIndex = atoi(argv[8]);
	int genSurfMode = atoi(argv[4]);
	double interfaceWidth = atof(argv[5]);
	double overlapRatioThreshold = atof(argv[4]);

	printf("interfaceWidth = %lf\n", interfaceWidth);
	printf("overlapRatioThreshold = %lf\n", overlapRatioThreshold);

	double areaOfReceptorsInterfaceWithStaticLigand = 0.0;
	double areaOfMovingLigandsInterfaceWithReceptor = 0.0;
	double areaOfinterfaceOverlap = 0.0;
	double ratio, ratio2;



	InterfaceFootprintOverlap* ifo = new InterfaceFootprintOverlap(std::string(argv[1]), std::string(argv[2]), std::string(argv[3]), interfaceWidth, false);
	areaOfReceptorsInterfaceWithStaticLigand = ifo->getAreaOfReceptorInterfaceWithStaticLigand();

	printf("Footprint of static ligand on receptor = %lf\n", areaOfReceptorsInterfaceWithStaticLigand);

	if(areaOfReceptorsInterfaceWithStaticLigand == 0.0)
	{
		printf("Static Ligand does not have interface with receptor. Cannot compute Overlap. Quitting");
		return -2;
	}


       	Matrix mtx;

	if(xformFile==NULL)
	{
		printf("Could not open xform file\n");
		return -3;
	}

	int numXForm;
	fscanf(xformFile, "%d", &numXForm);

	if(startIndex>=numXForm || startIndex<0 || endIndex>=numXForm || endIndex<0 || endIndex<startIndex)
	{
		printf("Xform Index out of bounds\n");
		return -4;
	}

	double dummy;

	int hitsBelow1 = 0;
	int hitsBelow10 = 0;
	int hitsBelow100 = 0;
	int hitsBelow1000 = 0;

	int hitsAbove10PercentOverlap = 0;
	int hitsAbove25PercentOverlap = 0;
	int hitsAbove50PercentOverlap = 0;
	int hitsAbove60PercentOverlap = 0;
	int hitsAbove70PercentOverlap = 0;
	int hitsAbove80PercentOverlap = 0;
	int hitsAbove90PercentOverlap = 0;

	int hits2Below1 = 0;
	int hits2Below10 = 0;
	int hits2Below100 = 0;
	int hits2Below1000 = 0;

	int hits2Above10PercentOverlap = 0;
	int hits2Above25PercentOverlap = 0;
	int hits2Above50PercentOverlap = 0;
	int hits2Above60PercentOverlap = 0;
	int hits2Above70PercentOverlap = 0;
	int hits2Above80PercentOverlap = 0;
	int hits2Above90PercentOverlap = 0;

	printf("\n\nRank \t\t OverlapRatio1 \t\t OverlapRatio2 \t\t OverlapArea\n");
	printf(    "---- \t\t ------------- \t\t ------------- \t\t -----------\n");
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

		fscanf(xformFile, "%lf", &dummy);
		fscanf(xformFile, "%lf", &dummy);

		if(i>=startIndex)
		{	 
			if(ifo->computeInterfaces(mtx))
			{
				areaOfinterfaceOverlap = ifo->getAreaOfOverlap();
				areaOfMovingLigandsInterfaceWithReceptor = ifo->getAreaOfMovingLigandInterface();
				ratio = areaOfinterfaceOverlap/areaOfReceptorsInterfaceWithStaticLigand;
				ratio2 = areaOfinterfaceOverlap/areaOfMovingLigandsInterfaceWithReceptor;

				if(ratio>0.9) hitsAbove90PercentOverlap++;
				else if(ratio>0.8) hitsAbove80PercentOverlap++;
				else if(ratio>0.7) hitsAbove70PercentOverlap++;
				else if(ratio>0.6) hitsAbove60PercentOverlap++;
				else if(ratio>0.5) hitsAbove50PercentOverlap++;
				else if(ratio>0.25) hitsAbove25PercentOverlap++;
				else if(ratio>0.1) hitsAbove10PercentOverlap++;

				if(ratio2>0.9) hits2Above90PercentOverlap++;
				else if(ratio2>0.8) hits2Above80PercentOverlap++;
				else if(ratio2>0.7) hits2Above70PercentOverlap++;
				else if(ratio2>0.6) hits2Above60PercentOverlap++;
				else if(ratio2>0.5) hits2Above50PercentOverlap++;
				else if(ratio2>0.25) hits2Above25PercentOverlap++;
				else if(ratio2>0.1) hits2Above10PercentOverlap++;

				if(ratio > overlapRatioThreshold || ratio2 > overlapRatioThreshold )
				{
					printf("%d \t\t %lf \t\t %lf \t\t %lf\n", i+1, ratio, ratio2, areaOfinterfaceOverlap);

					if(ratio > overlapRatioThreshold)
					{
						if(i<1) hitsBelow1 ++;
						if(i<10) hitsBelow10 ++;
						if(i<100) hitsBelow100 ++;
						if(i<10000) hitsBelow1000 ++;
					}
					if(ratio2 > overlapRatioThreshold)
					{				
						if(i<1) hits2Below1 ++;
						if(i<10) hits2Below10 ++;
						if(i<100) hits2Below100 ++;
						if(i<10000) hits2Below1000 ++;
					}
				}
			}
			else
			{
				printf("ERROR.\n ERROR.ERROR.\n ERROR.ERROR.ERROR.\n");
			}
		}
	}

	fclose(xformFile);

	printf("\n\n\n");
	printf("Ratio computed as overlap/(interface with static mol)\n\n");

	printf("Hits in [1,1]: %d\n", hitsBelow1);
	printf("Hits in [1,10]: %d\n", hitsBelow10);
	printf("Hits in [1,100]: %d\n", hitsBelow100);
	printf("Hits in [1,1000]: %d\n", hitsBelow1000);

	printf("\n\n");
	printf("Hits with more than 90 percent footprint overlap: %d\n", hitsAbove90PercentOverlap);
	printf("Hits with more than 80 percent footprint overlap: %d\n", hitsAbove80PercentOverlap);
	printf("Hits with more than 70 percent footprint overlap: %d\n", hitsAbove70PercentOverlap);
	printf("Hits with more than 60 percent footprint overlap: %d\n", hitsAbove60PercentOverlap);
	printf("Hits with more than 50 percent footprint overlap: %d\n", hitsAbove50PercentOverlap);
	printf("Hits with more than 25 percent footprint overlap: %d\n", hitsAbove25PercentOverlap);
	printf("Hits with more than 10 percent footprint overlap: %d\n", hitsAbove10PercentOverlap);


	printf("\n\n\n");
	printf("Ratio computed as overlap/(interface with moving mol)\n\n");

	printf("Hits in [1,1]: %d\n", hits2Below1);
	printf("Hits in [1,10]: %d\n", hits2Below10);
	printf("Hits in [1,100]: %d\n", hits2Below100);
	printf("Hits in [1,1000]: %d\n", hits2Below1000);

	printf("\n\n");
	printf("Hits with more than 90 percent footprint overlap: %d\n", hits2Above90PercentOverlap);
	printf("Hits with more than 80 percent footprint overlap: %d\n", hits2Above80PercentOverlap);
	printf("Hits with more than 70 percent footprint overlap: %d\n", hits2Above70PercentOverlap);
	printf("Hits with more than 60 percent footprint overlap: %d\n", hits2Above60PercentOverlap);
	printf("Hits with more than 50 percent footprint overlap: %d\n", hits2Above50PercentOverlap);
	printf("Hits with more than 25 percent footprint overlap: %d\n", hits2Above25PercentOverlap);
	printf("Hits with more than 10 percent footprint overlap: %d\n", hits2Above10PercentOverlap);
}
