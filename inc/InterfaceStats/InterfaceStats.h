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
#ifndef INT_STATS_H
#define INT_STATS_H

#include <cstdio>

using namespace std;

class InterfaceStats
{
	private:
		double interfaceWidth;

		double planarityIndex;
		double circularityIndex;
	
		double interfaceArea;
		double interfaceVolume;
		double gapIndex;

		double localDensityIndex;
		double globalDensityIndex;

		int resresContact[20][20];

		int numCoreResidue[2];
		int numRimResidue[2];

		int numPolarResidue[2];
		int numNonPolarResidue[2];

		int numAtoms[2];


	public:
		InterfaceStats(){};
		~InterfaceStats(){};

		void setInterfaceWidth(double d) {interfaceWidth = d;}
		void setPlanarityIndex(double d) {planarityIndex = d;}
		void setCircularityIndex(double d) {circularityIndex = d;}
		void setInterfaceArea(double d) {interfaceArea = d;}
		void setInterfaceVolume(double d) {interfaceVolume = d;}
		void setGapIndex(double d) {gapIndex = d;}
		void setLocalDensityIndex(double d) {localDensityIndex = d;}
		void setGlobalDensityIndex(double d) {globalDensityIndex = d;}
		void setResresContact(int i, int j, int a) {resresContact[i][j] = a;}
		void setNumCoreResidue(int i, int stc) {numCoreResidue[stc] = i;}
		void setNumRimResidue(int i, int stc) {numRimResidue[stc] = i;}
		void setNumPolarResidue(int i, int stc) {numPolarResidue[stc] = i;}
		void setNumNonPolarResidue(int i, int stc) {numNonPolarResidue[stc] = i;}
		void setNumAtoms(int i, int stc) {numAtoms[stc] = i;}

		double getInterfaceWidth() {return interfaceWidth;}
		double getPlanarityIndex() {return planarityIndex;}
		double getCircularityIndex() {return circularityIndex;}
		double getInterfaceArea() {return interfaceArea;}
		double getInterfaceVolume() {return interfaceVolume;}
		double getGapIndex() {return gapIndex;}
		double getLocalDensityIndex() {return localDensityIndex;}
		double getGlobalDensityIndex() {return globalDensityIndex;}
		int getResresContact(int i, int j) {return resresContact[i][j];}

		int getNumInterfaceAtoms() {return numAtoms[0] + numAtoms[1];}
		int getNumInterfaceAtoms(int stc) {return numAtoms[stc];}

		int getNumInterfaceResidues() {return getNumPolarResidue() + getNumNonPolarResidue();}
		int getNumInterfaceResidues(int stc) {return numPolarResidue[stc] + numNonPolarResidue[stc];}

		int getNumCoreResidue() {return numCoreResidue[0] + numCoreResidue[1];}
		int getNumCoreResidue(int stc) {return numCoreResidue[stc];}

		int getNumRimResidue() {return numRimResidue[0] + numRimResidue[1];}
		int getNumRimResidue(int stc) {return numRimResidue[stc];}

		int getNumPolarResidue() {return numPolarResidue[0] + numPolarResidue[1];}
		int getNumPolarResidue(int stc) {return numPolarResidue[stc];}

		int getNumNonPolarResidue() {return numNonPolarResidue[0] + numPolarResidue[1];}
		int getNumNonPolarResidue(int stc) {return numNonPolarResidue[stc];}

		void print()
		{
			printf("interfaceWidth: %lf \n", interfaceWidth);
			printf("planarityIndex: %lf \n", planarityIndex);
			printf("circularityIndex: %lf \n", circularityIndex);

			printf("interfaceArea: %lf \n", interfaceArea);
			printf("interfaceVolume: %lf \n", interfaceVolume);
			printf("gapIndex: %lf \n", gapIndex);

			printf("localDensityIndex: %lf \n", localDensityIndex);
			printf("globalDensityIndex: %lf \n", globalDensityIndex);

			printf("numCoreResidue 1: %d \n", numCoreResidue[0]);
			printf("numCoreResidue 2: %d \n", numCoreResidue[1]);

			printf("numRimResidue 1: %d \n", numRimResidue[0]);
			printf("numRimResidue 2: %d \n", numRimResidue[1]);

			printf("numPolarResidue 1: %d \n", numPolarResidue[0]);
			printf("numPolarResidue 2: %d \n", numPolarResidue[1]);

			printf("numNonPolarResidue 1: %d \n", numNonPolarResidue[0]);
			printf("numNonPolarResidue 2: %d \n", numNonPolarResidue[1]);

			printf("numAtoms 1: %d \n", numAtoms[0]);
			printf("numAtoms 2: %d \n", numAtoms[1]);
		}
};


#endif
