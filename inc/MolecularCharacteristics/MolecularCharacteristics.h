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
#ifndef CCV_MOLECULAR_CHARACTERISTICS_H
#define CCV_MOLECULAR_CHARACTERISTICS_H

#include <Utility/utility.h>
#include <Blurmaps/AreaVolume.h>
#include <Blurmaps/BlurMapsDataManager.h>
#include <Geometry/Geometry.h>
#include <PDBParser/Atom.h>
#include <PDBParser/ElementInformation.h>
#include <PDBParser/GroupOfAtoms.h>
#include <SimpleVolumeData/SimpleVolumeData.h>
#include <SimpleVolumeData/SimpleVolumeDataIsocontourer.h>
#include <VolumeFileTypes/VolumeLoader.h>

namespace PDBParser
{
	class GroupOfAtoms;
};
class SimpleVolumeData;

class MolecularCharacteristics
{
	public:
		MolecularCharacteristics(PDBParser::GroupOfAtoms* molecule, SimpleVolumeData* electrostatics, int numberOfSpheres, float* activeSite, double* radii, char* inputGOAFileName, char* inputElecVolumeFileName);
		virtual ~MolecularCharacteristics();
		bool getMoleculeCharacteristics(int* numberOfAtoms, int* numberOfHydrophobicAtoms, int* numberOfHydrophilicAtoms, double* numberOfResidues, double* numberOfHydrophobicResidues, double* numberOfHydrophilicResidues, double* area, double* volume, int* numberOfAtomsInSphere, int* numberOfHydrophobicAtomsInSphere, int* numberOfHydrophilicAtomsInSphere, double* numberOfResiduesInSphere, double* numberOfHydrophobicResiduesInSphere, double* numberOfHydrophilicResiduesInSphere, double* areaInSphere, double* volumeInSphere);
		bool appendCharacteristics(const char* filename);
		bool getAreaVolume(double* area, double* volume, double* areaInSphere, double* volumeInSphere);
		bool getElectrostaticCharacteristics(double* ratio, double*  sphereRatios);
		bool appendHeader(const char* filename);

	protected:
		double dist2(float* p1, float* p2);
		bool getNumberOfAtoms(int* numberOfAtoms);
		bool getNumberOfAtoms(int* numberOfAtoms, float* point, double radius, PDBParser::GroupOfAtoms* molecule);
		bool getNumberOfChargedAtoms(int* numberOfHydrophobicAtoms, int* numberOfHydrophilicAtoms, PDBParser::GroupOfAtoms* molecule);
		bool getNumberOfChargedAtoms(int* numberOfHydrophobicAtoms, int* numberOfHydrophilicAtoms, float* point, double radius, PDBParser::GroupOfAtoms* molecule);
		bool zeroOutVolumeShell(SimpleVolumeData* sData, float* center, double radius, char* filename) ;
		bool getElecRatio(float* data, int* dim, double* ratio);
		PDBParser::GroupOfAtoms* m_Molecule;
		SimpleVolumeData* m_Electrostatics;
		int m_NumberOfSpheres;
		float* m_ActiveSite;
		double* m_Radii;
		char m_GOAFileName[1024];
		char m_ElecFileName[1024];
};

#endif
