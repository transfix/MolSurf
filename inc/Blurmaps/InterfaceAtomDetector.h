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
#ifndef CCV_INTERFACE_ATOM_DETECTOR_H
#define CCV_INTERFACE_ATOM_DETECTOR_H

#include<string>
#include <Utility/utility.h>
#include <Blurmaps/BlurMapsDataManager.h>
#include <Blurmaps/GOABlur.h>
#include <PDBParser/Atom.h>
#include <PDBParser/CollectionData.h>
#include <PDBParser/ElementInformation.h>
#include <PDBParser/GroupOfAtoms.h>
#include <PDBParser/ResidueInformation.h>
#include <SimpleVolumeData/SimpleVolumeData.h>
#include <UsefulMath/TrilinearGrid.h>

using namespace std;

namespace PDBParser
{
	class GroupOfAtoms;
	class Atom;
};

class InterfaceAtomDetector
{
	public:
		InterfaceAtomDetector(PDBParser::GroupOfAtoms* molecule1, PDBParser::GroupOfAtoms* molecule2, string outputFileName, double interfaceWidth, string atomTypes, PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType);

		virtual ~InterfaceAtomDetector();

		bool detect();
		vector<int> getInterfaceAtoms(){return selectedIds;}

	private:
		vector<PDBParser::Atom*> m_AtomList1;
		vector<PDBParser::Atom*> m_AtomList2;
		vector <int> selectedIds;

		string outputFileName;
		string atomTypes;
		double interfaceWidth;
};


#endif
