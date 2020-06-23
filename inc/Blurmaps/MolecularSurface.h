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
#if !defined(AFX_MOLECULARSURFACE_H__EEA3AF45_ABF1_46D3_AA6B_121E2FDF6592__INCLUDED_)
#define AFX_MOLECULARSURFACE_H__EEA3AF45_ABF1_46D3_AA6B_121E2FDF6592__INCLUDED_

#include <Utility/utility.h>
#include <Blurmaps/AtomPopulator.h>
#include <Blurmaps/BlurMapsDataManager.h>
#include <Blurmaps/SkinRegion2.h>
#include <Geometry/Geometry.h>
#include <PDBParser/CollectionData.h>
#include <PDBParser/GroupOfAtoms.h>
#include <SimpleVolumeData/SimpleVolumeData.h>
#include <SimpleVolumeData/SimpleVolumeDataIsocontourer.h>
#include <VolumeFileTypes/VolumeLoader.h>

namespace PDBParser
{
	class GroupOfAtoms;
};

class Geometry;
class SimpleVolumeData;

class MolecularSurface
{
	public:
		MolecularSurface();
		virtual ~MolecularSurface();

		static Geometry* getMolecularSurface(PDBParser::GroupOfAtoms* molecule, int size);
		static Geometry* getMolecularSurface(PDBParser::GroupOfAtoms* molecule, int size,double enlargedRadius);
};

#endif
