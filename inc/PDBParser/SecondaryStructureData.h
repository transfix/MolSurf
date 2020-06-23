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
#if !defined(AFX_SECONDARYSTRUCTUREDATA_H__CFE99172_63D7_4387_8F80_1E686EB44A3D__INCLUDED_)
#define AFX_SECONDARYSTRUCTUREDATA_H__CFE99172_63D7_4387_8F80_1E686EB44A3D__INCLUDED_

#include <Utility/utility.h>
#include <UsefulMath/Vector.h>
#include <PDBParser/ResidueHeader.h>
#include <PDBParser/ResidueInformation.h>
#include <PDBParser/ElementInformation.h>
#include <PDBParser/ResidueData.h>
#include <PDBParser/SecondaryStructureData.h>
#include <PDBParser/ProteinData.h>
#include <PDBParser/CollectionData.h>
#include <PDBParser/BiomoleculeTransformations.h>

namespace PDBParser
{
	class SecondaryStructureData
	{
		public:
			SecondaryStructureData();
			// Copies all substructures recursively
			SecondaryStructureData(const SecondaryStructureData& secondaryStructureData);	
			virtual ~SecondaryStructureData();
			// Used only for sheets
			// We store a pair of ints for the residue sequence numbers.
			// Each pair defines a strand in a sheet
			std::vector<int> m_ResSequenceNumberPairList;
			// The geometry is defined as a 7 tuple, containing
			std::vector<double> m_StrandGeometry;
			// Used for Helices
			std::vector<double> m_HelixGeometry;
	};
};

#endif
