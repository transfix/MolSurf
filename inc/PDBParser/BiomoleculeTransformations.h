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
#if !defined(AFX_BIOMOLECULETRANSFORMATIONS_H__2C4DCE88_A8F6_4E08_BEC6_718DB814A308__INCLUDED_)
#define AFX_BIOMOLECULETRANSFORMATIONS_H__2C4DCE88_A8F6_4E08_BEC6_718DB814A308__INCLUDED_

#include <Utility/utility.h>

namespace PDBParser
{
	// We store a set of chain ids, and the transformations that needs to be applied to them.
	class BiomoleculeTransformations
	{
		public:
			BiomoleculeTransformations();
			BiomoleculeTransformations(const BiomoleculeTransformations& b);
			virtual ~BiomoleculeTransformations();
			bool print(FILE* fp);
			bool isTransformationPresent(char chainId);
			int getNumberOfTransformations(char chainId);
			bool getRotations(char chainId, double* rotations);
			bool getTranslations(char chainId, double* translation);
			int getNumberOfTransformations();
			bool getRotations(double* rotations);
			bool getTranslations(double* translation);
			void addChains(std::vector<char> chainIds);
			std::vector<char> chainIds;
			std::vector<double> rotation;
			std::vector<double> translation;
	};
};

#endif
