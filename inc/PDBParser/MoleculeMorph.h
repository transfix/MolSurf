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
#ifndef CCV_MOLECULE_MORPH_H
#define CCV_MOLECULE_MORPH_H

#include <Utility/utility.h>
#include <PDBParser/GroupOfAtoms.h>
#include <PDBParser/ResidueData.h>

namespace PDBParser
{
	class GroupOfAtoms;
	class MoleculeMorph
	{
		public:
			MoleculeMorph();
			virtual ~MoleculeMorph();
			bool morph(PDBParser::GroupOfAtoms* molecule1, PDBParser::GroupOfAtoms* molecule2, const char* outputFileNamePrefix, double resolution, int maxSteps);
		protected:
			double getSize();
			bool torsionAnglesPresent(PDBParser::GroupOfAtoms* molecule);
			bool getTorsionAngles(PDBParser::GroupOfAtoms* molecule, double* torsionAngles, int* count);
			bool getResidueSequenceIDs(PDBParser::GroupOfAtoms* molecule, int* residueSequenceNumbers, int* count);
			PDBParser::GroupOfAtoms* m_Molecule1;
			PDBParser::GroupOfAtoms* m_Molecule2;
			const char* m_OutputFileNamePrefix;
			double m_Resolution;
	};
};

#endif
