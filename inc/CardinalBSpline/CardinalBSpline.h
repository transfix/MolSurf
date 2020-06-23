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
#ifndef CARDINALBSPLINE_H
#define CARDINALBSPLINE_H

#include <Utility/utility.h>
#include <Blurmaps/BlurMapsDataManager.h>
#include <CardinalBSpline/Bbasis.h>
#include <PDBParser/Atom.h>
#include <SimpleVolumeData/SimpleVolumeData.h>

using namespace PDBParser;

class SimpleVolumeData;

namespace PDBParser
{
	class GroupOfAtoms;
	class Atom;
};

class CardinalBSpline
{
	public:
		CardinalBSpline();
		virtual ~CardinalBSpline();
		SimpleVolumeData* getCardinalBSplineSurface(PDBParser::GroupOfAtoms* molecule, unsigned int* dim, int order);
		bool getAtomListAndExtent(PDBParser::GroupOfAtoms* molecule, vector<PDBParser::Atom*> &atomList, float* minExtent, float* maxExtent);
		bool computeFunction(std::vector<PDBParser::Atom*> atomList, float* data, unsigned int* dim,float* minExtent,float* maxExtent,int order);
		void CardinalBSpline_One_Atom(float* position,float radius,int order, unsigned int* dim, float* minExtent, float* maxExtent, float* data);
	private:
		float minExtent[3],maxExtent[3];
};

#endif
