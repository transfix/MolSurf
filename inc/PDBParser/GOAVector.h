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
#ifndef __GOAVECTOR_H__
#define __GOAVECTOR_H__

#include <Utility/utility.h>
#include <PDBParser/Atom.h>

namespace PDBParser
{
	class Atom;
	// Implements 3-D vector ops
	class GOAVector
	{
		public:
			GOAVector();
			GOAVector(double ax, double ay, double az);
			GOAVector(Atom* atom);
			~GOAVector();

			void print();
			void norm();

			double x, y, z;
	};
	GOAVector* subtractAtoms(Atom* atomA, Atom* atomB);
	GOAVector* scalarProduct(double c, GOAVector* v1);
	GOAVector* crossProduct(GOAVector* v1, GOAVector* v2);
	double dotProduct(GOAVector* v1, GOAVector* v2);
	double dotProductAngle(GOAVector* v1, GOAVector* v2);
};

#endif
