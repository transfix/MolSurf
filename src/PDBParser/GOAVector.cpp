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
// GOAVector.cpp: Implements 3-D vector ops

#include <PDBParser/GOAVector.h>

using PDBParser::GOAVector;
using PDBParser::Atom;

GOAVector::GOAVector()
{
	x = 0;
	y = 0;
	z = 0;
}

GOAVector::GOAVector(double ax, double ay, double az)
{
	x = ax;
	y = ay;
	z = az;
}

GOAVector::GOAVector(Atom* atom)
{
	x = atom->m_Position[0];
	y = atom->m_Position[1];
	z = atom->m_Position[2];
}

GOAVector::~GOAVector()
{
}

void GOAVector::print(void)
{
	printf("GOAVector: ( %lf, %lf, %lf)\n", x, y, z);
}

// Destructively normalizes a non-0 vector.  Does nothing on 0 vectors.
void GOAVector::norm(void)
{
	double len = sqrt((x * x) + (y * y) + (z * z));
	if(len != 0)
	{
		x /= len;
		y /= len;
		z /= len;
	}
}

// Returns a vector equal to the center of atom_a minus the center of atom_b
GOAVector* PDBParser::subtractAtoms(Atom* atomA, Atom* atomB)
{
	if(atomA == 0 || atomB == 0)
	{
		return 0;
	}
	return new GOAVector(atomA->m_Position[0] - atomB->m_Position[0],
						 atomA->m_Position[1] - atomB->m_Position[1],
						 atomA->m_Position[2] - atomB->m_Position[2]);
}

GOAVector* PDBParser::scalarProduct(double c, GOAVector* v1)
{
	return new GOAVector(c * v1->x, c * v1->y, c * v1->z);
}

GOAVector* PDBParser::crossProduct(GOAVector* v1, GOAVector* v2)
{
	GOAVector* v3 = new GOAVector();
	v3->x = (v1->y * v2->z) - (v1->z * v2->y);
	v3->y = (v1->z * v2->x) - (v1->x * v2->z);
	v3->z = (v1->x * v2->y) - (v1->y * v2->x);
	return v3;
}

double PDBParser::dotProduct(GOAVector* v1, GOAVector* v2)
{
	return ((v1->x * v2->x) + (v1->y * v2->y) + (v1->z * v2->z));
}

// Returns the inverse cosin of the dot product, which is the angle between v1 and
// v2 if they are unit.
double PDBParser::dotProductAngle(GOAVector* v1, GOAVector* v2)
{
	return (acos((v1->x * v2->x) + (v1->y * v2->y) + (v1->z * v2->z)));
}
