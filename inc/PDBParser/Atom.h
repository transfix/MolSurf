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
#ifndef CCV_TEXMOL_PDBPARSER_ATOM_H
#define CCV_TEXMOL_PDBPARSER_ATOM_H

#include <Utility/utility.h>
#include <PDBParser/GroupOfAtoms.h>
#include <UsefulMath/Matrix.h>
#include <UsefulMath/Ray.h>
#include <DPG/PG.h>
#include <PDBParser/ResidueInformation.h>
#include <PDBParser/ElementInformation.h>
#include <PDBParser/MutatedAtomInformation.h>
#include <UsefulMath/Vector.h>


// arand: don't using the namespace everywhere since
//        DPG::Point conflicts with some CGAL stuff...
//using namespace DPG;

namespace PDBParser
{
	class MutatedAtomInformation;
	class Atom: public DPG::Ball
	{
		protected:

		public:
			Atom();
			Atom(const Atom& a);
			Atom& operator=(const Atom& a);
			virtual ~Atom();
			float getRadius() const;
			float getRadius(GroupOfAtoms::RADIUS_TYPE  radiusType) const;
			float* getPosition();
			int getElementIndex();
			int getResidueNum();
			int getAtomNum();
			const char* getName() const;
			unsigned char getResidueId() const;
			const char* getResidueName() const;
			char getChain();
			float getCharge();
			float getPerAtomHydrophobicity();
			float getPerResidueHydrophobicity();
			void setPosition(const float* position);
			void setPosition(float x, float y, float z);
			void setElementIndex(int elementIndex);
			void setAtomNum(int atNum);
			void setResidueNum(int resNum);
			void setChain(char chain);
			void setBoundaryFlag(const bool flag);
			// WARNING: these functions will mutate the atom
			void setRadius(const float radius);
			void setName(char* name);
			void setResidueId(unsigned char residueId);
			void setCharge(float charge);
			bool isBoundaryAtom();
			PDBParser::GroupOfAtoms::TYPE_OF_RESIDUE getResidueType() const;
			void getColor(float* red, float* green, float* blue, float* alpha);
			void getMinMax(double* minX, double* minY, double* minZ, double* maxX, double* maxY, double* maxZ);
			void getMaximumDistanceFromPoint(double xOrigin, double yOrigin, double zOrigin, double* distance);
			GroupOfAtoms* getGOAptr( void ) { return m_parent; }
			void setGOAptr( GroupOfAtoms *ptr ) { m_parent = ptr; }

			float getDistanceFromRay(CCVOpenGLMath::Ray targetVector);
			double getFunction(PDBParser::GroupOfAtoms::FUNCTIONS FUNCTION);
			bool transform(CCVOpenGLMath::Matrix transformationMatrix);
			MutatedAtomInformation* m_MutatedAtomInformation;

			float m_Position[3];			// x,y,z coordinates
			int m_AtomNum;				// Atom number in PDB
			int m_ResidueNum;			// Residue number in the PDB
			int m_ElementIndex;			// Used to define the atom type. See PDBParser::elementTable for details.
			unsigned int m_UserDefinedColor;
			char m_Chain;				// Chain ID (from PDB)

			bool m_Backbone;			// If true, the atom is a member of the macromolecule's backbone.
			bool m_IsBoundaryAtom;			// If true, atom contributes to the molecular surface

			bool m_donorOrAcceptor;			// false = not a hbond donor/accpetor, true = is a donor/acceptor
			int m_internalHbond;			// 1 = external hbond, 2 = internal, 3 = both internal and external
			GroupOfAtoms *m_parent;

			double m_interfacepropensity;
	};
};

#endif
