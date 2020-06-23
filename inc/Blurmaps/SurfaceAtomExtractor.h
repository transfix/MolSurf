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
#ifndef CCV_SURFACE_ATOM_EXTRACTOR_H
#define CCV_SURFACE_ATOM_EXTRACTOR_H

#include <Utility/utility.h>
#include <PDBParser/ResidueInformation.h>
#include <PDBParser/GroupOfAtoms.h>
#include <PDBParser/Atom.h>
#include <Blurmaps/GOABlur.h>
#include <PDBParser/ElementInformation.h>
#include <PDBParser/CollectionData.h>
#include <UsefulMath/TrilinearGrid.h>
#include <Blurmaps/BlurMapsDataManager.h>
#include <SphericalPatchIntersections/SphericalPatch.h>

namespace PDBParser
{
	class GroupOfAtoms;
	class Atom;
};

typedef vector<int> intvec;

class SurfaceAtomExtractor
{
	public:
		SurfaceAtomExtractor();
		virtual ~SurfaceAtomExtractor();

		bool getBoundaryAtoms(PDBParser::GroupOfAtoms* molecule, int dim1, int dim2, int dim3,
							  double probeRadius, PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType, bool** boundaryAtom);

	protected:

		void getSasInterior();
		void getBoundaryVerticesCells();
		void getSesSurfaceAndVolume();

		bool isBoundaryCell(bool* sasCellInterior);

		int fillBoundaryIndices();
		void getBoundarySphereIndices(bool** boundaryAtom);
		void getBoundingBox(int x, int y, int z, int* imin, int* jmin, int* kmin, int* imax, int* jmax, int* kmax,
							int iProbeRadius, int jProbeRadius, int kProbeRadius);
		bool intersects(double* center, double radius, int cell);
		bool intersects(double* center, double radius, int iIndex, int jIndex, int kIndex);

		void getSphereRadiusCenter(PDBParser::Atom* at, double* radius, double* center);

		bool* m_SasInteriorVertices;
		bool* m_SasBoundaryVertices;
		bool* m_SasBoundaryCells;

		int* m_IndexOfBoundaryCell;

		vector<PDBParser::Atom*> m_AtomList;
		float m_Min[3],m_Max[3];
		float m_Orig[3],m_Span[3];
		unsigned int m_Dim[3];
		double m_ProbeRadius;
		PDBParser::GroupOfAtoms::RADIUS_TYPE m_RadiusType;
		int m_NumBoundaryCells;

		int m_NumCells;
		int m_NumVerts;
};

#endif
