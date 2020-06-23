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
#ifndef CCV_CURVATURESGRID_H
#define CCV_CURVATURESGRID_H

#include <Utility/utility.h>
#include <Blurmaps/CurvaturesGridVoxel.h>
#include <PDBParser/Atom.h>
#include <PDBParser/GroupOfAtoms.h>
#include <UsefulMath/LinearAlgebra.h>

// arand commented...
//using PDBParser::Atom;
//using PDBParser::GroupOfAtoms;

namespace PDBParser
{
	class GroupOfAtoms;
};
class CurvaturesGridVoxel;

class CurvaturesGrid
{
	public:
		CurvaturesGrid();
		virtual ~CurvaturesGrid();

		bool getCurvatures(double* HandK, double* normals, double* k1Vec, double* k2Vec, PDBParser::GroupOfAtoms* molecule, int numberOfPoints, float* points, double blobbyness, int numberOfGridDivisions, double maxFunctionError, int radiusType, unsigned int level);
		bool write(int numberOfPoints, double* HandK, double* normals, double* k1Vec, double* k2Vec, const char* filename);

	protected:
		void createGrid(PDBParser::GroupOfAtoms* molecule, int numberOfPoints, float* points, double blobbyness, int numberOfGridDivisions, double maxFunctionError, int radiusType, unsigned int level);
		void fillAtomsRecursively(PDBParser::GroupOfAtoms* molecule, double overlapExtent, int radiusType, unsigned int level);
		void getMinMax(PDBParser::GroupOfAtoms* molecule, int numberOfPoints, float* points);
		int getIndex(double x, double y, double z);
		void getIndices(double x, double y, double z, int* gridx, int* gridy, int* gridz);
		void evalCurvature(double* phiX, double* phiY, double* phiZ, double* phiXX,
						   double* phiYY, double* phiZZ,  double* phiXY, double* phiXZ,
						   double* phiYZ,  double xc, double yc, double zc, double rad,
						   double x, double y, double z, int vertexNumber);

		double m_Min[3], m_Max[3];
		CurvaturesGridVoxel* grid;
		int m_NumberOfGridDivisions;
		double m_Blobbyness;
		int m_NumberOfPoints;
};

#endif
