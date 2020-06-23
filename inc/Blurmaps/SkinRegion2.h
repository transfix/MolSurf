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
#ifndef CCV_SKIN_REGION_2_H
#define CCV_SKIN_REGION_2_H

#include <Utility/utility.h>
#include <Blurmaps/BlurMapsDataManager.h>
#include <Blurmaps/GOABlur.h>
#include <PDBParser/Atom.h>
#include <PDBParser/CollectionData.h>
#include <PDBParser/ElementInformation.h>
#include <PDBParser/GOAColor.h>
#include <PDBParser/GroupOfAtoms.h>
#include <PDBParser/ResidueInformation.h>
#include <SimpleVolumeData/SimpleVolumeData.h>
#include <SphericalPatchIntersections/SphericalPatch.h>
#include <UsefulMath/TrilinearGrid.h>

using CCVOpenGLMath::TrilinearGrid;
extern double getTime();
namespace PDBParser
{
	class GroupOfAtoms;
	class Atom;
};
class SimpleVolumeData;
const unsigned char EXTERIOR_VOLUME =   0;
const unsigned char SAS_OUTER_SURFACE = 40;
const unsigned char SAS_SURFACE     =  50;
const unsigned char SAS_VOLUME      = 100;
const unsigned char SES_SURFACE     = 150;
const unsigned char SES_VOLUME      = 175;
const unsigned char VDW_SURFACE     = 200;
const unsigned char VDW_VOLUME      = 250;
typedef vector<int> intvec;

class SkinRegion2
{
	public:
		SkinRegion2();
		virtual ~SkinRegion2();
		void setEnlargedRadius(double r);
		bool getSkinRegion(PDBParser::GroupOfAtoms* molecule, int dim1, int dim2, int dim3,
						   double probeRadius, PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType, SimpleVolumeData* vol, int depth,
						   const char* cmapFile =0, bool withColor = false, PDBParser::GroupOfAtoms::GOA_TYPE colorLevel = PDBParser::GroupOfAtoms::ATOM);
	protected:
		bool getColoredAtomList(PDBParser::GroupOfAtoms* molecule, const char* cmapFile, PDBParser::GroupOfAtoms::GOA_TYPE colorLevel);
		bool allocateSpace(bool withColor, SimpleVolumeData* vol);
		void getSasInterior();
		void getBoundaryVerticesCells();
		void getSesSurfaceAndVolume();
		bool isBoundaryCell(bool* sasCellInterior);
		int fillBoundaryIndices();
		void getBoundarySphereIndices();
		void getSubdividedNodesSasInteriorVertices();
		void getSubdividedNodesSasBoundaryCells();
		double inline getDist(double x1, double y1, double z1, double x2, double y2, double z2);
		void getBoundingBox(int x, int y, int z, int* imin, int* jmin, int* kmin, int* imax, int* jmax, int* kmax,
							int iProbeRadius, int jProbeRadius, int kProbeRadius);
		bool intersects(double* center, double radius, int cell);
		bool intersects(double* center, double radius, int iIndex, int jIndex, int kIndex);
		void getSphereRadiusCenter(PDBParser::Atom* at, double* radius, double* center);
		double getDistToSas(int vertx, int verty, int vertz, int neighborx, int neighborj, int neighbork,
							unsigned char* red, unsigned char* green, unsigned char* blue, bool color);
		double getDistToSas2(int vertx, int verty, int vertz, int neighborx, int neighbory, int neighborz);
		void getColors(PDBParser::Atom* atom, unsigned char* red, unsigned char* green, unsigned char* blue);
		bool* m_SasInteriorVertices;
		bool* m_SasBoundaryVertices;
		bool* m_VDWInteriorVertices;
		bool* m_VDWBoundaryVertices;
		bool* m_SasBoundaryCells;
		int* m_IndexOfBoundaryCell;
		intvec* m_SphereIndices;
		bool* m_SubdividedNodeSasInteriorVertex;
		bool* m_SubdividedNodeSasBoundaryCell;
		////   This variable is going to be incredibly overloaded
		////   We need 5 on off bits for 5 boolean variables
		////   And 3 bits to represent 1 of 7 regions.
		////   Also, one of the boolean arrays is smaller than this ( contains cells, not vertices! )
		unsigned char* m_Regions;
		float* m_SES;
		unsigned char* m_Red, *m_Green, *m_Blue;
		vector<PDBParser::Atom*> m_AtomList;
		float m_Min[3],m_Max[3];
		float m_Orig[3],m_Span[3];
		unsigned int m_Dim[3];
		double m_ProbeRadius;
		PDBParser::GroupOfAtoms::RADIUS_TYPE m_RadiusType;
		double enlargedRadius;

		int m_Depth;
		int m_NumBoundaryCells;
		int m_NumCells;
		int m_NumVerts;
		int m_NumSmallCells;
		int m_NumSmallVerts;
};
#endif // CCV_SKIN_REGION_2_H
