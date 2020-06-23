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
#if !defined(AFX_GROUPOFATOMS_H__F71426B6_08DA_4EBE_BDA8_747B22B0902F__INCLUDED_)
#define AFX_GROUPOFATOMS_H__F71426B6_08DA_4EBE_BDA8_747B22B0902F__INCLUDED_

#include <Utility/utility.h>
#include <PDBParser/CollectionData.h>
#include <PDBParser/ElementInformation.h>
#include <PDBParser/ExpandableArray.h>
#include <PDBParser/ProteinData.h>
#include <PDBParser/ResidueData.h>
#include <PDBParser/ResidueHeader.h>
#include <PDBParser/ResidueInformation.h>
#include <PDBParser/SecondaryStructureData.h>
#include <UsefulMath/LinearAlgebra.h>
#include <UsefulMath/Matrix.h>
#include <UsefulMath/Quaternion.h>
#include <UsefulMath/Ray.h>
#include <UsefulMath/Vector.h>

namespace PDBParser
{
	class Atom;
	class ResidueHeader;
	class ResidueData;
	class ProteinData;
	class SecondaryStructureData;
	class CollectionData;

	// Different GroupOfAtom types
	const int GENERAL_TYPE      = 0;
	const int ATOM_TYPE         = 1;
	const int RESIDUE_TYPE      = 2;
	const int SS_TYPE           = 3;
	const int PROTEIN_TYPE      = 4;
	const int CHAIN_TYPE        = 5;
	const int COLLECTION_TYPE   = 6;	// the hierarchical level right above molecule

	// Different chain types
	static const int CHAIN_NULL    = 0;
	static const int CHAIN_PROTEIN = 1;
	static const int CHAIN_NUCLEIC = 2;

	// Types of multiresolution
	static const int MULTIRES_BOUNDING_SPHERE = 0;

	// Default colors for secondary structures
	static const int SS_HELIX      = 0;
	static const int SS_SHEET      = 1;
	static const int SS_TURN       = 2;
	static const int SS_NULL       = 3;		// The default null SS when an actual SS is not present.  Used to keep the hierarchy depth consistent.
	static const int TOTAL_SS_TYPE = 4;		// 4 secondary structures total

	// Index into these arrays are the types of secondary structures, such as SS_HELIX, etc.
	static const double ssColorR[TOTAL_SS_TYPE] = { 1.0, 0.3, 0.3, 0.3 };
	static const double ssColorG[TOTAL_SS_TYPE] = { 0.3, 1.0, 0.3, 0.3 };
	static const double ssColorB[TOTAL_SS_TYPE] = { 0.3, 0.3, 1.0, 0.3 };
	static const char* ssString[TOTAL_SS_TYPE] = { "HELIX", "SHEET", "TURN", "NULL" };

	// Default colors for the molecule
	// format: (red, green, blue).
	static const double moleculeColor[3] = { 0.9, 0.6, 0.9 };
	static const char* moleculeString = "MOLECULE";

	typedef struct _hbond {
		Atom *d, *a;
		float energy;
		bool internal;
	}HBOND;

	class GroupOfAtoms
	{
		public:
			GroupOfAtoms();
			GroupOfAtoms(const GroupOfAtoms& goa);	// Copies all substructures recursively
			virtual ~GroupOfAtoms();
			enum TYPE_OF_RESIDUE {AMINO_ACID_RESIDUE, RNA_RESIDUE, UNKNOWN_RESIDUE};
			enum POSITION {INSIDE, OUTSIDE, BOTH};
			enum RADIUS_TYPE { NUMBER_OF_RADII_TYPES = 2, VDW_RADIUS = 0, ROUX_RADIUS};
			static bool intToRadiusType(RADIUS_TYPE* radiusType, int intRadiusType);
			enum FUNCTIONS { NUMBER_OF_FUNCTIONS = 4, ELECTRON_DENSITY = 0, PER_ATOM_HYDROPHOBICITY, PER_RESIDUE_HYDROPHOBICITY, ELECTROSTATIC_POTENTIAL };
			static bool intToFunctionType(FUNCTIONS* functionType, int intFunctionType);
			enum GOA_TYPE { NUMBER_OF_GOA_TYPES = 6, ATOM = 0, RESIDUE, SECONDARY_STRUCTURE, CHAIN, PROTEIN, COLLECTION };
			static bool intToGOAType(GOA_TYPE* goaType, int intGOAType);
			void addAtom(int serial, int elementIndex, int resNum, char* name, unsigned char resId, char chain, double x, double y, double z, double radius, float charge, bool boundaryAtom);
			void buildHierarchicalInformation(double xmid, double ymid, double zmid);
			void setGOAColor(double r, double g, double b);	// Sets the color for the entire group of atoms
			void setResidueColor(int resID);
			void setSSColor(int ssID);
			void setMoleculeColor(void);
			void createMultires(double xmid, double ymid, double zmid, int multires_type);
			void createSingleBoundingSphere(double xmid, double ymid, double zmid);
			void allocateMemory(int numAtoms);
			bool getAngles(unsigned int n, double* psi, double* phi, unsigned char* resType);
			void normalizeResidues();
			void fillUpAtoms(double* x, double* y, double* z, int* index);
			bool fillUpAtoms(double* x, double* y, double* z, int* index, int strandNumber);
			bool getNumberOfAtomsRecursive(int* count);
			void getAtomsRecursive(double* posRadius, int maxLength, int* curPosition);
			bool getAttributes(double* points, double* radii, double* min, double* max, double* weights, double* maxRadius, PDBParser::GroupOfAtoms::FUNCTIONS FUNCTION, int GOA_level, PDBParser::GroupOfAtoms::RADIUS_TYPE radius_type, int maxLength);
			double getFunction(PDBParser::GroupOfAtoms::FUNCTIONS function);
			bool getNumberOfAtoms(int* count, int strandNumber);
			bool getEndPoints(CCVOpenGLMath::Vector* v1, CCVOpenGLMath::Vector* v2);
			bool getEndPoints(CCVOpenGLMath::Vector* v1, CCVOpenGLMath::Vector* v2, int strandNumber);
			bool getEndPoints(CCVOpenGLMath::Vector* v1, CCVOpenGLMath::Vector* v2, int startingResidue, int endingResidue);
			void transformAtoms(double* transformation);
			bool addChain(GroupOfAtoms* chain);
			//----------------------------------------------
			// Functions for build hierarchical information
			bool createHelicesForResidues(std::vector<double>* m_HelixGeometry, int startingResidue, int length );
			void createSS();
			void createSSs();
			//----------------------------------------------
			void getMinMax(double* minX, double* minY, double* minZ, double* maxX, double* maxY, double* maxZ);
			void getMaximumDistanceFromPoint(double xOrigin, double yOrigin, double zOrigin, double* distance);
			bool intersects(double* min, double* max);
			bool within(double* min, double* max);
			int getNumberOfGOAs(PDBParser::GroupOfAtoms::GOA_TYPE goaType);
			void getAtoms(double* min, double* max, std::vector<double>* interiorAtomList, std::vector<double>* exteriorAtomList, POSITION position);
			Atom* getSelectedAtom(CCVOpenGLMath::Ray targetVector);
			bool applyColorMap(const char* fileName);
			void printColors(FILE* fp);
			static GroupOfAtoms* getEmptyHierarchy();
			PDBParser::GroupOfAtoms::TYPE_OF_RESIDUE getResidueType()  const;
			Atom* getAtom(const char* atomName);
			bool transform(CCVOpenGLMath::Matrix transformationMatrix);
			static bool isValidGOA(PDBParser::GroupOfAtoms* molecule, double overlapAllowed);
			int getMassCenter(double* x, double* y, double* z);
			
			int getNumberOfReceptorAtoms(void) { return m_nAtomsReceptor; }
			void setNumberOfReceptorAtoms(int val) { m_nAtomsReceptor = val; }
			int getNumberOfLigandAtoms(void) { return m_nAtomsLigand; }
			void setNumberOfLigandAtoms(int val) { m_nAtomsLigand = val; }
			std::vector<HBOND> *getHBond(void) { return &m_HBond; }

			// Hierarchy components of the GOA
			std::vector<GroupOfAtoms*> m_SubGroups;
			std::vector<Atom*> m_Atoms;
			ResidueData* m_ResidueData;
			SecondaryStructureData* m_SecondaryStructureData;
			ProteinData* m_ProteinData;
			CollectionData* m_CollectionData;
			// Hydrogen-bond
			int m_nAtomsReceptor;
			int m_nAtomsLigand;
			std::vector<HBOND> m_HBond;

			// For CHAINS:
			// Not part of the hierarchy, but a list of all of the atoms that are members of
			// the backbone of a chain.  This array is used for chain type GOA's.
			// These atoms are not part of the hierarchy because they already exist in another
			// GOA as part of the hierarchy.  This array is just a quick access set of
			// pointers that speeds up computation.
			// This backbone does not contain the O atom for a protein, and contains only N, CA and C'.

			// Used for the multiresolution representation of the GroupOfAtoms
			std::vector<Atom*> m_MultiresAtoms;		// Multi-res. representaion, may contain 1 or more atoms
			// Color used for this group of atoms when it is represented in multiresolution
			float multi_r;
			float multi_g;
			float multi_b;

			// WARNING: the parser may not be up to date yet to reflect the comments below
			int type;				// GENERAL_TYPE, ATOM_TYPE, RESIDUE_TYPE, SS_TYPE, CHAIN_TYPE, PROTEIN_TYPE

			int resID;				// For a residue, the integer index for the residue in the table 'residueTable'
			// For a secondary structure, it is one of the SS types.
			// For a chain, it is the character identifier of the chain, cast as an int

			int m_ResSeq;			// The actual residue sequence number. It was strangly missing in the first few versions
			// of this code !
			int resPos;				// For a residue, contains the atom serial number of the first atom in the current
			// residue (the serial number is the 5-digit number that appears right after ATOM).
			// These serial numbers appear to be unique identifiers throughout a PDB file, so
			// they can be used to  quickly linked back to the part of the file thatt created
			// them.
			//
			// For a chain, it is one of the CHAIN types (CHAIN_NULL, CHAIN_PROTEIN, or CHAIN_NUCLEIC)

			ResidueHeader* resHeader;	// The header info for a residue


			// Used to indicate the boundaries for a secondary structure
			ResidueHeader* startRes;
			ResidueHeader* endRes;

			float m_Min[3], m_Max[3];

		protected:
			bool add(std::vector<double>* list, double x, double y, double z, double r);
			void getMinMaxRecursive(double* minX, double* minY, double* minZ, double* maxX, double* maxY, double* maxZ);
			void getAtomsRecursive(double* min, double* max, std::vector<double>* interiorAtomList, std::vector<double>* exteriorAtomList, POSITION position);
			void getAttributesRecursive(double* points, double* radii, double* min, double* max, double* weights, double* maxRadius, PDBParser::GroupOfAtoms::FUNCTIONS FUNCTION, int GOA_level, PDBParser::GroupOfAtoms::RADIUS_TYPE radius_type, int* current_atom_position, int maxLength);
			void recursiveGetClosestAtom(Atom** atom, double* distance, CCVOpenGLMath::Ray targetVector);
			bool intersects(CCVOpenGLMath::Ray targetVector);
			bool getNumberOfGOAsRecursive(int* count, PDBParser::GroupOfAtoms::GOA_TYPE goaType);
	};

	static double virusTransformations[] =
	{
		-0.80902, -0.50000, 0.30902, -0.50000, 0.30902, -0.80902, 0.30902, -0.80902, -0.50000, 0.00000, 1.00000, 0.00000, 0.00000, 0.00000, -1.00000, -1.00000, 0.00000, 0.00000, 0.80902, -0.50000, -0.30902, -0.50000, -0.30902, -0.80902, 0.30902, 0.80902, -0.50000, 0.50000, 0.30902, -0.80902, -0.30902, -0.80902, -0.50000, -0.80902, 0.50000, -0.30902, -0.30902, -0.80902, -0.50000, -0.80902, 0.50000, -0.30902, 0.50000, 0.30902, -0.80902, -0.80902, 0.50000, -0.30902, 0.50000, 0.30902, -0.80902, -0.30902, -0.80902, -0.50000, -0.80902, -0.50000, -0.30902, 0.50000, -0.30902, -0.80902, 0.30902, -0.80902, 0.50000, -0.30902, 0.80902, -0.50000, -0.80902, -0.50000, -0.30902, -0.50000, 0.30902, 0.80902, 0.50000, -0.30902, -0.80902, -0.30902, 0.80902, -0.50000, 0.80902, 0.50000, 0.30902, 0.00000, 0.00000, -1.00000, -1.00000, 0.00000, 0.00000, 0.00000, 1.00000, 0.00000, -0.50000, -0.30902, -0.80902, 0.30902, 0.80902, -0.50000, 0.80902, -0.50000, -0.30902, -0.50000, 0.30902, -0.80902, 0.30902, -0.80902, -0.50000, -0.80902, -0.50000, 0.30902, -0.30902, -0.80902, -0.50000, 0.80902, -0.50000, 0.30902, -0.50000, -0.30902, 0.80902, -0.80902, 0.50000, -0.30902, -0.50000, -0.30902, 0.80902, 0.30902, 0.80902, 0.50000, 0.50000, 0.30902, -0.80902, 0.30902, 0.80902, 0.50000, 0.80902, -0.50000, 0.30902, -0.50000, 0.30902, -0.80902, -0.30902, 0.80902, 0.50000, 0.80902, 0.50000, -0.30902, 0.00000, 0.00000, -1.00000, 1.00000, 0.00000, 0.00000, 0.00000, -1.00000, 0.00000, -0.50000, -0.30902, -0.80902, -0.30902, -0.80902, 0.50000, -0.80902, 0.50000, 0.30902, 0.00000, -1.00000, 0.00000, 0.00000, 0.00000, 1.00000, -1.00000, 0.00000, 0.00000, -0.80902, 0.50000, 0.30902, 0.50000, 0.30902, 0.80902, 0.30902, 0.80902, -0.50000, 0.80902, 0.50000, -0.30902, 0.50000, -0.30902, 0.80902, 0.30902, -0.80902, -0.50000, -0.30902, 0.80902, -0.50000, 0.80902, 0.50000, 0.30902, 0.50000, -0.30902, -0.80902, 0.50000, -0.30902, -0.80902, 0.30902, -0.80902, 0.50000, -0.80902, -0.50000, -0.30902, -0.80902, -0.50000, -0.30902, -0.50000, 0.30902, 0.80902, -0.30902, 0.80902, -0.50000, -0.30902, -0.80902, 0.50000, -0.80902, 0.50000, 0.30902, -0.50000, -0.30902, -0.80902, -0.30902, 0.80902, 0.50000, 0.80902, 0.50000, -0.30902, -0.50000, 0.30902, -0.80902, 1.00000, 0.00000, 0.00000, 0.00000, -1.00000, 0.00000, 0.00000, 0.00000, -1.00000, 0.30902, 0.80902, -0.50000, 0.80902, -0.50000, -0.30902, -0.50000, -0.30902, -0.80902, 0.30902, -0.80902, -0.50000, -0.80902, -0.50000, 0.30902, -0.50000, 0.30902, -0.80902, -1.00000, 0.00000, 0.00000, 0.00000, 1.00000, 0.00000, 0.00000, 0.00000, -1.00000, 0.80902, 0.50000, -0.30902, -0.50000, 0.30902, -0.80902, -0.30902, 0.80902, 0.50000, 0.00000, -1.00000, 0.00000, 0.00000, 0.00000, -1.00000, 1.00000, 0.00000, 0.00000, -0.80902, 0.50000, 0.30902, -0.50000, -0.30902, -0.80902, -0.30902, -0.80902, 0.50000, -0.50000, -0.30902, 0.80902, -0.30902, -0.80902, -0.50000, 0.80902, -0.50000, 0.30902, 0.30902, 0.80902, 0.50000, -0.80902, 0.50000, -0.30902, -0.50000, -0.30902, 0.80902, 0.80902, -0.50000, 0.30902, 0.50000, 0.30902, -0.80902, 0.30902, 0.80902, 0.50000, 0.80902, 0.50000, 0.30902, 0.50000, -0.30902, -0.80902, -0.30902, 0.80902, -0.50000, 0.30902, -0.80902, 0.50000, -0.80902, -0.50000, -0.30902, 0.50000, -0.30902, -0.80902, -0.50000, 0.30902, 0.80902, -0.30902, 0.80902, -0.50000, -0.80902, -0.50000, -0.30902, 0.00000, 0.00000, 1.00000, -1.00000, 0.00000, 0.00000, 0.00000, -1.00000, 0.00000, 0.50000, 0.30902, 0.80902, 0.30902, 0.80902, -0.50000, -0.80902, 0.50000, 0.30902, 0.50000, -0.30902, 0.80902, 0.30902, -0.80902, -0.50000, 0.80902, 0.50000, -0.30902, 0.30902, 0.80902, 0.50000, 0.80902, -0.50000, 0.30902, 0.50000, 0.30902, -0.80902, 0.80902, -0.50000, 0.30902, -0.50000, -0.30902, 0.80902, -0.30902, -0.80902, -0.50000, -0.50000, -0.30902, 0.80902, 0.30902, 0.80902, 0.50000, -0.80902, 0.50000, -0.30902, 0.50000, -0.30902, 0.80902, -0.30902, 0.80902, 0.50000, -0.80902, -0.50000, 0.30902, 0.00000, 0.00000, 1.00000, 1.00000, 0.00000, 0.00000, 0.00000, 1.00000, 0.00000, 0.50000, 0.30902, 0.80902, -0.30902, -0.80902, 0.50000, 0.80902, -0.50000, -0.30902, 0.00000, 1.00000, 0.00000, 0.00000, 0.00000, 1.00000, 1.00000, 0.00000, 0.00000, 0.80902, -0.50000, -0.30902, 0.50000, 0.30902, 0.80902, -0.30902, -0.80902, 0.50000, -0.80902, -0.50000, 0.30902, 0.50000, -0.30902, 0.80902, -0.30902, 0.80902, 0.50000, 0.30902, -0.80902, 0.50000, 0.80902, 0.50000, 0.30902, -0.50000, 0.30902, 0.80902, -0.50000, 0.30902, 0.80902, 0.30902, -0.80902, 0.50000, 0.80902, 0.50000, 0.30902, 0.80902, 0.50000, 0.30902, -0.50000, 0.30902, 0.80902, 0.30902, -0.80902, 0.50000, 0.30902, 0.80902, -0.50000, -0.80902, 0.50000, 0.30902, 0.50000, 0.30902, 0.80902, 0.30902, -0.80902, -0.50000, 0.80902, 0.50000, -0.30902, 0.50000, -0.30902, 0.80902, -1.00000, 0.00000, 0.00000, 0.00000, -1.00000, 0.00000, 0.00000, 0.00000, 1.00000, -0.30902, -0.80902, 0.50000, 0.80902, -0.50000, -0.30902, 0.50000, 0.30902, 0.80902, -0.30902, 0.80902, 0.50000, -0.80902, -0.50000, 0.30902, 0.50000, -0.30902, 0.80902, 1.00000, 0.00000, 0.00000, 0.00000, 1.00000, 0.00000, 0.00000, 0.00000, 1.00000
	};

};

#endif
