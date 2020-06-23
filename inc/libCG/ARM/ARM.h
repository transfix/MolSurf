#ifndef ARM_H
#define ARM_H

#include <Blurmaps/BlurMapsDataManager.h>
#include <Geometry/Geometry.h>
#include <GOAFileTypes/GOALoader.h>
#include <PDBParser/Atom.h>
#include <PDBParser/GOAFileIO.h>
#include <PDBParser/GroupOfAtoms.h>
#include <PDBParser/parserPDBtoGOA.h>

#define ARMLevel 2
#define ARMUpperLevel 1
#define MAXUNITS 100


typedef struct BeadStruct
{
	int index;
	int numchilds;
	int childs[MAXUNITS];
} BeadStructure;

class ARM
{
	private:
		std::vector<BeadStruct> armHierarchy;
		std::vector<PDBParser::Atom*> arm[2]; // the atomic level and the first level of adaptive resolution model

	public:
		ARM(char* pqrfiles[ARMLevel]);
		// first level: wenqi's coarse-grained model
		// second level: amino-acid (future work)
		~ARM();
		static PDBParser::GroupOfAtoms*  read_pdb(const char* fname);
		static void getAtomList(std::vector<PDBParser::Atom*>& atomList, const char* fname);
		void setHierarchicalInformation(char* hinfofile);
		void getAtomListOfBead(int index, std::vector<PDBParser::Atom*>& atomList);
		Geometry* getSurfaceOfBead(int index, int size);
		void printAtomList(std::vector<PDBParser::Atom*>& atomList, char* outfilename);
		void getOptimizedCG(std::vector<PDBParser::Atom*>& cg, int size, int iteration, char* outfilename);
};

#endif
