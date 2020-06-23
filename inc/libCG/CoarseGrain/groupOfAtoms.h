#ifndef GROUPOFATOMS_H
#define GROUPOFATOMS_H

#include <vector>
#include <libCG/CoarseGrain/surface.h>

namespace MOLECULE
{
	class Atom;
	class GroupOfAtoms
	{
		public:
			GroupOfAtoms();
			GroupOfAtoms(GroupOfAtoms* m);
			virtual ~GroupOfAtoms();
			void addAtom(Atom* currentAtom);
			void getBornRadii(GEOMETRY::Surface* surface);
			void getNeighbors();
			void assignBRtoChild(GroupOfAtoms* mol);
			std::vector<Atom*> m_Atoms;
			int numOfAtoms;
			int numOfRes;
			double total_charge;
			bool* neighbor;
		private:
			bool overlap(int i, int j);
			double dotProduct(double* u, double* v);
	};
}

#endif //GROUPOFATOMS_H
