#include <Utility/utility.h>
#include <libCG/CoarseGrain/coarseGrain.h>
#include <libCG/CoarseGrain/atom2.h>
#include <libCG/CoarseGrain/groupOfAtoms.h>

using namespace MOLECULE;

CoarseGrain::CoarseGrain()
{
}

CoarseGrain::CoarseGrain(GroupOfAtoms* m)
{
	molecule = m;
}

GroupOfAtoms* CoarseGrain::CoarseGrainMolecule(GroupOfAtoms* m)
{
	int resIndex;
	int head, tail;
	int i;
	GroupOfAtoms* groupOfBeads = new GroupOfAtoms();
	molecule = m;
	head = 0;
	tail = 0;
	resIndex = molecule->m_Atoms[0]->resSeq;
	for(i = 1; i < molecule->numOfAtoms; i++)
	{
		if(molecule->m_Atoms[i]->resSeq == resIndex)
		{
			tail = i;
		}
		else
		{
			atomClustering(groupOfBeads, head, tail);
			resIndex = molecule->m_Atoms[i]->resSeq;
			head = i;
			tail = i;
		}
	}
	atomClustering(groupOfBeads, head, tail);
	groupOfBeads->numOfAtoms = groupOfBeads->m_Atoms.size();
	groupOfBeads->numOfRes = molecule->numOfRes;
	CoarseGrainGeometry(groupOfBeads);
	return groupOfBeads;
}
