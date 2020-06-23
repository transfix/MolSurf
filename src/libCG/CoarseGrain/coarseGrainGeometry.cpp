#include <Utility/utility.h>
#include <libCG/CoarseGrain/coarseGrain.h>
#include <libCG/CoarseGrain/atom2.h>
#include <libCG/CoarseGrain/groupOfAtoms.h>

using namespace MOLECULE;

void CoarseGrain::atomClustering(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	generateBackboneBead(GroupOfBeads, head, tail);
	if(strcmp(molecule->m_Atoms[head]->resName, "ARG") == 0)
	{
		generateARGBead(GroupOfBeads, head, tail);
		return;
	}
	if(strcmp(molecule->m_Atoms[head]->resName, "GLN") == 0)
	{
		generateGLNBead(GroupOfBeads, head, tail);
		return;
	}
	if(strcmp(molecule->m_Atoms[head]->resName, "GLU") == 0)
	{
		generateGLUBead(GroupOfBeads, head, tail);
		return;
	}
	if(strcmp(molecule->m_Atoms[head]->resName, "HIS") == 0)
	{
		generateHISBead(GroupOfBeads, head, tail);
		return;
	}
	if(strcmp(molecule->m_Atoms[head]->resName, "LYS") == 0)
	{
		generateLYSBead(GroupOfBeads, head, tail);
		return;
	}
	if(strcmp(molecule->m_Atoms[head]->resName, "MET") == 0)
	{
		generateMETBead(GroupOfBeads, head, tail);
		return;
	}
	if(strcmp(molecule->m_Atoms[head]->resName, "PHE") == 0)
	{
		generatePHEBead(GroupOfBeads, head, tail);
		return;
	}
	if(strcmp(molecule->m_Atoms[head]->resName, "TRP") == 0)
	{
		generateTRPBead(GroupOfBeads, head, tail);
		return;
	}
	if(strcmp(molecule->m_Atoms[head]->resName, "TYR") == 0)
	{
		generateTYRBead(GroupOfBeads, head, tail);
		return;
	}
	if(strcmp(molecule->m_Atoms[head]->resName, "GLY") == 0)
	{
		generateGLYBead(GroupOfBeads, head, tail);
		return;
	}
	// new clustering
	if(strcmp(molecule->m_Atoms[head]->resName, "VAL") == 0)
	{
		generateVALBead(GroupOfBeads, head, tail);
		return;
	}
	if(strcmp(molecule->m_Atoms[head]->resName, "LEU") == 0)
	{
		generateLEUBead(GroupOfBeads, head, tail);
		return;
	}
	if(strcmp(molecule->m_Atoms[head]->resName, "ILE") == 0)
	{
		generateILEBead(GroupOfBeads, head, tail);
		return;
	}
	if(strcmp(molecule->m_Atoms[head]->resName, "ASN") == 0)
	{
		generateASNBead(GroupOfBeads, head, tail);
		return;
	}
	if(strcmp(molecule->m_Atoms[head]->resName, "ASP") == 0)
	{
		generateASPBead(GroupOfBeads, head, tail);
		return;
	}
	if(strcmp(molecule->m_Atoms[head]->resName, "ALA") == 0)
	{
		generateALABead(GroupOfBeads, head, tail);
		return;
	}
	if(strcmp(molecule->m_Atoms[head]->resName, "SER") == 0)
	{
		generateSERBead(GroupOfBeads, head, tail);
		return;
	}
	if(strcmp(molecule->m_Atoms[head]->resName, "THR") == 0)
	{
		generateTHRBead(GroupOfBeads, head, tail);
		return;
	}
	if(strcmp(molecule->m_Atoms[head]->resName, "CYS") == 0)
	{
		generateCYSBead(GroupOfBeads, head, tail);
		return;
	}
	if(strcmp(molecule->m_Atoms[head]->resName, "PRO") == 0)
	{
		generatePROBead(GroupOfBeads, head, tail);
		return;
	}
}

void CoarseGrain::generateBackboneBead(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	char resName[5];
	strcpy(resName, molecule->m_Atoms[head]->resName);
	int resSeq = molecule->m_Atoms[head]->resSeq;
	char chainID = molecule->m_Atoms[head]->chainID;
	int i;
	Atom* bead_bb1 = new Atom();
	Atom* bead_bb2 = new Atom();
	Atom* bead_bb3 = new Atom();
	vector<int> BackboneList1;
	vector<int> BackboneList2;
	vector<int> BackboneList3;
	if(strstr(resName, "GLY") !=NULL)
	{
		for(i = head; i <= tail; i++)
		{
			BackboneList1.push_back(i);
		}
	}
	else
		for(i = head; i <= tail; i++)
		{
			if(isBackboneAtom(i))
			{
				BackboneList1.push_back(i);
			}
		}
	bead_bb1->child = new int[BackboneList1.size()];
	for(i = 0; i < BackboneList1.size(); i++)
	{
		bead_bb1->child[i] = BackboneList1[i];
	}
	bead_bb1->numOfChildren = BackboneList1.size();
	bead_bb1->setResName(resName);
	bead_bb1->setResIndex(resSeq);
	bead_bb1->setChainID(chainID);
	GroupOfBeads->m_Atoms.push_back(bead_bb1);
}

void CoarseGrain::generateARGBead(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	char resName[5];
	strcpy(resName, "ARG");
	int resSeq = molecule->m_Atoms[head]->resSeq;
	char chainID = molecule->m_Atoms[head]->chainID;
	int i;
	Atom* bead_sc1 = new Atom();
	Atom* bead_sc2 = new Atom();
	vector<int> SidechainList1;
	vector<int> SidechainList2;
	for(i = head; i <= tail; i++)
	{
		if(!isBackboneAtom(i))
		{
			char atomName[10];
			strcpy(atomName, molecule->m_Atoms[i]->atomName);
			if(strstr(atomName, "B") != NULL ||
					strstr(atomName, "G") != NULL)
			{
				SidechainList1.push_back(i);
			}
			else
			{
				SidechainList2.push_back(i);
			}
		}
	}
	bead_sc1->child = new int[SidechainList1.size()];
	bead_sc2->child = new int[SidechainList2.size()];
	for(i = 0; i < SidechainList1.size(); i++)
	{
		bead_sc1->child[i] = SidechainList1[i];
	}
	for(i = 0; i < SidechainList2.size(); i++)
	{
		bead_sc2->child[i] = SidechainList2[i];
	}
	bead_sc1->numOfChildren = SidechainList1.size();
	bead_sc2->numOfChildren = SidechainList2.size();
	bead_sc1->setResName(resName);
	bead_sc2->setResName(resName);
	bead_sc1->setResIndex(resSeq);
	bead_sc2->setResIndex(resSeq);
	bead_sc1->setChainID(chainID);
	bead_sc2->setChainID(chainID);
	GroupOfBeads->m_Atoms.push_back(bead_sc1);
	GroupOfBeads->m_Atoms.push_back(bead_sc2);
}

void CoarseGrain::generateGLNBead(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	char resName[5];
	strcpy(resName, "GLN");
	int resSeq = molecule->m_Atoms[head]->resSeq;
	char chainID = molecule->m_Atoms[head]->chainID;
	int i;
	Atom* bead_sc1 = new Atom();
	Atom* bead_sc2 = new Atom();
	vector<int> SidechainList1;
	vector<int> SidechainList2;
	for(i = head; i <= tail; i++)
	{
		if(!isBackboneAtom(i))
		{
			char atomName[10];
			strcpy(atomName, molecule->m_Atoms[i]->atomName);
			if(strstr(atomName, "B") != NULL ||
					strstr(atomName, "G") != NULL)
			{
				SidechainList1.push_back(i);
			}
			else
			{
				SidechainList2.push_back(i);
			}
		}
	}
	bead_sc1->child = new int[SidechainList1.size()];
	bead_sc2->child = new int[SidechainList2.size()];
	for(i = 0; i < SidechainList1.size(); i++)
	{
		bead_sc1->child[i] = SidechainList1[i];
	}
	for(i = 0; i < SidechainList2.size(); i++)
	{
		bead_sc2->child[i] = SidechainList2[i];
	}
	bead_sc1->numOfChildren = SidechainList1.size();
	bead_sc2->numOfChildren = SidechainList2.size();
	bead_sc1->setResName(resName);
	bead_sc2->setResName(resName);
	bead_sc1->setResIndex(resSeq);
	bead_sc2->setResIndex(resSeq);
	bead_sc1->setChainID(chainID);
	bead_sc2->setChainID(chainID);
	GroupOfBeads->m_Atoms.push_back(bead_sc1);
	GroupOfBeads->m_Atoms.push_back(bead_sc2);
}

void CoarseGrain::generateGLUBead(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	char resName[5];
	strcpy(resName, "GLU");
	int resSeq = molecule->m_Atoms[head]->resSeq;
	char chainID = molecule->m_Atoms[head]->chainID;
	int i;
	Atom* bead_sc1 = new Atom();
	Atom* bead_sc2 = new Atom();
	vector<int> SidechainList1;
	vector<int> SidechainList2;
	for(i = head; i <= tail; i++)
	{
		if(!isBackboneAtom(i))
		{
			char atomName[10];
			strcpy(atomName, molecule->m_Atoms[i]->atomName);
			if(strstr(atomName, "B") != NULL ||
					strstr(atomName, "G") != NULL)
			{
				SidechainList1.push_back(i);
			}
			else
			{
				SidechainList2.push_back(i);
			}
		}
	}
	bead_sc1->child = new int[SidechainList1.size()];
	bead_sc2->child = new int[SidechainList2.size()];
	for(i = 0; i < SidechainList1.size(); i++)
	{
		bead_sc1->child[i] = SidechainList1[i];
	}
	for(i = 0; i < SidechainList2.size(); i++)
	{
		bead_sc2->child[i] = SidechainList2[i];
	}
	bead_sc1->numOfChildren = SidechainList1.size();
	bead_sc2->numOfChildren = SidechainList2.size();
	bead_sc1->setResName(resName);
	bead_sc2->setResName(resName);
	bead_sc1->setResIndex(resSeq);
	bead_sc2->setResIndex(resSeq);
	bead_sc1->setChainID(chainID);
	bead_sc2->setChainID(chainID);
	GroupOfBeads->m_Atoms.push_back(bead_sc1);
	GroupOfBeads->m_Atoms.push_back(bead_sc2);
}

void CoarseGrain::generateHISBead(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	char resName[5];
	strcpy(resName, "HIS");
	int resSeq = molecule->m_Atoms[head]->resSeq;
	char chainID = molecule->m_Atoms[head]->chainID;
	int i;
	Atom* bead_sc1 = new Atom();
	Atom* bead_sc2 = new Atom();
	vector<int> SidechainList1;
	vector<int> SidechainList2;
	for(i = head; i <= tail; i++)
	{
		if(!isBackboneAtom(i))
		{
			char atomName[10];
			strcpy(atomName, molecule->m_Atoms[i]->atomName);
			if(strstr(atomName, "B") != NULL ||
					strstr(atomName, "G") != NULL ||
					strstr(atomName, "ND1") != NULL ||
					strstr(atomName, "HD1") != NULL)
			{
				SidechainList1.push_back(i);
			}
			else
			{
				SidechainList2.push_back(i);
			}
		}
	}
	bead_sc1->child = new int[SidechainList1.size()];
	bead_sc2->child = new int[SidechainList2.size()];
	for(i = 0; i < SidechainList1.size(); i++)
	{
		bead_sc1->child[i] = SidechainList1[i];
	}
	for(i = 0; i < SidechainList2.size(); i++)
	{
		bead_sc2->child[i] = SidechainList2[i];
	}
	bead_sc1->numOfChildren = SidechainList1.size();
	bead_sc2->numOfChildren = SidechainList2.size();
	bead_sc1->setResName(resName);
	bead_sc2->setResName(resName);
	bead_sc1->setResIndex(resSeq);
	bead_sc2->setResIndex(resSeq);
	bead_sc1->setChainID(chainID);
	bead_sc2->setChainID(chainID);
	GroupOfBeads->m_Atoms.push_back(bead_sc1);
	GroupOfBeads->m_Atoms.push_back(bead_sc2);
}

void CoarseGrain::generateLYSBead(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	char resName[5];
	strcpy(resName, "LYS");
	int resSeq = molecule->m_Atoms[head]->resSeq;
	char chainID = molecule->m_Atoms[head]->chainID;
	int i;
	Atom* bead_sc1 = new Atom();
	Atom* bead_sc2 = new Atom();
	vector<int> SidechainList1;
	vector<int> SidechainList2;
	for(i = head; i <= tail; i++)
	{
		if(!isBackboneAtom(i))
		{
			char atomName[10];
			strcpy(atomName, molecule->m_Atoms[i]->atomName);
			if(strstr(atomName, "B") != NULL ||
					strstr(atomName, "G") != NULL)
			{
				SidechainList1.push_back(i);
			}
			else
			{
				SidechainList2.push_back(i);
			}
		}
	}
	bead_sc1->child = new int[SidechainList1.size()];
	bead_sc2->child = new int[SidechainList2.size()];
	for(i = 0; i < SidechainList1.size(); i++)
	{
		bead_sc1->child[i] = SidechainList1[i];
	}
	for(i = 0; i < SidechainList2.size(); i++)
	{
		bead_sc2->child[i] = SidechainList2[i];
	}
	bead_sc1->numOfChildren = SidechainList1.size();
	bead_sc2->numOfChildren = SidechainList2.size();
	bead_sc1->setResName(resName);
	bead_sc2->setResName(resName);
	bead_sc1->setResIndex(resSeq);
	bead_sc2->setResIndex(resSeq);
	bead_sc1->setChainID(chainID);
	bead_sc2->setChainID(chainID);
	GroupOfBeads->m_Atoms.push_back(bead_sc1);
	GroupOfBeads->m_Atoms.push_back(bead_sc2);
}

void CoarseGrain::generateMETBead(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	char resName[5];
	strcpy(resName, "MET");
	int resSeq = molecule->m_Atoms[head]->resSeq;
	char chainID = molecule->m_Atoms[head]->chainID;
	int i;
	Atom* bead_sc1 = new Atom();
	Atom* bead_sc2 = new Atom();
	vector<int> SidechainList1;
	vector<int> SidechainList2;
	for(i = head; i <= tail; i++)
	{
		if(!isBackboneAtom(i))
		{
			char atomName[10];
			strcpy(atomName, molecule->m_Atoms[i]->atomName);
			if(strstr(atomName, "B") != NULL ||
					strstr(atomName, "G") != NULL)
			{
				SidechainList1.push_back(i);
			}
			else
			{
				SidechainList2.push_back(i);
			}
		}
	}
	bead_sc1->child = new int[SidechainList1.size()];
	bead_sc2->child = new int[SidechainList2.size()];
	for(i = 0; i < SidechainList1.size(); i++)
	{
		bead_sc1->child[i] = SidechainList1[i];
	}
	for(i = 0; i < SidechainList2.size(); i++)
	{
		bead_sc2->child[i] = SidechainList2[i];
	}
	bead_sc1->numOfChildren = SidechainList1.size();
	bead_sc2->numOfChildren = SidechainList2.size();
	bead_sc1->setResName(resName);
	bead_sc2->setResName(resName);
	bead_sc1->setResIndex(resSeq);
	bead_sc2->setResIndex(resSeq);
	bead_sc1->setChainID(chainID);
	bead_sc2->setChainID(chainID);
	GroupOfBeads->m_Atoms.push_back(bead_sc1);
	GroupOfBeads->m_Atoms.push_back(bead_sc2);
}

void CoarseGrain::generatePHEBead(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	char resName[5];
	strcpy(resName, "PHE");
	int resSeq = molecule->m_Atoms[head]->resSeq;
	char chainID = molecule->m_Atoms[head]->chainID;
	int i;
	Atom* bead_sc1 = new Atom();
	Atom* bead_sc2 = new Atom();
	vector<int> SidechainList1;
	vector<int> SidechainList2;
	for(i = head; i <= tail; i++)
	{
		if(!isBackboneAtom(i))
		{
			char atomName[10];
			strcpy(atomName, molecule->m_Atoms[i]->atomName);
			if(strstr(atomName, "B") != NULL ||
					strstr(atomName, "G") != NULL ||
					strstr(atomName, "D") != NULL)
			{
				SidechainList1.push_back(i);
			}
			else
			{
				SidechainList2.push_back(i);
			}
		}
	}
	bead_sc1->child = new int[SidechainList1.size()];
	bead_sc2->child = new int[SidechainList2.size()];
	for(i = 0; i < SidechainList1.size(); i++)
	{
		bead_sc1->child[i] = SidechainList1[i];
	}
	for(i = 0; i < SidechainList2.size(); i++)
	{
		bead_sc2->child[i] = SidechainList2[i];
	}
	bead_sc1->numOfChildren = SidechainList1.size();
	bead_sc2->numOfChildren = SidechainList2.size();
	bead_sc1->setResName(resName);
	bead_sc2->setResName(resName);
	bead_sc1->setResIndex(resSeq);
	bead_sc2->setResIndex(resSeq);
	bead_sc1->setChainID(chainID);
	bead_sc2->setChainID(chainID);
	GroupOfBeads->m_Atoms.push_back(bead_sc1);
	GroupOfBeads->m_Atoms.push_back(bead_sc2);
}

void CoarseGrain::generateTRPBead(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	char resName[5];
	strcpy(resName, "TRP");
	int resSeq = molecule->m_Atoms[head]->resSeq;
	char chainID = molecule->m_Atoms[head]->chainID;
	int i;
	Atom* bead_sc1 = new Atom();
	Atom* bead_sc2 = new Atom();
	vector<int> SidechainList1;
	vector<int> SidechainList2;
	for(i = head; i <= tail; i++)
	{
		if(!isBackboneAtom(i))
		{
			char atomName[10];
			strcpy(atomName, molecule->m_Atoms[i]->atomName);
			if(strstr(atomName, "B") != NULL ||
					strstr(atomName, "G") != NULL ||
					strstr(atomName, "D1") != NULL ||
					strstr(atomName, "E1") != NULL)
			{
				SidechainList1.push_back(i);
			}
			else
			{
				SidechainList2.push_back(i);
			}
		}
	}
	bead_sc1->child = new int[SidechainList1.size()];
	bead_sc2->child = new int[SidechainList2.size()];
	for(i = 0; i < SidechainList1.size(); i++)
	{
		bead_sc1->child[i] = SidechainList1[i];
	}
	for(i = 0; i < SidechainList2.size(); i++)
	{
		bead_sc2->child[i] = SidechainList2[i];
	}
	bead_sc1->numOfChildren = SidechainList1.size();
	bead_sc2->numOfChildren = SidechainList2.size();
	bead_sc1->setResName(resName);
	bead_sc2->setResName(resName);
	bead_sc1->setResIndex(resSeq);
	bead_sc2->setResIndex(resSeq);
	bead_sc1->setChainID(chainID);
	bead_sc2->setChainID(chainID);
	GroupOfBeads->m_Atoms.push_back(bead_sc1);
	GroupOfBeads->m_Atoms.push_back(bead_sc2);
}

void CoarseGrain::generateTYRBead(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	char resName[5];
	strcpy(resName, "TYR");
	int resSeq = molecule->m_Atoms[head]->resSeq;
	char chainID = molecule->m_Atoms[head]->chainID;
	int i;
	Atom* bead_sc1 = new Atom();
	Atom* bead_sc2 = new Atom();
	vector<int> SidechainList1;
	vector<int> SidechainList2;
	for(i = head; i <= tail; i++)
	{
		if(!isBackboneAtom(i))
		{
			char atomName[10];
			strcpy(atomName, molecule->m_Atoms[i]->atomName);
			if(strstr(atomName, "B") != NULL ||
					strstr(atomName, "G") != NULL ||
					strstr(atomName, "D") != NULL)
			{
				SidechainList1.push_back(i);
			}
			else
			{
				SidechainList2.push_back(i);
			}
		}
	}
	bead_sc1->child = new int[SidechainList1.size()];
	bead_sc2->child = new int[SidechainList2.size()];
	for(i = 0; i < SidechainList1.size(); i++)
	{
		bead_sc1->child[i] = SidechainList1[i];
	}
	for(i = 0; i < SidechainList2.size(); i++)
	{
		bead_sc2->child[i] = SidechainList2[i];
	}
	bead_sc1->numOfChildren = SidechainList1.size();
	bead_sc2->numOfChildren = SidechainList2.size();
	bead_sc1->setResName(resName);
	bead_sc2->setResName(resName);
	bead_sc1->setResIndex(resSeq);
	bead_sc2->setResIndex(resSeq);
	bead_sc1->setChainID(chainID);
	bead_sc2->setChainID(chainID);
	GroupOfBeads->m_Atoms.push_back(bead_sc1);
	GroupOfBeads->m_Atoms.push_back(bead_sc2);
}

void CoarseGrain::generateGLYBead(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	char resName[5];
	strcpy(resName, "GLY");
	int resSeq = molecule->m_Atoms[head]->resSeq;
	char chainID = molecule->m_Atoms[head]->chainID;
}

/* new clustering begins */
void CoarseGrain::generateVALBead(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	char resName[5];
	strcpy(resName, "VAL");
	int resSeq = molecule->m_Atoms[head]->resSeq;
	char chainID = molecule->m_Atoms[head]->chainID;
	int i;
	Atom* bead_sc1 = new Atom();
	Atom* bead_sc2 = new Atom();
	vector<int> SidechainList1;
	vector<int> SidechainList2;
	for(i = head; i <= tail; i++)
	{
		if(!isBackboneAtom(i))
		{
			char atomName[10];
			strcpy(atomName, molecule->m_Atoms[i]->atomName);
			if(strstr(atomName, "B") != NULL)
			{
				SidechainList1.push_back(i);
			}
			else
			{
				SidechainList2.push_back(i);
			}
		}
	}
	bead_sc1->child = new int[SidechainList1.size()];
	bead_sc2->child = new int[SidechainList2.size()];
	for(i = 0; i < SidechainList1.size(); i++)
	{
		bead_sc1->child[i] = SidechainList1[i];
	}
	for(i = 0; i < SidechainList2.size(); i++)
	{
		bead_sc2->child[i] = SidechainList2[i];
	}
	bead_sc1->numOfChildren = SidechainList1.size();
	bead_sc2->numOfChildren = SidechainList2.size();
	bead_sc1->setResName(resName);
	bead_sc2->setResName(resName);
	bead_sc1->setResIndex(resSeq);
	bead_sc2->setResIndex(resSeq);
	bead_sc1->setChainID(chainID);
	bead_sc2->setChainID(chainID);
	GroupOfBeads->m_Atoms.push_back(bead_sc1);
	GroupOfBeads->m_Atoms.push_back(bead_sc2);
}


void CoarseGrain::generateLEUBead(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	char resName[5];
	strcpy(resName, "LEU");
	int resSeq = molecule->m_Atoms[head]->resSeq;
	char chainID = molecule->m_Atoms[head]->chainID;
	int i;
	Atom* bead_sc1 = new Atom();
	Atom* bead_sc2 = new Atom();
	vector<int> SidechainList1;
	vector<int> SidechainList2;
	for(i = head; i <= tail; i++)
	{
		if(!isBackboneAtom(i))
		{
			char atomName[10];
			strcpy(atomName, molecule->m_Atoms[i]->atomName);
			if(strstr(atomName, "B") != NULL ||
					strstr(atomName, "G") != NULL)
			{
				SidechainList1.push_back(i);
			}
			else
			{
				SidechainList2.push_back(i);
			}
		}
	}
	bead_sc1->child = new int[SidechainList1.size()];
	bead_sc2->child = new int[SidechainList2.size()];
	for(i = 0; i < SidechainList1.size(); i++)
	{
		bead_sc1->child[i] = SidechainList1[i];
	}
	for(i = 0; i < SidechainList2.size(); i++)
	{
		bead_sc2->child[i] = SidechainList2[i];
	}
	bead_sc1->numOfChildren = SidechainList1.size();
	bead_sc2->numOfChildren = SidechainList2.size();
	bead_sc1->setResName(resName);
	bead_sc2->setResName(resName);
	bead_sc1->setResIndex(resSeq);
	bead_sc2->setResIndex(resSeq);
	bead_sc1->setChainID(chainID);
	bead_sc2->setChainID(chainID);
	GroupOfBeads->m_Atoms.push_back(bead_sc1);
	GroupOfBeads->m_Atoms.push_back(bead_sc2);
}


void CoarseGrain::generateILEBead(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	char resName[5];
	strcpy(resName, "ILE");
	int resSeq = molecule->m_Atoms[head]->resSeq;
	char chainID = molecule->m_Atoms[head]->chainID;
	int i;
	Atom* bead_sc1 = new Atom();
	Atom* bead_sc2 = new Atom();
	vector<int> SidechainList1;
	vector<int> SidechainList2;
	for(i = head; i <= tail; i++)
	{
		if(!isBackboneAtom(i))
		{
			char atomName[10];
			strcpy(atomName, molecule->m_Atoms[i]->atomName);
			if(strstr(atomName, "B") != NULL ||
					strstr(atomName, "G2") != NULL)
			{
				SidechainList1.push_back(i);
			}
			else
			{
				SidechainList2.push_back(i);
			}
		}
	}
	bead_sc1->child = new int[SidechainList1.size()];
	bead_sc2->child = new int[SidechainList2.size()];
	for(i = 0; i < SidechainList1.size(); i++)
	{
		bead_sc1->child[i] = SidechainList1[i];
	}
	for(i = 0; i < SidechainList2.size(); i++)
	{
		bead_sc2->child[i] = SidechainList2[i];
	}
	bead_sc1->numOfChildren = SidechainList1.size();
	bead_sc2->numOfChildren = SidechainList2.size();
	bead_sc1->setResName(resName);
	bead_sc2->setResName(resName);
	bead_sc1->setResIndex(resSeq);
	bead_sc2->setResIndex(resSeq);
	bead_sc1->setChainID(chainID);
	bead_sc2->setChainID(chainID);
	GroupOfBeads->m_Atoms.push_back(bead_sc1);
	GroupOfBeads->m_Atoms.push_back(bead_sc2);
}

void CoarseGrain::generateASNBead(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	char resName[5];
	strcpy(resName, "ASN");
	int resSeq = molecule->m_Atoms[head]->resSeq;
	char chainID = molecule->m_Atoms[head]->chainID;
	int i;
	Atom* bead_sc1 = new Atom();
	Atom* bead_sc2 = new Atom();
	vector<int> SidechainList1;
	vector<int> SidechainList2;
	for(i = head; i <= tail; i++)
	{
		if(!isBackboneAtom(i))
		{
			char atomName[10];
			strcpy(atomName, molecule->m_Atoms[i]->atomName);
			if(strstr(atomName, "B") != NULL)
			{
				SidechainList1.push_back(i);
			}
			else
			{
				SidechainList2.push_back(i);
			}
		}
	}
	bead_sc1->child = new int[SidechainList1.size()];
	bead_sc2->child = new int[SidechainList2.size()];
	for(i = 0; i < SidechainList1.size(); i++)
	{
		bead_sc1->child[i] = SidechainList1[i];
	}
	for(i = 0; i < SidechainList2.size(); i++)
	{
		bead_sc2->child[i] = SidechainList2[i];
	}
	bead_sc1->numOfChildren = SidechainList1.size();
	bead_sc2->numOfChildren = SidechainList2.size();
	bead_sc1->setResName(resName);
	bead_sc2->setResName(resName);
	bead_sc1->setResIndex(resSeq);
	bead_sc2->setResIndex(resSeq);
	bead_sc1->setChainID(chainID);
	bead_sc2->setChainID(chainID);
	GroupOfBeads->m_Atoms.push_back(bead_sc1);
	GroupOfBeads->m_Atoms.push_back(bead_sc2);
}

void CoarseGrain::generateASPBead(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	char resName[5];
	strcpy(resName, "ASP");
	int resSeq = molecule->m_Atoms[head]->resSeq;
	char chainID = molecule->m_Atoms[head]->chainID;
	int i;
	Atom* bead_sc1 = new Atom();
	Atom* bead_sc2 = new Atom();
	vector<int> SidechainList1;
	vector<int> SidechainList2;
	for(i = head; i <= tail; i++)
	{
		if(!isBackboneAtom(i))
		{
			char atomName[10];
			strcpy(atomName, molecule->m_Atoms[i]->atomName);
			if(strstr(atomName, "B") != NULL)
			{
				SidechainList1.push_back(i);
			}
			else
			{
				SidechainList2.push_back(i);
			}
		}
	}
	bead_sc1->child = new int[SidechainList1.size()];
	bead_sc2->child = new int[SidechainList2.size()];
	for(i = 0; i < SidechainList1.size(); i++)
	{
		bead_sc1->child[i] = SidechainList1[i];
	}
	for(i = 0; i < SidechainList2.size(); i++)
	{
		bead_sc2->child[i] = SidechainList2[i];
	}
	bead_sc1->numOfChildren = SidechainList1.size();
	bead_sc2->numOfChildren = SidechainList2.size();
	bead_sc1->setResName(resName);
	bead_sc2->setResName(resName);
	bead_sc1->setResIndex(resSeq);
	bead_sc2->setResIndex(resSeq);
	bead_sc1->setChainID(chainID);
	bead_sc2->setChainID(chainID);
	GroupOfBeads->m_Atoms.push_back(bead_sc1);
	GroupOfBeads->m_Atoms.push_back(bead_sc2);
}


void CoarseGrain::generateALABead(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	char resName[5];
	strcpy(resName, "ALA");
	int resSeq = molecule->m_Atoms[head]->resSeq;
	char chainID = molecule->m_Atoms[head]->chainID;
	int i;
	Atom* bead_sc1 = new Atom();
	vector<int> SidechainList1;
	for(i = head; i <= tail; i++)
	{
		if(!isBackboneAtom(i))
		{
			char atomName[10];
			strcpy(atomName, molecule->m_Atoms[i]->atomName);
			SidechainList1.push_back(i);
		}
	}
	bead_sc1->child = new int[SidechainList1.size()];
	for(i = 0; i < SidechainList1.size(); i++)
	{
		bead_sc1->child[i] = SidechainList1[i];
	}
	bead_sc1->numOfChildren = SidechainList1.size();
	bead_sc1->setResName(resName);
	bead_sc1->setResIndex(resSeq);
	bead_sc1->setChainID(chainID);
	GroupOfBeads->m_Atoms.push_back(bead_sc1);
}


void CoarseGrain::generateSERBead(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	char resName[5];
	strcpy(resName, "SER");
	int resSeq = molecule->m_Atoms[head]->resSeq;
	char chainID = molecule->m_Atoms[head]->chainID;
	int i;
	Atom* bead_sc1 = new Atom();
	Atom* bead_sc2 = new Atom();
	vector<int> SidechainList1;
	vector<int> SidechainList2;
	for(i = head; i <= tail; i++)
	{
		if(!isBackboneAtom(i))
		{
			char atomName[10];
			strcpy(atomName, molecule->m_Atoms[i]->atomName);
			if(strstr(atomName, "G") != NULL)
			{
				SidechainList1.push_back(i);
			}
			else
			{
				SidechainList2.push_back(i);
			}
		}
	}
	bead_sc1->child = new int[SidechainList1.size()];
	bead_sc2->child = new int[SidechainList2.size()];
	for(i = 0; i < SidechainList1.size(); i++)
	{
		bead_sc1->child[i] = SidechainList1[i];
	}
	for(i = 0; i < SidechainList2.size(); i++)
	{
		bead_sc2->child[i] = SidechainList2[i];
	}
	bead_sc1->numOfChildren = SidechainList1.size();
	bead_sc2->numOfChildren = SidechainList2.size();
	bead_sc1->setResName(resName);
	bead_sc2->setResName(resName);
	bead_sc1->setResIndex(resSeq);
	bead_sc2->setResIndex(resSeq);
	bead_sc1->setChainID(chainID);
	bead_sc2->setChainID(chainID);
	GroupOfBeads->m_Atoms.push_back(bead_sc1);
	GroupOfBeads->m_Atoms.push_back(bead_sc2);
}


void CoarseGrain::generateTHRBead(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	char resName[5];
	strcpy(resName, "THR");
	int resSeq = molecule->m_Atoms[head]->resSeq;
	char chainID = molecule->m_Atoms[head]->chainID;
	int i;
	Atom* bead_sc1 = new Atom();
	Atom* bead_sc2 = new Atom();
	vector<int> SidechainList1;
	vector<int> SidechainList2;
	for(i = head; i <= tail; i++)
	{
		if(!isBackboneAtom(i))
		{
			char atomName[10];
			strcpy(atomName, molecule->m_Atoms[i]->atomName);
			if(strstr(atomName, "G1") != NULL)
			{
				SidechainList1.push_back(i);
			}
			else
			{
				SidechainList2.push_back(i);
			}
		}
	}
	bead_sc1->child = new int[SidechainList1.size()];
	bead_sc2->child = new int[SidechainList2.size()];
	for(i = 0; i < SidechainList1.size(); i++)
	{
		bead_sc1->child[i] = SidechainList1[i];
	}
	for(i = 0; i < SidechainList2.size(); i++)
	{
		bead_sc2->child[i] = SidechainList2[i];
	}
	bead_sc1->numOfChildren = SidechainList1.size();
	bead_sc2->numOfChildren = SidechainList2.size();
	bead_sc1->setResName(resName);
	bead_sc2->setResName(resName);
	bead_sc1->setResIndex(resSeq);
	bead_sc2->setResIndex(resSeq);
	bead_sc1->setChainID(chainID);
	bead_sc2->setChainID(chainID);
	GroupOfBeads->m_Atoms.push_back(bead_sc1);
	GroupOfBeads->m_Atoms.push_back(bead_sc2);
}


void CoarseGrain::generateCYSBead(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	char resName[5];
	strcpy(resName, "CYS");
	int resSeq = molecule->m_Atoms[head]->resSeq;
	char chainID = molecule->m_Atoms[head]->chainID;
	int i;
	Atom* bead_sc1 = new Atom();
	Atom* bead_sc2 = new Atom();
	vector<int> SidechainList1;
	vector<int> SidechainList2;
	for(i = head; i <= tail; i++)
	{
		if(!isBackboneAtom(i))
		{
			char atomName[10];
			strcpy(atomName, molecule->m_Atoms[i]->atomName);
			if(strstr(atomName, "S") != NULL)
			{
				SidechainList1.push_back(i);
			}
			else
			{
				SidechainList2.push_back(i);
			}
		}
	}
	bead_sc1->child = new int[SidechainList1.size()];
	bead_sc2->child = new int[SidechainList2.size()];
	for(i = 0; i < SidechainList1.size(); i++)
	{
		bead_sc1->child[i] = SidechainList1[i];
	}
	for(i = 0; i < SidechainList2.size(); i++)
	{
		bead_sc2->child[i] = SidechainList2[i];
	}
	bead_sc1->numOfChildren = SidechainList1.size();
	bead_sc2->numOfChildren = SidechainList2.size();
	bead_sc1->setResName(resName);
	bead_sc2->setResName(resName);
	bead_sc1->setResIndex(resSeq);
	bead_sc2->setResIndex(resSeq);
	bead_sc1->setChainID(chainID);
	bead_sc2->setChainID(chainID);
	GroupOfBeads->m_Atoms.push_back(bead_sc1);
	GroupOfBeads->m_Atoms.push_back(bead_sc2);
}


void CoarseGrain::generatePROBead(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	char resName[5];
	strcpy(resName, "PRO");
	int resSeq = molecule->m_Atoms[head]->resSeq;
	char chainID = molecule->m_Atoms[head]->chainID;
	int i;
	Atom* bead_sc1 = new Atom();
	Atom* bead_sc2 = new Atom();
	vector<int> SidechainList1;
	vector<int> SidechainList2;
	for(i = head; i <= tail; i++)
	{
		if(!isBackboneAtom(i))
		{
			char atomName[10];
			strcpy(atomName, molecule->m_Atoms[i]->atomName);
			if(strstr(atomName, "B") != NULL ||
					strstr(atomName, "G") != NULL)
			{
				SidechainList1.push_back(i);
			}
			else
			{
				SidechainList2.push_back(i);
			}
		}
	}
	bead_sc1->child = new int[SidechainList1.size()];
	bead_sc2->child = new int[SidechainList2.size()];
	for(i = 0; i < SidechainList1.size(); i++)
	{
		bead_sc1->child[i] = SidechainList1[i];
	}
	for(i = 0; i < SidechainList2.size(); i++)
	{
		bead_sc2->child[i] = SidechainList2[i];
	}
	bead_sc1->numOfChildren = SidechainList1.size();
	bead_sc2->numOfChildren = SidechainList2.size();
	bead_sc1->setResName(resName);
	bead_sc2->setResName(resName);
	bead_sc1->setResIndex(resSeq);
	bead_sc2->setResIndex(resSeq);
	bead_sc1->setChainID(chainID);
	bead_sc2->setChainID(chainID);
	GroupOfBeads->m_Atoms.push_back(bead_sc1);
	GroupOfBeads->m_Atoms.push_back(bead_sc2);
}

void CoarseGrain::generateNormalBead(GroupOfAtoms* GroupOfBeads, int head, int tail)
{
	char resName[5];
	strcpy(resName, molecule->m_Atoms[head]->resName);
	int resSeq = molecule->m_Atoms[head]->resSeq;
	char chainID = molecule->m_Atoms[head]->chainID;
	int i;
	Atom* bead_bb = new Atom();
	Atom* bead_sc = new Atom();
	vector<int> BackboneList;
	vector<int> SidechainList;
	for(i = head; i <= tail; i++)
	{
		if(isBackboneAtom(i))
		{
			BackboneList.push_back(i);
		}
		else
		{
			SidechainList.push_back(i);
		}
	}
	bead_bb->child = new int[BackboneList.size()];
	bead_sc->child = new int[SidechainList.size()];
	for(i = 0; i < BackboneList.size(); i++)
	{
		bead_bb->child[i] = BackboneList[i];
	}
	for(i = 0; i < SidechainList.size(); i++)
	{
		bead_sc->child[i] = SidechainList[i];
	}
	bead_bb->numOfChildren = BackboneList.size();
	bead_sc->numOfChildren = SidechainList.size();
	bead_bb->setResName(resName);
	bead_sc->setResName(resName);
	bead_bb->setResIndex(resSeq);
	bead_sc->setResIndex(resSeq);
	bead_bb->setChainID(chainID);
	bead_sc->setChainID(chainID);
	GroupOfBeads->m_Atoms.push_back(bead_bb);
	GroupOfBeads->m_Atoms.push_back(bead_sc);
}

bool CoarseGrain::isBackboneAtom(int a)
{
	char atomName[10];
	strcpy(atomName, molecule->m_Atoms[a]->atomName);
	if(strcmp(atomName, " N  ") == 0 ||
			strcmp(atomName, " CA ") == 0 ||
			strcmp(atomName, " C  ") == 0 ||
			strcmp(atomName, " O  ") == 0 ||
			strcmp(atomName, " H  ") == 0 ||
			strcmp(atomName, " HA ") == 0 ||
			strcmp(atomName, " H2 ") == 0 ||
			strcmp(atomName, " H3 ") == 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void CoarseGrain::CoarseGrainGeometry(GroupOfAtoms* GroupOfBeads)
{
	int i;
	for(i = 0; i < GroupOfBeads->numOfAtoms; i++)
	{
		double center[3];
		double radius;
		int* childlist = GroupOfBeads->m_Atoms[i]->child;
		int numOfChild = GroupOfBeads->m_Atoms[i]->numOfChildren;
		getGeometryCenter(childlist, center, numOfChild);
		radius = getInitialRadius(childlist, center, numOfChild);
		GroupOfBeads->m_Atoms[i]->setPosition(center);
		GroupOfBeads->m_Atoms[i]->setDefinedRadius(radius);
	}
}

void CoarseGrain::getMassCenter(int* AtomList, double* coord, int numOfChild)
{
	int i;
	double c[3], total_mass;
	c[0] = c[1] = c[2] = total_mass = 0.0;
	for(i = 0; i < numOfChild; i++)
	{
		double p[3], mass;
		int k = AtomList[i];
		p[0] = molecule->m_Atoms[k]->position[0];
		p[1] = molecule->m_Atoms[k]->position[1];
		p[2] = molecule->m_Atoms[k]->position[2];
		mass  = molecule->m_Atoms[k]->mass;
		c[0] += mass*p[0];
		c[1] += mass*p[1];
		c[2] += mass*p[2];
		total_mass += mass;
	}
	coord[0] = c[0]/total_mass;
	coord[1] = c[1]/total_mass;
	coord[2] = c[2]/total_mass;
}

void CoarseGrain::getGeometryCenter(int* AtomList, double* coord, int numOfChild)
{
	int i;
	double c[3];
	c[0] = c[1] = c[2] = 0.0;
	for(i = 0; i < numOfChild; i++)
	{
		double p[3], mass;
		int k = AtomList[i];
		p[0] = molecule->m_Atoms[k]->position[0];
		p[1] = molecule->m_Atoms[k]->position[1];
		p[2] = molecule->m_Atoms[k]->position[2];
		c[0] += p[0];
		c[1] += p[1];
		c[2] += p[2];
	}
	coord[0] = c[0]/numOfChild;
	coord[1] = c[1]/numOfChild;
	coord[2] = c[2]/numOfChild;
}

double CoarseGrain::getInitialRadius(int* AtomList, double* center, int numOfChild)
{
	int i;
	double p[3], r, d, max_d;
	max_d = 0.0;
	for(i = 0; i < numOfChild; i++)
	{
		p[0] = molecule->m_Atoms[AtomList[i]]->position[0];
		p[1] = molecule->m_Atoms[AtomList[i]]->position[1];
		p[2] = molecule->m_Atoms[AtomList[i]]->position[2];
		r = molecule->m_Atoms[AtomList[i]]->radius;
		d = sqrt((p[0] - center[0])*(p[0] - center[0]) +
				 (p[1] - center[1])*(p[1] - center[1]) +
				 (p[2] - center[2])*(p[2] - center[2])) + r;
		if(d > max_d)
		{
			max_d = d;
		}
	}
	return max_d;
}
