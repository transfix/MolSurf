#include <Utility/utility.h>
#include <libCG/CoarseGrain/coarseGrain.h>
#include <libCG/CoarseGrain/atom2.h>
#include <libCG/CoarseGrain/groupOfAtoms.h>
#include <libCG/CoarseGrain/optimizer.h>

using namespace MOLECULE;
using namespace CGOptim;

void CoarseGrain::CoarseGrainChargeGB(GroupOfAtoms* groupOfBeads)
{
	// total charge constraint
	int N = groupOfBeads->numOfAtoms;
	int atomIndex;
	int i, j;
	double Q_total = 0.0;
	double E[N], R[N], Q[N];
	for(i = 0; i < N; i++)
	{
		Q[i] = 0.0;
		for(j = 0; j < groupOfBeads->m_Atoms[i]->numOfChildren; j++)
		{
			atomIndex  = groupOfBeads->m_Atoms[i]->child[j];
			Q[i] += molecule->m_Atoms[atomIndex]->charge;
			Q_total += molecule->m_Atoms[atomIndex]->charge;
		}
		R[i] = groupOfBeads->m_Atoms[i]->Born_radius;
		E[i] = getBeadSelfEnergy(groupOfBeads, i);
	}
	Optimizer* opt = new Optimizer(E, R, Q, N, Q_total);
	opt->optimizeCharge();
	// assign charges
	for(i = 0; i < N; i++)
	{
		groupOfBeads->m_Atoms[i]->setCharge(opt->Q[i]);
	}
	// check
	double error = 0.0;
	for(i = 0; i < N; i++)
	{
		error += (E[i]-Q[i]*Q[i]/R[i])*(E[i]-Q[i]*Q[i]/R[i]);
	}
	printf("%f \n", sqrt(error)); //E[i], Q[i]*Q[i]/R[i]);
}

void CoarseGrain::CoarseGrainChargeGB2(GroupOfAtoms* groupOfBeads)
{
	// total charge constraint
	int N = groupOfBeads->numOfAtoms;
	int atomIndex;
	int i, j;
	double Q_total = 0.0;
	double AE, R[N], Q[N], X[N], Y[N], Z[N];
	for(i = 0; i < N; i++)
	{
		Q[i] = 0.0;
		for(j = 0; j < groupOfBeads->m_Atoms[i]->numOfChildren; j++)
		{
			atomIndex  = groupOfBeads->m_Atoms[i]->child[j];
			Q[i] += molecule->m_Atoms[atomIndex]->charge;
			Q_total += molecule->m_Atoms[atomIndex]->charge;
		}
		X[i] = groupOfBeads->m_Atoms[i]->position[0];
		Y[i] = groupOfBeads->m_Atoms[i]->position[1];
		Z[i] = groupOfBeads->m_Atoms[i]->position[2];
		R[i] = groupOfBeads->m_Atoms[i]->Born_radius;
	}
	AE = getAtomicGBEnergy();
	Optimizer* opt = new Optimizer(AE, R, Q, X, Y, Z, N, Q_total);
	opt->optimizeCharge();
	// assign charges
	for(i = 0; i < N; i++)
	{
		groupOfBeads->m_Atoms[i]->setCharge(opt->Q[i]);
	}
	// check
	double ecol = getAtomicColumbEnergy();
	double err1 = 0.0;
	double err2 = 0.0;
	for(i = 0; i < N; i++)
		for(j = 0; j < N; j++)
		{
			double rij2 = (X[i]-X[j])*(X[i]-X[j]) + (Y[i]-Y[j])*(Y[i]-Y[j]) + (Z[i]-Z[j])*(Z[i]-Z[j]);
			err1 += groupOfBeads->m_Atoms[i]->charge*groupOfBeads->m_Atoms[j]->charge / sqrt(rij2 + R[i]*R[j]*exp(-1.0*rij2/(4.0*R[i]*R[j])));
		}
	printf("AE GB = %f CE GB = %f\n", AE, err1); //E[i], Q[i]*Q[i]/R[i]);
	for(i = 0; i < N; i++)
		for(j = i+1; j < N; j++)
		{
			double rij2 = (X[i]-X[j])*(X[i]-X[j]) + (Y[i]-Y[j])*(Y[i]-Y[j]) + (Z[i]-Z[j])*(Z[i]-Z[j]);
			err2 += groupOfBeads->m_Atoms[i]->charge*groupOfBeads->m_Atoms[j]->charge / sqrt(rij2);
		}
	printf("AE Col = %f CE Col = %f\n", ecol, err2);
}

double CoarseGrain::getBeadSelfEnergy(GroupOfAtoms* groupOfBeads, int index)
{
	int i, j, atomi, atomj;
	double qi, qj, xi[3], xj[3], Ri, Rj;
	double rij2;
	int N = groupOfBeads->m_Atoms[index]->numOfChildren;
	double f = 0.0;
	for(i = 0; i < N; i++)
	{
		atomi = groupOfBeads->m_Atoms[index]->child[i];
		qi = molecule->m_Atoms[atomi]->charge;
		//		printf("q%d %f\n", i, qi);
		Ri = molecule->m_Atoms[atomi]->Born_radius;
		xi[0] = molecule->m_Atoms[atomi]->position[0];
		xi[1] = molecule->m_Atoms[atomi]->position[1];
		xi[2] = molecule->m_Atoms[atomi]->position[2];
		for(j = 0; j < N; j++)
		{
			atomj = groupOfBeads->m_Atoms[index]->child[j];
			qj = molecule->m_Atoms[atomj]->charge;
			//			printf("q%d %f\n", j, qj);
			Rj = molecule->m_Atoms[atomj]->Born_radius;
			xj[0] = molecule->m_Atoms[atomj]->position[0];
			xj[1] = molecule->m_Atoms[atomj]->position[1];
			xj[2] = molecule->m_Atoms[atomj]->position[2];
			rij2 = (xi[0] - xj[0])*(xi[0] - xj[0]) + (xi[1] - xj[1])*(xi[1] - xj[1]) + (xi[2] - xj[2])*(xi[2] - xj[2]);
			f += qi*qj / sqrt(rij2 + Ri*Rj * exp(-1.0*rij2 / (4.0*Ri*Rj)));
		}
	}
	return f;
}

double CoarseGrain::getAtomicGBEnergy()
{
	int i, j;
	double qi, qj, xi[3], xj[3], Ri, Rj;
	double rij2;
	double f= 0.0;
	for(i = 0; i < molecule->numOfAtoms; i++)
	{
		qi = molecule->m_Atoms[i]->charge;
		Ri = molecule->m_Atoms[i]->Born_radius;
		xi[0] = molecule->m_Atoms[i]->position[0];
		xi[1] = molecule->m_Atoms[i]->position[1];
		xi[2] = molecule->m_Atoms[i]->position[2];
		for(j = 0; j < molecule->numOfAtoms; j++)
		{
			qj = molecule->m_Atoms[j]->charge;
			Rj = molecule->m_Atoms[j]->Born_radius;
			xj[0] = molecule->m_Atoms[j]->position[0];
			xj[1] = molecule->m_Atoms[j]->position[1];
			xj[2] = molecule->m_Atoms[j]->position[2];
			rij2 = (xi[0] - xj[0])*(xi[0] - xj[0]) + (xi[1] - xj[1])*(xi[1] - xj[1]) + (xi[2] - xj[2])*(xi[2] - xj[2]);
			f += qi*qj / sqrt(rij2 + Ri*Rj * exp(-1.0*rij2 / (4.0*Ri*Rj)));
		}
	}
	return f;
}

double CoarseGrain::getAtomicColumbEnergy()
{
	int i, j;
	double qi, qj, xi[3], xj[3], Ri, Rj;
	double rij2;
	double f= 0.0;
	for(i = 0; i < molecule->numOfAtoms; i++)
	{
		qi = molecule->m_Atoms[i]->charge;
		Ri = molecule->m_Atoms[i]->Born_radius;
		xi[0] = molecule->m_Atoms[i]->position[0];
		xi[1] = molecule->m_Atoms[i]->position[1];
		xi[2] = molecule->m_Atoms[i]->position[2];
		for(j = i+1; j < molecule->numOfAtoms; j++)
		{
			qj = molecule->m_Atoms[j]->charge;
			Rj = molecule->m_Atoms[j]->Born_radius;
			xj[0] = molecule->m_Atoms[j]->position[0];
			xj[1] = molecule->m_Atoms[j]->position[1];
			xj[2] = molecule->m_Atoms[j]->position[2];
			rij2 = (xi[0] - xj[0])*(xi[0] - xj[0]) + (xi[1] - xj[1])*(xi[1] - xj[1]) + (xi[2] - xj[2])*(xi[2] - xj[2]);
			f += qi*qj / sqrt(rij2);
		}
	}
	return f;
}
