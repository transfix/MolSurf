#ifndef COARSEGRAIN_H
#define COARSEGRAIN_H

#include <vector>

namespace MOLECULE
{
	class GroupOfAtoms;
	class CoarseGrain
	{
		public:
			CoarseGrain();
			CoarseGrain(GroupOfAtoms* m);
			GroupOfAtoms* CoarseGrainMolecule(GroupOfAtoms* m);
			void CoarseGrainCharge(GroupOfAtoms* GroupOfBeads);
			void CoarseGrainChargeGB(GroupOfAtoms* groupOfBeads);
			void CoarseGrainChargeGB2(GroupOfAtoms* groupOfBeads);
			void CoarseGrainGeometry(GroupOfAtoms* GroupOfBeads);
		private:
			GroupOfAtoms* molecule;
			typedef struct _indexRange
			{
				int min[3];
				int max[3];
			} IndexRange;

			void atomClustering(GroupOfAtoms* GroupOfBeads, int head, int tail);

			void generateNormalBead(GroupOfAtoms* GroupOfBeads, int head, int tail);
			void generateBackboneBead(GroupOfAtoms* GroupOfBeads, int head, int tail);
			void generateARGBead(GroupOfAtoms* GroupOfBeads, int head, int tail);
			void generateGLNBead(GroupOfAtoms* GroupOfBeads, int head, int tail);
			void generateGLUBead(GroupOfAtoms* GroupOfBeads, int head, int tail);
			void generateHISBead(GroupOfAtoms* GroupOfBeads, int head, int tail);
			void generateLYSBead(GroupOfAtoms* GroupOfBeads, int head, int tail);
			void generateMETBead(GroupOfAtoms* GroupOfBeads, int head, int tail);
			void generatePHEBead(GroupOfAtoms* GroupOfBeads, int head, int tail);
			void generateTRPBead(GroupOfAtoms* GroupOfBeads, int head, int tail);
			void generateTYRBead(GroupOfAtoms* GroupOfBeads, int head, int tail);
			void generateGLYBead(GroupOfAtoms* GroupOfBeads, int head, int tail);
			// new clustering
			void generateVALBead(GroupOfAtoms* GroupOfBeads, int head, int tail);
			void generateLEUBead(GroupOfAtoms* GroupOfBeads, int head, int tail);
			void generateILEBead(GroupOfAtoms* GroupOfBeads, int head, int tail);
			void generateASNBead(GroupOfAtoms* GroupOfBeads, int head, int tail);
			void generateASPBead(GroupOfAtoms* GroupOfBeads, int head, int tail);
			void generateALABead(GroupOfAtoms* GroupOfBeads, int head, int tail);
			void generateSERBead(GroupOfAtoms* GroupOfBeads, int head, int tail);
			void generateTHRBead(GroupOfAtoms* GroupOfBeads, int head, int tail);
			void generateCYSBead(GroupOfAtoms* GroupOfBeads, int head, int tail);
			void generatePROBead(GroupOfAtoms* GroupOfBeads, int head, int tail);

			bool isBackboneAtom(int a);
			void getMassCenter(int* AtomList, double* coord, int numOfChild);
			void getGeometryCenter(int* AtomList, double* coord, int numOfChild);
			double getInitialRadius(int* AtomList, double* center, int numOfChild);

			void getBoundingBox(GroupOfAtoms* groupOfBeads, double cutoff, int* Min, int* Max);
			void getBeadRange(int* min, int* max, double cutoff, GroupOfAtoms* groupOfBeads, IndexRange* beadRamge);
			void getResidueEnd(GroupOfAtoms* groupOfBeads, int* residueEnd);
			void setupLinearSystem(double* A, double* b, GroupOfAtoms* groupOfBeads, IndexRange* beadRange, IndexRange* atomRange);
			void setupLinearSystem2bead(double* A, double* b, GroupOfAtoms* groupOfBeads, IndexRange* beadRange, IndexRange* atomRange, double* Q);
			void setupLinearSystemAnybead2(double* A, double* b, GroupOfAtoms* groupOfBeads, IndexRange* beadRange, IndexRange* atomRange, double* Q, int* residueEnd, int* min, int* max);
			void setupLinearSystemAnybead1(double* A, double* b, GroupOfAtoms* groupOfBeads, IndexRange* beadRange, IndexRange* atomRange, double* Q, int* residueEnd);
			double getDistanceSum(int M, int* Mmin, int* Mmax, GroupOfAtoms* groupOfBeads);
			double getDistanceSum(int i, int j, int* imin, int* imax, int* jmin, int* jmax, GroupOfAtoms* groupOfBeads);
			double getyDistanceSum(int M, int* Mmin, int* Mmax, GroupOfAtoms* groupOfBeads, IndexRange* atomRange);
			double getYk(int kx, int ky, int kz, GroupOfAtoms* goa, IndexRange* atomRange);
			void LinearSystemSolver(double* AA, double* bb, int N);
			void totalChargeContained(GroupOfAtoms* groupOfBeads, double* Q, int N);
			void assignCharge(GroupOfAtoms* groupOfBeads, double* b, double* Q, int* residueEnd);
			void check(GroupOfAtoms* groupOfBeads, int* min, int* max, IndexRange* beadRange, IndexRange* atomRange);
			bool isInAtomRange(int i, int j, int k, int l, IndexRange* atomRange);
			int getResidueIndex(int beadIndex, int* residueEnd, int Nr);
			double getBeadSelfEnergy(GroupOfAtoms* groupOfBeads, int index);
			double getAtomicGBEnergy();
			double getAtomicColumbEnergy();
	};
}

#endif // COARSEGRAIN_H
