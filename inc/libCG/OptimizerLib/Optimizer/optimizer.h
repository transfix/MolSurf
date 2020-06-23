#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include <libCG/CoarseGrain/surface.h>
#include <libCG/OptimizerLib/Optimizer/optimizationFunction.h>
#include <PDBParser/Atom.h>
#include <PDBParser/GroupOfAtoms.h>
#include <vector>

class Optimizer
{
	public:
		Optimizer() {}
		Optimizer(std::vector<PDBParser::Atom*>& al, PDBParser::Atom* bd, GEOMETRY::Surface* surface);
		~Optimizer();
		void radiusOptimize(double* p, double* x, int m, int n, void* data);
		void centerOptimize(double* c, double* x, int mm, int n, void* data);
		void radiuscenterOptimize(double* c, double* x, int m, int n, void* data);
		void optimizeCenterAndRadius(int iteration);
		void optimizeCenterAndRadiusN(int iteration);
		atomData adata[1];
		double** samplePoints;
		int numOfSamples;
		double* blurValue;
		double m_blobby;
	private:
		double getMyTime();
		void getSamplePoints(GEOMETRY::Surface* surface);
		double evalDensity(double blobby, double ctr[3], double pnt[3], double radius);
		void blurAtomList(std::vector<PDBParser::Atom*> atomList, double* blurValue);
		void initializeFunctionData(FunctionData* fdata);
		void updateFunctionData(FunctionData* fdata, double* c, double p);

		void updateFunctionDataN(FunctionData* fdata, double* p);
		std::vector<PDBParser::Atom*> atomList;
		PDBParser::Atom* bead;
};
#endif // OPTIMIZER_H
