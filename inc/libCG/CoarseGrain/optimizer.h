#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include <libCG/CoarseGrain/optimizationFunction.h>
#include <PDBParser/Atom.h>
#include <PDBParser/GroupOfAtoms.h>
#include <libCG/CoarseGrain/surface.h>
#include <vector>

namespace CGOptim {

class Optimizer {
public:
  Optimizer() {}
  Optimizer(std::vector<PDBParser::Atom *> al, PDBParser::Atom *bd,
            GEOMETRY::Surface *surface);
  Optimizer(double *e, double *r, double *q, int n, double q_total);
  Optimizer(double ae, double *r, double *q, double *x, double *y, double *z,
            int n, double q_total);
  void radiusOptimize(double *p, double *x, int m, int n, void *data);
  void centerOptimize(double *c, double *x, int mm, int n, void *data);
  void radiuscenterOptimize(double *c, double *x, int m, int n, void *data);
  void chargeOptimize(double *r, double *x, int m, int n, void *data);
  void chargeOptimize2(double *r, double *x, int m, int n, void *data);
  void optimizeCenterAndRadius(int iteration);
  void optimizeCharge();

  // radius and center optimization
  atomData adata[1];
  double **samplePoints;
  int numOfSamples;
  double *blurValue;
  double m_blobby;

  // charge optimization
  int numOfUnknown;
  double *E;
  double *R;
  double Q_total;
  double *Q;
  double AE;
  double *X, *Y, *Z;

private:
  double getMyTime();
  std::vector<PDBParser::Atom *> getAtomList();
  void getSamplePoints(GEOMETRY::Surface *surface);
  double evalDensity(double blobby, double ctr[3], double pnt[3],
                     double radius);
  void blurAtomList(std::vector<PDBParser::Atom *> atomList, double *blurValue);
  void initializeFunctionData(FunctionData *fdata);
  void updateFunctionData(FunctionData *fdata, double *c, double p);
  void initializeChargeData(FunctionData *fdata);
  void initializeChargeData2(FunctionData *fdata);
  void updateChargeData(FunctionData *fdata, double *q);

  std::vector<PDBParser::Atom *> atomList;
  PDBParser::Atom *bead;
};

} // namespace CGOptim
#endif // OPTIMIZER_H
