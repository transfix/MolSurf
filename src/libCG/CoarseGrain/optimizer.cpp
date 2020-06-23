#include <Utility/utility.h>
#include <float.h>
#include <libCG/CoarseGrain/optimizer.h>
#include <libCG/CoarseGrain/optimizationFunction.h>
#include <libCG/CoarseGrain/groupOfAtoms.h>
#include <libCG/CoarseGrain/atom2.h>
#include <libCG/CoarseGrain/geometry.h>
#include <libCG/CoarseGrain/surface.h>
#include "../src/levmar-2.5/levmar.h"
#include <time.h>

#ifdef _WIN32
#include <sys/types.h>
#include <sys/timeb.h>
#else
#include <sys/time.h>
#endif

#define SUPOPT_THRESHOLD 10.0
using namespace MOLECULE;
using namespace GEOMETRY;
using namespace CGOptim;

Optimizer::Optimizer(std::vector<PDBParser::Atom*> al, PDBParser::Atom* bd, Surface* surface)
{
	vector<PDBParser::Atom*>::iterator it;
	for(it = al.begin(); it != al.end(); it++)
	{
		atomList.push_back(*it);
	}
	bead = bd;
	getSamplePoints(surface);
	numOfSamples = surface->numbpts;
	blurValue = new double[numOfSamples];
	adata[0].center[0] = bead->m_Position[0];
	adata[0].center[1] = bead->m_Position[1];
	adata[0].center[2] = bead->m_Position[2];
	adata[0].radius = bead->getRadius();
	m_blobby = -2.3442;
}

Optimizer::Optimizer(double* e, double* r, double* q, int n, double q_total)
{
	numOfUnknown = n;
	E = new double[n];
	R = new double[n];
	Q = new double[n];
	for(int i = 0; i < n; i++)
	{
		E[i] = e[i];
		R[i] = r[i];
		Q[i] = q[i]; // initial value
	}
	Q_total = q_total;
}

Optimizer::Optimizer(double ae, double* r, double* q, double* x, double* y, double* z, int n, double q_total)
{
	numOfUnknown = n;
	AE = ae;
	R = new double[n];
	Q = new double[n];
	X = new double[n];
	Y = new double[n];
	Z = new double[n];
	for(int i = 0; i < n; i++)
	{
		R[i] = r[i];
		Q[i] = q[i]; // initial value
		X[i] = x[i];
		Y[i] = y[i];
		Z[i] = z[i];
	}
	Q_total = q_total;
}

double Optimizer::getMyTime()
{
#ifdef _WIN32
	time_t ltime;
	time(&ltime);
	return (double) ltime;
#else
	struct timeval t;
	gettimeofday(&t, NULL);
	return (double)(t.tv_sec + 1e-6*t.tv_usec);
#endif
}

void Optimizer::getSamplePoints(Surface* surface)
{
	int i, j, N;
	N = surface->numbpts;
	samplePoints = new double*[N];
	for(i = 0; i < N; i++)
	{
		samplePoints[i] = new double[3];
	}
	for(i = 0; i < N; i++)
		for(j = 0; j < 3; j++)
		{
			samplePoints[i][j] = surface->m_Vertices[i]->Point[j];
		}
}

double Optimizer::evalDensity(double blobby, double ctr[3], double pnt[3], double radius)
{
	double r = (ctr[0]-pnt[0])*(ctr[0]-pnt[0]) +
			   (ctr[1]-pnt[1])*(ctr[1]-pnt[1]) +
			   (ctr[2]-pnt[2])*(ctr[2]-pnt[2]);
	double r0 = radius * radius;
	double expval = blobby*r/r0 - blobby;
	return exp(expval);
}

void Optimizer::blurAtomList(std::vector<PDBParser::Atom*> atomList, double* blurValue)
{
	double x[3];
	int i, j;
	for(i = 0; i < numOfSamples; i++)
	{
		blurValue[i] = 0.0;
		for(j = 0; j < atomList.size(); j++)
		{
			PDBParser::Atom* atom = atomList[j];
			if(!atom)
			{
				continue;
			}
			// center and radii in angstroms
			x[0] = atomList[j]->m_Position[0];
			x[1] = atomList[j]->m_Position[1];
			x[2] = atomList[j]->m_Position[2];
			double r = atomList[j]->getRadius();
			blurValue[i] += evalDensity(m_blobby, x, samplePoints[i], r);
		}
	}
}

void Optimizer::initializeFunctionData(FunctionData* fdata)
{
	fdata->adata[0].center[0] = adata[0].center[0];
	fdata->adata[0].center[1] = adata[0].center[1];
	fdata->adata[0].center[2] = adata[0].center[2];
	fdata->adata[0].radius = adata[0].radius;
	fdata->samplePoints = new double*[numOfSamples];
	fdata->blurValue = new double[numOfSamples];
	for(int i = 0; i < numOfSamples; i++)
	{
		fdata->samplePoints[i] = new double[3];
		fdata->samplePoints[i][0] = samplePoints[i][0];
		fdata->samplePoints[i][1] = samplePoints[i][1];
		fdata->samplePoints[i][2] = samplePoints[i][2];
		fdata->blurValue[i] = blurValue[i];
	}
	fdata->m_blobby = m_blobby;
}

void Optimizer::initializeChargeData(FunctionData* fdata)
{
	fdata->E = new double[numOfUnknown];
	fdata->Q = new double[numOfUnknown];
	fdata->R = new double[numOfUnknown];
	for(int i = 0; i < numOfUnknown; i++)
	{
		fdata->E[i] = E[i];
		fdata->R[i] = R[i];
		fdata->Q[i] = Q[i];
	}
	fdata->Q_total = Q_total;
	fdata->numOfUnknown = numOfUnknown;
}

void Optimizer::initializeChargeData2(FunctionData* fdata)
{
	fdata->AE = AE;
	fdata->Q = new double[numOfUnknown];
	fdata->R = new double[numOfUnknown];
	fdata->X = new double[numOfUnknown];
	fdata->Y = new double[numOfUnknown];
	fdata->Z = new double[numOfUnknown];
	for(int i = 0; i < numOfUnknown; i++)
	{
		fdata->R[i] = R[i];
		fdata->Q[i] = Q[i];
		fdata->X[i] = X[i];
		fdata->Y[i] = Y[i];
		fdata->Z[i] = Z[i];
	}
	fdata->Q_total = Q_total;
	fdata->numOfUnknown = numOfUnknown;
}

void Optimizer::updateFunctionData(FunctionData* fdata, double* c, double p)
{
	fdata->adata[0].center[0] = c[0];
	fdata->adata[0].center[1] = c[1];
	fdata->adata[0].center[2] = c[2];
	fdata->adata[0].radius = p;
}
void Optimizer::updateChargeData(FunctionData* fdata, double* q)
{
	for(int i = 0; i < numOfUnknown; i++)
	{
		fdata->Q[i] = q[i];
	}
}

void Optimizer::radiusOptimize(double* p, double* x, int m, int n, void* data)
{
	register int i,j,k;
	atomData* dptr = (atomData*)data;
	/* 1D least square version (fast) */
	for(k = 0; k < n; k++)
	{
		x[k] = 0.0;
		for(j = 0; j < numOfSamples; j++)
		{
			double tmp = 0.0;
			for(i=0; i<m; i++)
			{
				tmp += evalDensity(m_blobby, dptr[i].center, samplePoints[j], p[i]);    // - blurValue[j];
			}
			x[k] += (tmp-blurValue[j])*(tmp-blurValue[j]);
		}
	}
}

void Optimizer::centerOptimize(double* c, double* x, int mm, int n, void* data)
{
	register int i,j,k;
	atomData* dptr = (atomData*)data;
	double center[3];
	int m = mm/3;
	/* 1D least square version (fast)*/
	for(k = 0; k < n; k++)
	{
		x[k] = 0.0;
		for(j = 0; j < numOfSamples; j++)
		{
			double tmp = 0.0;
			for(i=0; i<m; i++)
			{
				center[0] = c[3*i];
				center[1] = c[3*i+1];
				center[2] = c[3*i+2];
				tmp += evalDensity(m_blobby, center, samplePoints[j], dptr[i].radius);// - blurValue[j];
			}
			x[k] += (tmp-blurValue[j])*(tmp-blurValue[j]);
		}
	}
}

void Optimizer::radiuscenterOptimize(double* c, double* x, int m, int n, void* data)
{
	register int i,j,k;
	//  atomData* dptr = (atomData *)data;
	double center[3];
	double radius;
	/* 1D least square version (fast)*/
	for(k = 0; k < n; k++)
	{
		x[k] = 0.0;
		for(j = 0; j < numOfSamples; j++)
		{
			double tmp = 0.0;
			for(i=0; i<m; i++)
			{
				center[0] = c[3*i];
				center[1] = c[3*i+1];
				center[2] = c[3*i+2];
				radius = c[3*i+3];
				tmp += evalDensity(m_blobby, center, samplePoints[j], radius);// - blurValue[j];
			}
			x[k] += (tmp-blurValue[j])*(tmp-blurValue[j]);
		}
	}
}

void Optimizer::chargeOptimize(double* q, double* x, int m, int n, void* data)
{
	register int i,j;
	double tmp = 0.0;
	for(j = 0; j < n; j++)
	{
		x[j] = 0.0;
		for(i = 0; i < m; i++)
		{
			tmp = (E[i]- q[i]*q[i]/R[i]);
			x[j] += tmp*tmp;
		}
	}
}

void Optimizer::chargeOptimize2(double* q, double* x, int m, int n, void* data)
{
	register int i,j,k;
	double rik2;
	double tmp = 0.0;
	for(j = 0; j < n; j++)
	{
		x[j] = 0.0;
		for(i = 0; i < m; i++)
		{
			for(k = 0; k < m; k++)
			{
				rik2 = (X[i]-X[k])*(X[i]-X[k]) + (Y[i]-Y[k])*(Y[i]-Y[k]) + (Z[i]-Z[k])*(Z[i]-Z[k]);
				//	if(rik2 > SUPOPT_THRESHOLD)
				//		continue;
				tmp = q[i]*q[k]/sqrt(rik2 + R[i]*R[k] * exp(-1.0*rik2 / (4.0*R[i]*R[k])));
				//		tmp = q[i]*q[k]/sqrt(rik2);
				x[j] += tmp;
			}
		}
		x[j] -= AE;
	}
}

void Optimizer::optimizeCenterAndRadius(int iteration)
{
	// m - num of CG atoms
	// n - num of surface pts
	// p - CG radii
	// x - f(CG) - f(AA), the chi function
	double p[1], c[3], cp[4];
	int m, n;
	register int i, iter;
	double start, end;
	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	int ret;
	opts[0]=LM_INIT_MU;
	opts[1]=1E-15;
	opts[2]=1E-15;
	opts[3]=1E-20;
	opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference jacobian version is used
	//	std::vector<PDBParser::Atom*> atomList = getAtomList();
	blurAtomList(atomList,blurValue);
	n = numOfSamples;
	double x[n]; // chi function
	srand48((long)time(NULL));
	p[0] = adata[0].radius;
	c[0] = adata[0].center[0];
	c[1] = adata[0].center[1];
	c[2] = adata[0].center[2];
	cp[0] = c[0];
	cp[1] = c[1];
	cp[2] = c[2];
	cp[3] = p[0];
	for(i = 0; i < n; i++)
	{
		x[i] = 0.0;
	}
	FunctionData* fdata = new FunctionData;
	initializeFunctionData(fdata);
	start = getMyTime();
	for(iter = 0; iter < iteration; iter++)
	{
		// optimize radius first since the initial centers are close to the optimal centers
		m = 1;
		ret=dlevmar_dif(radiuscenterOptimizeFunction, p, x, m, n, 1000, opts, info, NULL, NULL, (void*)fdata);
		updateFunctionData(fdata, c, p[0]);
		int mm = 3*m;
		ret = dlevmar_dif(radiuscenterOptimizeFunction, c, x, mm, n, 1000, opts, info, NULL, NULL, (void*)fdata);
		updateFunctionData(fdata, c, p[0]);
		/*
				m = 1;
				int mmm = 4*m;
				ret = dlevmar_dif(radiuscenterOptimizeFunction, cp, x, mmm, n, 1000, opts, info, NULL, NULL, (void*)fdata);
				updateFunctionData(fdata, cp, cp[3]);
		*/
	}
	end = getMyTime();
	adata[0].center[0] = fdata->adata[0].center[0];
	adata[0].center[1] = fdata->adata[0].center[1];
	adata[0].center[2] = fdata->adata[0].center[2];
	adata[0].radius = fdata->adata[0].radius;
	//  ret = dlevmar_dif(radiuscenterOptimizeFunction, c, x, m, n, 50, opts, info, NULL, NULL, (void*)adata);
	// no jacobian
	/*
	  printf("Levenberg-Marquardt returned %d in %g iter, takes %f seconds, reason %g\nSolution: \n", ret, info[5], end-start, info[6]);
	  for(i=0; i<m; ++i)
	  {
	      if (p[i] < 0.0)	p[i] = -p[i];
	      printf("%.7g %.7g %.7g %.7g\n", adata[i].center[0], adata[i].center[1], adata[i].center[2], adata[i].radius);
	  }
	  printf("\nMinimization info:\n");
	  for(i=0; i<LM_INFO_SZ; ++i)
	    printf("%g ", info[i]);
	  printf("\n");
	*/
	/*
	 double info[LM_INFO_SZ],
	 O: information regarding the minimization. Set to NULL if don't care
	 info[0]= ||e||_2 at initial p.
	 info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||Dp||_2, \mu/max[J^T J]_ii ], all computed at estimated p.
	 info[5]= # iterations,
	 info[6]=reason for terminating: 1 - stopped by small gradient J^T e
	 2 - stopped by small Dp
	 3 - stopped by itmax
	 4 - singular matrix. Restart from current p with increased \mu
	 5 - no further error reduction is possible. Restart with increased mu
	 6 - stopped by small ||e||_2
	 7 - stopped by invalid (i.e. NaN or Inf) "func" values. This is a user error
	 info[7]= # function evaluations
	 info[8]= # Jacobian evaluations
	 */
}


void Optimizer::optimizeCharge()
{
	// m - num of CG atoms
	// n - num of surface pts
	// p - CG radii
	// x - f(CG) - f(AA), the chi function
	double q[numOfUnknown];
	int m, n;
	register int i;
	double start, end;
	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	int ret;
	opts[0]=LM_INIT_MU;
	opts[1]=1E-15;
	opts[2]=1E-15;
	opts[3]=1E-20;
	opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference jacobian version is used
	m = numOfUnknown;
	n = m;
	double x[n]; // chi function
	srand48((long)time(NULL));
	for(i = 0; i < numOfUnknown; i++)
	{
		q[i] = Q[i];
	}
	for(i = 0; i < n; i++)
	{
		x[i] = 0.0;
	}
	FunctionData* fdata = new FunctionData;
	//	initializeChargeData(fdata);  // optimize self energy only
	initializeChargeData2(fdata); // optimize GB energy
	start = getMyTime();
	// constrained condition
	double A[numOfUnknown], b[1];
	for(i = 0; i < numOfUnknown; i++)
	{
		A[i] = 1.0;
	}
	b[0] = Q_total;
	printf("%d %f \n", numOfUnknown, Q_total);
	/*
	      double A[3*5]={1.0, 3.0, 0.0, 0.0, 0.0,  0.0, 0.0, 1.0, 1.0, -2.0,  0.0, 1.0, 0.0, 0.0, -1.0},
	             b[3]={0.0, 0.0, 0.0};
	*/
	//	ret = dlevmar_lec_dif(chargeOptimizeFunction, q, x, m, n, A, b, 1, 1000, opts, info, NULL, NULL, (void*)fdata); // optimize self energy only
	//	ret = dlevmar_lec_dif(chargeOptimizeFunction2, q, x, m, n, A, b, 1, 1000, opts, info, NULL, NULL, (void*)fdata); // optimize GB energy
	ret = dlevmar_dif(chargeOptimizeFunction2, q, x, m, n, 6/*000*/, opts, info, NULL, NULL, (void*)fdata);
	updateChargeData(fdata, q);
	end = getMyTime();
	for(i = 0; i < numOfUnknown; i++)
	{
		Q[i] = fdata->Q[i];
	}
	{
		int j,k;
		double rik2;
		double tmp = 0.0;
		for(j = 0; j < m; j++)
		{
			for(k = 0; k < m; k++)
			{
				rik2 = (X[j]-X[k])*(X[j]-X[k]) + (Y[j]-Y[k])*(Y[j]-Y[k]) + (Z[j]-Z[k])*(Z[j]-Z[k]);
				tmp += Q[j]*Q[k]/sqrt(rik2 + R[j]*R[k] * exp(-1.0*rik2 / (4.0*R[j]*R[k])));
			}
		}
		printf("check: %f\n", tmp);
	}
	printf("Levenberg-Marquardt returned %d in %g iter, takes %f seconds, reason %g\nSolution: \n", ret, info[5], end-start, info[6]);
	/*
	  for(i=0; i<m; ++i)
	  {
	      if (p[i] < 0.0)   p[i] = -p[i];
	      printf("%.7g %.7g %.7g %.7g\n", adata[i].center[0], adata[i].center[1], adata[i].center[2], adata[i].radius);
	  }
	*/
	printf("\nMinimization info:\n");
	for(i=0; i<LM_INFO_SZ; ++i)
	{
		printf("%g ", info[i]);
	}
	printf("\n");
}
