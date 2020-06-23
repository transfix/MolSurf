#include <Utility/utility.h>
#include <libCG/CoarseGrain/coarseGrain.h>

extern "C" void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);

using namespace MOLECULE;

void CoarseGrain::LinearSystemSolver(double* AA, double* bb, int N)
{
	int n = N;
	int k = 1;
	double A[n*n], b[n];
	int ipiv[n], info;
	int i, j;
	for(i = 0; i < n; i++)
	{
		for(j = i+1; j < n; j++)
		{
			A[i* n+j] = AA[i*n+j];
			A[j* n+i] = A[i*n+j];
		}
		A[i* n+i] = AA[i*n+i];
		b[i] = bb[i];
	}
	dgesv_(&n, &k, A, &n, ipiv, b, &n, &info);
	for(i = 0; i < n; i++)
	{
		bb[i] = b[i];
	}
}
