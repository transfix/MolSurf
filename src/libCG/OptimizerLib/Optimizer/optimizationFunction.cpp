#include <libCG/OptimizerLib/Optimizer/optimizationFunction.h>
#include <libCG/OptimizerLib/Optimizer/optimizer.h>

void radiusOptimizeFunction(double* p, double* x, int m, int n, void* data)
{
	FunctionData* fdata = (FunctionData*)data;
	Optimizer* opt = new Optimizer();
	opt->adata[0].center[0] = fdata->adata[0].center[0];
	opt->adata[0].center[1] = fdata->adata[0].center[1];
	opt->adata[0].center[2] = fdata->adata[0].center[2];
	opt->adata[0].radius = fdata->adata[0].radius;
	opt->samplePoints = new double*[fdata->numOfSamples];
	opt->blurValue = new double[fdata->numOfSamples];
	for(int i = 0; i < fdata->numOfSamples; i++)
	{
		opt->samplePoints[i] = new double[3];
		opt->samplePoints[i][0] = fdata->samplePoints[i][0];
		opt->samplePoints[i][1] = fdata->samplePoints[i][1];
		opt->samplePoints[i][2] = fdata->samplePoints[i][2];
		opt->blurValue[i] = fdata->blurValue[i];
	}
	opt->m_blobby = fdata->m_blobby;
	opt->numOfSamples = fdata->numOfSamples;
	opt->radiusOptimize(p, x, m, n, (void*)opt->adata);
	for(int i = 0; i < fdata->numOfSamples; i++)
	{
		delete opt->samplePoints[i];
	}
	delete opt->samplePoints;
	delete opt->blurValue;
	delete opt;
}


void centerOptimizeFunction(double* p, double* x, int m, int n, void* data)
{
	FunctionData* fdata = (FunctionData*)data;
	Optimizer* opt = new Optimizer();
	opt->adata[0].center[0] = fdata->adata[0].center[0];
	opt->adata[0].center[1] = fdata->adata[0].center[1];
	opt->adata[0].center[2] = fdata->adata[0].center[2];
	opt->adata[0].radius = fdata->adata[0].radius;
	opt->samplePoints = new double*[fdata->numOfSamples];
	opt->blurValue = new double[fdata->numOfSamples];
	for(int i = 0; i < fdata->numOfSamples; i++)
	{
		opt->samplePoints[i] = new double[3];
		opt->samplePoints[i][0] = fdata->samplePoints[i][0];
		opt->samplePoints[i][1] = fdata->samplePoints[i][1];
		opt->samplePoints[i][2] = fdata->samplePoints[i][2];
		opt->blurValue[i] = fdata->blurValue[i];
	}
	opt->m_blobby = fdata->m_blobby;
	opt->numOfSamples = fdata->numOfSamples;
	opt->centerOptimize(p, x, m, n, (void*)opt->adata);
	for(int i = 0; i < fdata->numOfSamples; i++)
	{
		delete opt->samplePoints[i];
	}
	delete opt->samplePoints;
	delete opt->blurValue;
	delete opt;
}

void radiuscenterOptimizeFunction(double* p, double* x, int m, int n, void* data)
{
	FunctionData* fdata = (FunctionData*)data;
	Optimizer* opt = new Optimizer();
	opt->adata[0].center[0] = fdata->adata[0].center[0];
	opt->adata[0].center[1] = fdata->adata[0].center[1];
	opt->adata[0].center[2] = fdata->adata[0].center[2];
	opt->adata[0].radius = fdata->adata[0].radius;
	opt->samplePoints = new double*[fdata->numOfSamples];
	opt->blurValue = new double[fdata->numOfSamples];
	for(int i = 0; i < fdata->numOfSamples; i++)
	{
		opt->samplePoints[i] = new double[3];
		opt->samplePoints[i][0] = fdata->samplePoints[i][0];
		opt->samplePoints[i][1] = fdata->samplePoints[i][1];
		opt->samplePoints[i][2] = fdata->samplePoints[i][2];
		opt->blurValue[i] = fdata->blurValue[i];
	}
	opt->m_blobby = fdata->m_blobby;
	opt->numOfSamples = fdata->numOfSamples;
	opt->radiuscenterOptimize(p, x, m, n, (void*)opt->adata);
	for(int i = 0; i < fdata->numOfSamples; i++)
	{
		delete opt->samplePoints[i];
	}
	delete opt->samplePoints;
	delete opt->blurValue;
	delete opt;
}




/*
void radiusOptimizeFunction(double *p, double *x, int m, int n, void *data)
{

  register int i,j;
  double tmp = 0.0;
  FunctionData* dptr = (FunctionData *)data;
  for(j=0; j<n; j++)
  {
    x[j] = 0.0;
    for(i=0; i<m; i++)
    {
      x[j] += evalDensity(m_blobby, dptr[i].center, samplePoints[j], p[i]);// - blurValue[j];
    //x[i]=100.0*(p[1]-p[0]*p[0])*(p[1]-p[0]*p[0]) + (p[0]-1.0)*(p[0]-1.0);
    }
    //printf("x[%d] = %f\n",j,x[j]);
    x[j] -= blurValue[j];
  }
  //printf("x[71]=%f\n",x[71]);
}

void centerOptimizeFunction(double *c, double *x, int mm, int n, void *data)
{
  register int i,j;
  double tmp = 0.0;
  FunctionData* dptr = (FunctionData *)data;
  double center[3];
  int m = mm/3;
  for(j=0; j<n; j++)
  {
    x[j] = 0.0;
    for(i=0; i<m; i++)
    {
      center[0] = c[3*i];	center[1] = c[3*i+1];	center[2] = c[3*i+2];
      x[j] += evalDensity(m_blobby, center, samplePoints[j], dptr[i].radius);// - blurValue[j];
    }
    x[j] -= blurValue[j];
  }
}


*/
