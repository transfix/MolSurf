#ifndef OPTIMIZATIONFUNCTION_H
#define OPTIMIZATIONFUNCTION_H

typedef struct atomDataType {
  double radius;
  double center[3];
} atomData;

typedef struct DataType {
  atomData adata[1];
  double **samplePoints;
  double *blurValue;
  double m_blobby;
  int numOfSamples;
} FunctionData;

void radiusOptimizeFunction(double *p, double *x, int m, int n, void *data);
void centerOptimizeFunction(double *c, double *x, int mm, int n, void *data);
void radiuscenterOptimizeFunction(double *c, double *x, int mm, int n,
                                  void *data);

#endif // OPITMIZATIONFUNCTION_H
